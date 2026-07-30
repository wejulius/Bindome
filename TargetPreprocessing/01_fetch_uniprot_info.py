import csv
import os
import requests
from pathlib import Path
import argparse
import json

import tqdm
from classify_uniprot_entry import classify_uniprot_id, get_subcellular_localizations

DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent / "Data" / "UniProt"
DEFAULT_ORGANISM_ID = 9606 # Human

def _reviewed_query_clause(reviewed: bool | None) -> str:
    if reviewed is True:
        return "reviewed:true"
    if reviewed is False:
        return ""
    return ""


def fetch_uniprot_nonmembrane(out_csv: Path, organism_id: int = DEFAULT_ORGANISM_ID, reviewed: bool | None = True) -> None:
    url = "https://rest.uniprot.org/uniprotkb/stream"
    query_parts = []
    clause = _reviewed_query_clause(reviewed)
    if clause:
        query_parts.append(clause)
    query_parts.append(f"organism_id:{organism_id}")
    query_parts.append("(NOT keyword:KW-0472)")
    params = {
        "query": " AND ".join(query_parts),
        "format": "tsv",
        "fields": "accession,id,gene_names,protein_name",
    }
    resp = requests.get(url, params=params, headers={"Accept": "text/tab-separated-values"}, timeout=60)
    if not resp.ok:
        raise requests.HTTPError(f"{resp.status_code} {resp.reason}: {resp.text}", response=resp)

    lines = resp.text.strip().splitlines()
    rows = [line.split("\t") for line in lines[1:]]  # skip header

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["uniprot_id", "uniprot_name", "gene_names", "protein_name"])
        writer.writerows(rows)


def fetch_uniprot_nonmembrane_reviewed(out_csv: Path, organism_id: int = DEFAULT_ORGANISM_ID) -> None:
    return fetch_uniprot_nonmembrane(out_csv, organism_id=organism_id, reviewed=True)

def fetch_uniprot_entries(out_csv: Path, organism_id: int = DEFAULT_ORGANISM_ID, reviewed: bool | None = True) -> None:
    url = "https://rest.uniprot.org/uniprotkb/stream"
    query_parts = []
    clause = _reviewed_query_clause(reviewed)
    if clause:
        query_parts.append(clause)
    query_parts.append(f"organism_id:{organism_id}")
    params = {
        "query": " AND ".join(query_parts),
        "format": "tsv",
        "fields": "accession,id,gene_names,protein_name,sequence",
    }
    resp = requests.get(url, params=params, headers={"Accept": "text/tab-separated-values"}, timeout=60)
    if not resp.ok:
        raise requests.HTTPError(f"{resp.status_code} {resp.reason}: {resp.text}", response=resp)

    lines = resp.text.strip().splitlines()
    rows = [line.split("\t") for line in lines[1:]]  # skip header

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["uniprot_id", "uniprot_name", "gene_names", "protein_name", "sequence"])
        writer.writerows(rows)


def fetch_uniprot_reviewed(out_csv: Path, organism_id: int = DEFAULT_ORGANISM_ID) -> None:
    return fetch_uniprot_entries(out_csv, organism_id=organism_id, reviewed=True)

def _get_uniprot_json(accession: str, session: requests.Session) -> dict:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    r = session.get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def _safe_uniprot_feature_boundary(location: dict, boundary: str) -> int | None:
    """Return integer boundary when available, otherwise None.

    UniProt feature locations are not always complete; in rare cases a
    start/end value is missing or non-numeric.
    """
    value = (location.get(boundary) or {}).get("value")
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

def uniprot_membrane_embedded_residues(accession: str, data: dict | None = None, session: requests.Session | None = None) -> dict:
    if data is None:
        if session is None:
            session = requests.Session()
        data = _get_uniprot_json(accession, session)
    seq = data["sequence"]["value"]
    L = len(seq)

    membrane_embedded_ranges = []  # list of (start, end)
    for ft in data.get("features", []):
        ftype = ft.get("type")
        if ftype in ("Transmembrane", "Intramembrane"):
            location = ft.get("location") or {}
            start = _safe_uniprot_feature_boundary(location, "start")
            end = _safe_uniprot_feature_boundary(location, "end")
            if start is None or end is None:
                continue
            membrane_embedded_ranges.append((start, end))

    return {
        "sequence": seq,
        "length": L,
        "membrane_embedded_ranges": membrane_embedded_ranges,
    }

def uniprot_has_membrane_keyword(accession: str, data: dict | None = None, session: requests.Session | None = None) -> bool:
    if data is None:
        if session is None:
            session = requests.Session()
        data = _get_uniprot_json(accession, session)
    for kw in data.get("keywords", []):
        if kw.get("id") == "KW-0472":
            return True
    return False

def update_annotations(out_csv: Path | None = None) -> None:
    in_csv = out_csv
    with in_csv.open("r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    with requests.Session() as session:
        for row in tqdm.tqdm(rows, desc="Updating membrane annotations", total=len(rows)):
            accession = row.get("uniprot_id", "").strip()
            if not accession:
                continue
            seq_in_csv = (row.get("sequence") or "").strip()

            data = _get_uniprot_json(accession, session)

            subcellular_localizations = get_subcellular_localizations(data)
            row["subcellular_localization"] = json.dumps(subcellular_localizations)

            is_membrane = uniprot_has_membrane_keyword(accession, data=data)
            row["membrane_protein"] = is_membrane

            info = uniprot_membrane_embedded_residues(accession, data=data)
            seq_from_json = info.get("sequence", "")
            ranges = info.get("membrane_embedded_ranges", [])
            row["membrane_embedded_ranges"] = json.dumps(ranges if seq_in_csv == seq_from_json else [])

            try:
                classification = classify_uniprot_id(accession, entry=data)
                row["class_14"] = classification.class_name
                row["class_rule_id"] = classification.rule_id
            except Exception:
                row["class_14"] = "Unassigned"
                row["class_rule_id"] = "CLASSIFY_ERROR"

    fieldnames = list(rows[0].keys()) if rows else ["uniprot_id", "uniprot_name", "gene_names", "protein_name", "sequence", "membrane_embedded_ranges", "membrane_protein", "subcellular_localization", "class_14", "class_rule_id"]
    if "membrane_embedded_ranges" not in fieldnames:
        fieldnames.append("membrane_embedded_ranges")
    if "membrane_protein" not in fieldnames:
        fieldnames.append("membrane_protein")
    if "subcellular_localization" not in fieldnames:
        fieldnames.append("subcellular_localization")
    if "class_14" not in fieldnames:
        fieldnames.append("class_14")
    if "class_rule_id" not in fieldnames:
        fieldnames.append("class_rule_id")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch reviewed UniProt entries to CSV.")
    parser.add_argument("--output-dir", type=str, default=str(DEFAULT_OUTPUT_DIR), help="Directory to write the CSV file.")
    parser.add_argument("--organism-id", type=int, default=DEFAULT_ORGANISM_ID, help="Organism ID to fetch UniProt entries for.")
    parser.add_argument(
        "--reviewed",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Restrict query to reviewed UniProt entries (use --no-reviewed to include unreviewed entries).",
    )
    parser.add_argument("--print-args", action="store_true", help="Print the parsed arguments to stdout before running.")
    return parser.parse_args()

if __name__ == "__main__":
    args = _parse_args()
    if args.print_args:
        print("Parsed arguments:")
        for arg, value in vars(args).items():
            print(f"  {arg}: {value}")
    out_dir = Path(args.output_dir)
    out_path = out_dir / "uniprot_reviewed_annotated.csv"
    os.makedirs(out_path.parent, exist_ok=True)
    fetch_uniprot_entries(out_path, args.organism_id, reviewed=args.reviewed)
    update_annotations(out_path)
    print(f"Written updated UniProt reviewed entries to {out_path}")
    print("Done.")