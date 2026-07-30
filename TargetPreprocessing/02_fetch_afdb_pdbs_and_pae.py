"""Download AlphaFold DB structure and confidence files for UniProt IDs.

This script:
1) Reads UniProt accessions from a CSV column.
2) Tries default AFDB IDs and isoforms (-1 to -10) to find available models.
3) Downloads matching PDB and PAE JSON files into AFDB_<version>/PDB and AFDB_<version>/PAE.

Default input: Data/UniProt/uniprot_human_reviewed.csv
Default output root: Data/AFDB_<version>/
"""

import argparse
import csv
import os
from pathlib import Path

import requests
from tqdm import tqdm

AFDB_VERSION = "v6"
DEFAULT_OUTPUT_ROOT = Path(__file__).resolve().parent.parent / "Data"
DEFAULT_TIMEOUT = 15

def _parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent.parent
    default_csv = script_dir / "Data" / "UniProt" / "uniprot_reviewed_annotated.csv"
    default_output_dir = DEFAULT_OUTPUT_ROOT
    parser = argparse.ArgumentParser(description="Download AlphaFold PDB and PAE by UniProt IDs from CSV.")
    parser.add_argument("--csv", default=str(default_csv), help="Path to CSV containing UniProt IDs.")
    parser.add_argument("--id-column", default="uniprot_id", help="CSV column name for UniProt IDs.")
    parser.add_argument("--output-dir", default=str(default_output_dir), help="Base directory for outputs (AFDB_<version> will be created).")
    parser.add_argument("--afdb-version", default=AFDB_VERSION, help="AlphaFold DB version (e.g., v6).")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT, help="HTTP timeout in seconds.")
    return parser.parse_args()

def _read_uniprot_ids(csv_path: Path, column: str) -> list[str]:
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError("CSV has no header.")
        col_map = {c.lower(): c for c in reader.fieldnames}
        if column.lower() not in col_map:
            raise ValueError(f"Column '{column}' not found in CSV.")
        key = col_map[column.lower()]
        return [row[key].strip() for row in reader if row.get(key, "").strip()]

def _resolve_output_dir(base_dir: Path | None, afdb_version: str) -> Path:
    root = base_dir if base_dir else DEFAULT_OUTPUT_ROOT
    out_dir = root / f"AFDB_{afdb_version}"
    (out_dir / "PDB").mkdir(parents=True, exist_ok=True)
    (out_dir / "PAE").mkdir(parents=True, exist_ok=True)
    return out_dir

def _candidate_af_ids(uniprot_id: str) -> list[str]:
    return [uniprot_id] + [f"{uniprot_id}-{i}" for i in range(1, 11)]

def _download(url: str, save_path: Path, label: str, uniprot_id: str, timeout: int) -> bool:
    try:
        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
    except requests.RequestException as e:
        #print(f"[!] {label} failed for {uniprot_id}: {e}")
        return False
    save_path.write_text(resp.text)
    return True

def download_pdb(uniprot_id: str, out_dir: Path, afdb_version: str, timeout: int) -> tuple[bool, str | None]:
    for af_id in _candidate_af_ids(uniprot_id):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{af_id}-F1-model_{afdb_version}.pdb"
        if _download(url, out_dir / "PDB" / f"{uniprot_id}.pdb", f"PDB ({af_id})", uniprot_id, timeout):
            return True, af_id
    print(f"No AFDB PDB found for {uniprot_id} (default or isoforms 1-10).")
    return False, None

def download_pae(uniprot_id: str, out_dir: Path, af_id: str | None, afdb_version: str, timeout: int) -> bool:
    if not af_id:
        print(f"[!] Skipping PAE for {uniprot_id}: no available pdb.")
        return False
    url = f"https://alphafold.ebi.ac.uk/files/AF-{af_id}-F1-predicted_aligned_error_{afdb_version}.json"
    return _download(url, out_dir / "PAE" / f"{uniprot_id}.pae.json", f"PAE ({af_id})", uniprot_id, timeout)

if __name__ == "__main__":
    args = _parse_args()
    ids = _read_uniprot_ids(Path(args.csv), args.id_column)
    out_dir = _resolve_output_dir(Path(args.output_dir) if args.output_dir else None, args.afdb_version)

    for uniprot_id in tqdm(ids, desc="Downloading PDB/PAE"):
        pdb_ok, af_id = download_pdb(uniprot_id, out_dir, args.afdb_version, args.timeout)
        pae_ok = download_pae(uniprot_id, out_dir, af_id, args.afdb_version, args.timeout)
        if not pdb_ok or not pae_ok:
            print(f"[!] Incomplete download for {uniprot_id}: PDB={pdb_ok}, PAE={pae_ok}")
    print("Done.")