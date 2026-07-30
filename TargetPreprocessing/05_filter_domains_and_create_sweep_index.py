"""Filter domain statistics and create sweep index TSV batches.

This script:
1) Loads per-domain statistics from `domain_stats.csv`.
2) Applies quality/structure/membrane filters based on CLI thresholds.
3) Writes the filtered stats table for traceability.
4) Builds one or more sweep index TSV files for downstream design runs.

Default outputs:
- Filtered stats: Data/Filtered/filtered_domain_stats.csv
- Sweep index batches: Data/Filtered/sweep_index_batch*.tsv
"""

import argparse
import pandas as pd
from pathlib import Path
import csv

DEFAULT_STATS_CSV_PATH = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Statistics" / "domain_stats.csv"
DEFAULT_IDP_STATS_CSV_PATH = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Statistics" / "idp_crops_stats.csv"
DEFAULT_DOMAIN_SWEEP_INDEX_PATH = Path(__file__).resolve().parent.parent / "Data" / "Filtered" / "sweep_index.tsv"
DEFAULT_IDP_SWEEP_INDEX_PATH = Path(__file__).resolve().parent.parent / "Data" / "Filtered" / "idp_sweep_index.tsv"
DEFAULT_DOMAIN_FILTERED_STATS_CSV_PATH = Path(__file__).resolve().parent.parent / "Data" / "Filtered" / "filtered_domain_stats.csv"
DEFAULT_IDP_FILTERED_STATS_CSV_PATH = Path(__file__).resolve().parent.parent / "Data" / "Filtered" / "filtered_idp_crops_stats.csv"
DEFAULT_NUMBER_OF_SWEEP_BATCHES = 4

# Filters (structured mode)
DEFAULT_MIN_PLDDT = 70.0
DEFAULT_MIN_DOMAIN_LENGTH = 30
DEFAULT_MAX_DOMAIN_LENGTH = 550
DEFAULT_MAX_RADIUS_OF_GYRATION_NORMALIZED = 5.0
DEFAULT_MIN_CONTACTS_PER_RESIDUE = 2.0
DEFAULT_DROP_ONLY_SINGLE_HELIX_DOMAINS = True
DEFAULT_DROP_ONLY_DOUBLE_HELIX_DOMAINS = True
DEFAULT_DROP_ONLY_LOOPY_DOMAINS = True
DEFAULT_DROP_UNIPORT_AFDB_SEQUENCE_MISMATCHES = True
DEFAULT_DROP_MEMBRANE_EMBEDDED_RESIDUES_OVERLAP_DOMAIN = True

# Filters (IDP mode — inverted/relaxed relative to structured mode)
DEFAULT_IDP_MAX_PLDDT = 50.0        # keep only low-confidence (disordered) regions
DEFAULT_IDP_MIN_AF_BIND_W = 0.0     # minimum mean AF-Bind weighted score
DEFAULT_IDP_MIN_DOMAIN_LENGTH = 15  # shortest sensible IDP crop
DEFAULT_IDP_MAX_DOMAIN_LENGTH = 30  # longest sensible IDP crop (window-based)

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter domains/IDP crops and create sweep index.")
    parser.add_argument("--idp", action="store_true", default=False,
                        help="IDP mode: read idp_crops_stats.csv, invert structural filters, "
                             "apply AF-Bind score filter. Defaults for --stats-csv-path, "
                             "--sweep-index-path, and --filtered-stats-csv-path change accordingly.")
    parser.add_argument("--stats-csv-path", type=str, default=None,
                        help="Path to domain stats CSV. Defaults to domain_stats.csv (structured) "
                             "or idp_crops_stats.csv (--idp).")
    parser.add_argument("--sweep-index-path", type=str, default=None,
                        help="Path to write the domain index TSV. Defaults differ by mode.")
    parser.add_argument("--filtered-stats-csv-path", type=str, default=None,
                        help="Path to write the filtered stats CSV. Defaults differ by mode.")
    parser.add_argument("--min-plddt", type=float, default=DEFAULT_MIN_PLDDT,
                        help="Minimum average pLDDT to keep a domain.")
    parser.add_argument("--max-radius-of-gyration-normalized", type=float, default=DEFAULT_MAX_RADIUS_OF_GYRATION_NORMALIZED,
                        help="Maximum normalized radius of gyration to keep a domain.")
    parser.add_argument("--min-contacts-per-residue", type=float, default=DEFAULT_MIN_CONTACTS_PER_RESIDUE,
                        help="Minimum average CA contacts per residue to keep a domain.")
    parser.add_argument("--min-domain-length", type=int, default=DEFAULT_MIN_DOMAIN_LENGTH,
                        help="Minimum domain length to keep.")
    parser.add_argument("--max-domain-length", type=int, default=DEFAULT_MAX_DOMAIN_LENGTH,
                        help="Maximum domain length to keep.")
    parser.add_argument("--drop-only-single-helix-domains", action=argparse.BooleanOptionalAction, default=DEFAULT_DROP_ONLY_SINGLE_HELIX_DOMAINS,
                        help="Exclude domains that are only single helix.")
    parser.add_argument("--drop-only-double-helix-domains", action=argparse.BooleanOptionalAction, default=DEFAULT_DROP_ONLY_DOUBLE_HELIX_DOMAINS,
                        help="Exclude domains that are only double helix.")
    parser.add_argument("--drop-only-loopy-domains", action=argparse.BooleanOptionalAction, default=DEFAULT_DROP_ONLY_LOOPY_DOMAINS,
                        help="Exclude domains that are only loopy.")
    parser.add_argument("--drop-uniprot-afdb-sequence-mismatches", action=argparse.BooleanOptionalAction, default=DEFAULT_DROP_UNIPORT_AFDB_SEQUENCE_MISMATCHES,
                        help="Exclude domains where the AFDB sequence does not match the UniProt sequence.")
    parser.add_argument("--drop-membrane-embedded-residues-overlap-domain", action=argparse.BooleanOptionalAction, default=DEFAULT_DROP_MEMBRANE_EMBEDDED_RESIDUES_OVERLAP_DOMAIN,
                        help="Exclude domains that overlap with membrane embedded residues.")
    parser.add_argument("--num-sweep-batches", type=int, default=DEFAULT_NUMBER_OF_SWEEP_BATCHES,
                        help="Number of sweep index TSV batches to write.")
    # IDP-specific filters (only active when --idp is set)
    parser.add_argument("--idp-max-plddt", type=float, default=DEFAULT_IDP_MAX_PLDDT,
                        help="[--idp] Maximum mean pLDDT to keep a crop (default: 50).")
    parser.add_argument("--idp-min-af-bind-w", type=float, default=DEFAULT_IDP_MIN_AF_BIND_W,
                        help="[--idp] Minimum mean AF-Bind weighted score to keep a crop.")
    return parser.parse_args()


def write_domain_index(sweep_index_path, config, rows):
    """
    Write a TSV index with a config summary header.
    rows: list of dicts with keys: uniprot, domain_name, domain_pdb_path, full_pdb_path, json_path, domain_length, design_time, design_runs, num_accepted_designs, num_trajectories
    """
    sweep_index_path = Path(sweep_index_path)
    sweep_index_path.parent.mkdir(parents=True, exist_ok=True)
    header_lines = ["# Domain segmentation config:"]
    for k, v in config.items():
        header_lines.append(f"# {k}: {v}")
    header_lines.append("# columns: uniprot domain_name domain_pdb_path full_pdb_path json_path domain_length design_time design_runs num_accepted_designs num_trajectories")
    with sweep_index_path.open("w", newline="") as handle:
        handle.write("\n".join(header_lines) + "\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for r in rows:
            writer.writerow([
                r["uniprot"],
                r["domain_name"],
                r["domain_pdb_path"],
                r["full_pdb_path"],
                r["json_path"],
                r["domain_length"],
                0,
                0,
                0,
                0,
            ])


def _split_rows(rows, num_batches):
    if num_batches < 1:
        raise ValueError("num_sweep_batches must be >= 1")
    batch_size = (len(rows) + num_batches - 1) // num_batches or 1
    for batch_idx in range(num_batches):
        start = batch_idx * batch_size
        end = start + batch_size
        if start >= len(rows):
            break
        yield rows[start:end]

def write_domain_index_batches(sweep_index_path, config, rows, num_batches):
    sweep_index_path = Path(sweep_index_path)
    for idx, batch_rows in enumerate(_split_rows(rows, num_batches), start=1):
        batch_path = sweep_index_path.with_name(f"{sweep_index_path.stem}_batch{idx}{sweep_index_path.suffix}")
        batch_config = {**config, "batch": f"{idx}/{num_batches}"}
        write_domain_index(batch_path, batch_config, batch_rows)

def main() -> None:
    args = _parse_args()

    IDP_MODE = args.idp

    # Resolve mode-dependent path defaults
    stats_csv_path = Path(
        args.stats_csv_path if args.stats_csv_path
        else (DEFAULT_IDP_STATS_CSV_PATH if IDP_MODE else DEFAULT_STATS_CSV_PATH)
    )
    sweep_index_path = Path(
        args.sweep_index_path if args.sweep_index_path
        else (DEFAULT_IDP_SWEEP_INDEX_PATH if IDP_MODE else DEFAULT_DOMAIN_SWEEP_INDEX_PATH)
    )
    filtered_stats_csv_path = Path(
        args.filtered_stats_csv_path if args.filtered_stats_csv_path
        else (DEFAULT_IDP_FILTERED_STATS_CSV_PATH if IDP_MODE else DEFAULT_DOMAIN_FILTERED_STATS_CSV_PATH)
    )

    print("Parsed arguments:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")
    print(f"  [resolved] stats_csv_path: {stats_csv_path}")
    print(f"  [resolved] sweep_index_path: {sweep_index_path}")
    print(f"  [resolved] filtered_stats_csv_path: {filtered_stats_csv_path}")
    print(f"  mode: {'IDP' if IDP_MODE else 'structured'}")

    NUM_SWEEP_BATCHES = args.num_sweep_batches

    df = pd.read_csv(stats_csv_path)
    print(f"\nLoaded {len(df)} rows from {stats_csv_path}")

    def _apply_filter(df: pd.DataFrame, mask, label: str) -> pd.DataFrame:
        before = len(df)
        df = df[mask]
        after = len(df)
        removed = before - after
        print(f"  {label}: {before} -> {after} (removed {removed})")
        return df

    if IDP_MODE:
        # ── IDP mode: select low-pLDDT disordered crops ──────────────────────
        IDP_MAX_PLDDT = args.idp_max_plddt
        IDP_MIN_AF_BIND_W = args.idp_min_af_bind_w
        MIN_DOMAIN_LENGTH = args.min_domain_length if args.min_domain_length != DEFAULT_MIN_DOMAIN_LENGTH else DEFAULT_IDP_MIN_DOMAIN_LENGTH
        MAX_DOMAIN_LENGTH = args.max_domain_length if args.max_domain_length != DEFAULT_MAX_DOMAIN_LENGTH else DEFAULT_IDP_MAX_DOMAIN_LENGTH

        print("\nApplying IDP filters:")
        if IDP_MAX_PLDDT > 0.0:
            df = _apply_filter(df, df["mean_plddt"] <= IDP_MAX_PLDDT, f"max_plddt <= {IDP_MAX_PLDDT}")
        if MIN_DOMAIN_LENGTH > 0:
            df = _apply_filter(df, df["domain_length"] >= MIN_DOMAIN_LENGTH, f"min_domain_length >= {MIN_DOMAIN_LENGTH}")
        if MAX_DOMAIN_LENGTH > 0:
            df = _apply_filter(df, df["domain_length"] <= MAX_DOMAIN_LENGTH, f"max_domain_length <= {MAX_DOMAIN_LENGTH}")
        if IDP_MIN_AF_BIND_W > 0.0 and "mean_af_bind_w" in df.columns:
            df = _apply_filter(df, df["mean_af_bind_w"] >= IDP_MIN_AF_BIND_W,
                               f"min_af_bind_w >= {IDP_MIN_AF_BIND_W}")

        config_summary = {
            "mode": "IDP",
            "idp_max_plddt": IDP_MAX_PLDDT,
            "idp_min_af_bind_w": IDP_MIN_AF_BIND_W,
            "min_domain_length": MIN_DOMAIN_LENGTH,
            "max_domain_length": MAX_DOMAIN_LENGTH,
            "num_sweep_batches": NUM_SWEEP_BATCHES,
        }

    else:
        # ── Structured mode: original filter set ─────────────────────────────
        MIN_PLDDT = args.min_plddt
        MIN_DOMAIN_LENGTH = args.min_domain_length
        MAX_DOMAIN_LENGTH = args.max_domain_length
        MAX_RADIUS_OF_GYRATION_NORMALIZED = args.max_radius_of_gyration_normalized
        MIN_CONTACTS_PER_RESIDUE = args.min_contacts_per_residue
        DROP_ONLY_SINGLE_HELIX_DOMAINS = args.drop_only_single_helix_domains
        DROP_ONLY_DOUBLE_HELIX_DOMAINS = args.drop_only_double_helix_domains
        DROP_ONLY_LOOPY_DOMAINS = args.drop_only_loopy_domains
        DROP_UNIPROT_AFDB_SEQUENCE_MISMATCHES = args.drop_uniprot_afdb_sequence_mismatches
        DROP_MEMBRANE_EMBEDDED_RESIDUES_OVERLAP_DOMAIN = args.drop_membrane_embedded_residues_overlap_domain

        print("\nApplying structured-domain filters:")
        if MIN_PLDDT > 0.0:
            df = _apply_filter(df, df["mean_plddt"] >= MIN_PLDDT, f"min_plddt >= {MIN_PLDDT}")
        if MIN_DOMAIN_LENGTH > 0:
            df = _apply_filter(df, df["domain_length"] >= MIN_DOMAIN_LENGTH, f"min_domain_length >= {MIN_DOMAIN_LENGTH}")
        if MAX_DOMAIN_LENGTH > 0:
            df = _apply_filter(df, df["domain_length"] <= MAX_DOMAIN_LENGTH, f"max_domain_length <= {MAX_DOMAIN_LENGTH}")
        if MAX_RADIUS_OF_GYRATION_NORMALIZED > 0.0:
            df = _apply_filter(df, df["radius_of_gyration_normalized"] <= MAX_RADIUS_OF_GYRATION_NORMALIZED,
                               f"max_radius_of_gyration_normalized <= {MAX_RADIUS_OF_GYRATION_NORMALIZED}")
        if MIN_CONTACTS_PER_RESIDUE > 0.0:
            df = _apply_filter(df, df["average_ca_contacts_per_residue"] >= MIN_CONTACTS_PER_RESIDUE,
                               f"min_contacts_per_residue >= {MIN_CONTACTS_PER_RESIDUE}")
        if DROP_ONLY_SINGLE_HELIX_DOMAINS:
            df = _apply_filter(df, ~df["is_only_single_helix"], "drop_only_single_helix_domains")
        if DROP_ONLY_DOUBLE_HELIX_DOMAINS:
            df = _apply_filter(df, ~df["is_only_double_helix"], "drop_only_double_helix_domains")
        if DROP_ONLY_LOOPY_DOMAINS:
            df = _apply_filter(df, ~df["is_only_loopy"], "drop_only_loopy_domains")
        if DROP_UNIPROT_AFDB_SEQUENCE_MISMATCHES:
            df = _apply_filter(df, df["uniprot_afdb_sequences_match"] == True, "drop_uniprot_afdb_sequence_mismatches")
        if DROP_MEMBRANE_EMBEDDED_RESIDUES_OVERLAP_DOMAIN:
            df = _apply_filter(df, df["domain_overlaps_membrane_embedded_residues"] == False,
                               "drop_membrane_embedded_residues_overlap_domain")

        config_summary = {
            "mode": "structured",
            "min_plddt": MIN_PLDDT,
            "min_domain_length": MIN_DOMAIN_LENGTH,
            "max_domain_length": MAX_DOMAIN_LENGTH,
            "max_radius_of_gyration_normalized": MAX_RADIUS_OF_GYRATION_NORMALIZED,
            "min_contacts_per_residue": MIN_CONTACTS_PER_RESIDUE,
            "drop_only_single_helix_domains": DROP_ONLY_SINGLE_HELIX_DOMAINS,
            "drop_only_double_helix_domains": DROP_ONLY_DOUBLE_HELIX_DOMAINS,
            "drop_only_loopy_domains": DROP_ONLY_LOOPY_DOMAINS,
            "drop_uniprot_afdb_sequence_mismatches": DROP_UNIPROT_AFDB_SEQUENCE_MISMATCHES,
            "drop_membrane_embedded_residues_overlap_domain": DROP_MEMBRANE_EMBEDDED_RESIDUES_OVERLAP_DOMAIN,
            "num_sweep_batches": NUM_SWEEP_BATCHES,
        }

    # Save filtered stats CSV
    filtered_stats_csv_path.parent.mkdir(parents=True, exist_ok=True)
    df = df.sort_values(by="domain_length", inplace=False)
    df.to_csv(filtered_stats_csv_path, index=False)
    print(f"\nWrote filtered stats: {filtered_stats_csv_path}  ({len(df)} rows)")

    # Build sweep index rows
    saved_domains = df.to_dict(orient="records")
    saved_domains = sorted(saved_domains, key=lambda r: r["domain_length"])
    rows = [
        {
            "uniprot": r["uniprot"],
            "domain_name": r["domain_name"],
            "domain_pdb_path": r["domain_pdb_path"],
            "full_pdb_path": r["full_pdb_path"],
            "json_path": r["json_path"],
            "domain_length": r["domain_length"],
        }
        for r in saved_domains
    ]
    write_domain_index_batches(sweep_index_path, config_summary, rows, NUM_SWEEP_BATCHES)

    mode_label = "IDP crops" if IDP_MODE else "domains"
    print(f"Sweep index written ({len(rows)} {mode_label}): {sweep_index_path}")
    print("Filtering completed successfully.")

if __name__ == "__main__":
    main()