#!/usr/bin/env python3
"""
Step 04b: Generate IDP-specific crops for BindCraft.

Reads domain_stats.csv (from step 04), identifies IDP-like domains
(mean_plddt < --idp-max-plddt AND coil_frac >= --min-coil-frac), and for each:
  1. Runs AF-Bind window scoring (Piovesan et al. 2022) on the full-length PDB
  2. Selects non-overlapping crops greedily by AF-Bind score
  3. Writes cropped PDB files to --crops-dir
  4. Writes BindCraft target settings JSONs to --settings-dir

Outputs:
  --idp-crops-stats-path  (CSV, analogous to domain_stats.csv — input for step 05 --idp)
  --bindcraft-jobs-file   (text file, one settings JSON path per line)

Usage:
  python 04b_generate_idp_crops.py \\
      --stats-csv /path/to/domain_stats.csv \\
      --crops-dir /path/to/PDB_Crops \\
      --settings-dir /path/to/Settings_IDP \\
      --design-root /path/to/output_idp \\
      --bindcraft-advanced-settings /path/to/idp_4stage_multimer_iface_mv_strongiface.json
"""

import argparse
import csv
import json
import os
import tempfile
import concurrent.futures
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.DSSP import DSSP
from tqdm.auto import tqdm

# ── defaults ──────────────────────────────────────────────────────────────────
DEFAULT_DOMAIN_STATS_PATH = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Statistics" / "domain_stats.csv"
DEFAULT_CROPS_DIR = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "PDB_IDP_Crops"
DEFAULT_SETTINGS_DIR = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Settings_IDP"
DEFAULT_IDP_CROPS_STATS_PATH = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Statistics" / "idp_crops_stats.csv"
DEFAULT_BINDCRAFT_JOBS_FILE = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "idp_bindcraft_jobs.txt"
DEFAULT_DSSP_BINARY_PATH = Path(__file__).resolve().parent.parent / "BindCraft" / "functions" / "dssp"
DEFAULT_LOCAL_DATA_BASE_DIR = Path(__file__).resolve().parent.parent / "Data"
DEFAULT_REMOTE_DATA_BASE_DIR = Path("/leonardo_work/AIFAC_S03_099/config/scripts/Bindome/Data")

# IDP domain selection thresholds
DEFAULT_IDP_MAX_PLDDT = 50.0
DEFAULT_MIN_COIL_FRAC = 0.5

# Crop generation parameters
DEFAULT_WINDOW_SIZES = [20, 25]
DEFAULT_STRIDE_FRAC = 0.5
DEFAULT_MIN_SPACING = 5
DEFAULT_BINDER_LENGTH = 90
DEFAULT_N_FINAL_DESIGNS = 5

# AF-Bind constants (Piovesan et al. 2022)
AF_BIND_RSA_THRESHOLD = 0.581
AF_BIND_PLDDT_LOW = 25.0
AF_BIND_PLDDT_HIGH = 68.8
AF_BIND_RSA_SMOOTH_WIN = 5


# ── argument parsing ──────────────────────────────────────────────────────────
def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stats-csv", type=str, default=str(DEFAULT_DOMAIN_STATS_PATH),
                        help="domain_stats.csv from step 04.")
    parser.add_argument("--crops-dir", type=str, default=str(DEFAULT_CROPS_DIR),
                        help="Output directory for cropped PDB files.")
    parser.add_argument("--settings-dir", type=str, default=str(DEFAULT_SETTINGS_DIR),
                        help="Output directory for BindCraft target settings JSONs.")
    parser.add_argument("--design-root", type=str, required=True,
                        help="Root directory for BindCraft design output (design_path in settings JSON).")
    parser.add_argument("--bindcraft-advanced-settings", type=str, default="",
                        help="Path to BindCraft advanced settings JSON (e.g. strongiface+mv). "
                             "Recorded in crops stats for downstream use.")
    parser.add_argument("--idp-crops-stats-path", type=str, default=str(DEFAULT_IDP_CROPS_STATS_PATH),
                        help="Output CSV path for IDP crop stats (input for step 05 --idp).")
    parser.add_argument("--bindcraft-jobs-file", type=str, default=str(DEFAULT_BINDCRAFT_JOBS_FILE),
                        help="Output text file listing settings JSON paths (one per line).")
    parser.add_argument("--dssp-binary-path", type=str, default=str(DEFAULT_DSSP_BINARY_PATH),
                        help="Path to the DSSP binary.")
    parser.add_argument("--idp-max-plddt", type=float, default=DEFAULT_IDP_MAX_PLDDT,
                        help="Max mean pLDDT for a domain to be treated as IDP.")
    parser.add_argument("--min-coil-frac", type=float, default=DEFAULT_MIN_COIL_FRAC,
                        help="Min coil fraction for a domain to be treated as IDP.")
    parser.add_argument("--window-sizes", nargs="+", type=int, default=DEFAULT_WINDOW_SIZES,
                        help="Window sizes to score (default: 20 25).")
    parser.add_argument("--stride-frac", type=float, default=DEFAULT_STRIDE_FRAC,
                        help="Stride as fraction of window size.")
    parser.add_argument("--min-spacing", type=int, default=DEFAULT_MIN_SPACING,
                        help="Minimum residue gap between selected crops.")
    parser.add_argument("--binder-length", type=int, default=DEFAULT_BINDER_LENGTH,
                        help="Binder length for BindCraft settings JSON.")
    parser.add_argument("--n-final-designs", type=int, default=DEFAULT_N_FINAL_DESIGNS,
                        help="number_of_final_designs in BindCraft settings JSON.")
    parser.add_argument("--local-data-base-dir", type=str, default=str(DEFAULT_LOCAL_DATA_BASE_DIR),
                        help="Base directory for local paths (for local→remote mapping).")
    parser.add_argument("--remote-data-base-dir", type=str, default=str(DEFAULT_REMOTE_DATA_BASE_DIR),
                        help="Base directory for remote paths.")
    parser.add_argument("--max-workers", type=int, default=None,
                        help="Max parallel workers (default: SLURM, then CPU count, else 8).")
    parser.add_argument("--debug", action=argparse.BooleanOptionalAction, default=False,
                        help="Process only the first 10 IDP domains.")
    return parser.parse_args()


# ── helpers ───────────────────────────────────────────────────────────────────
def _slurm_cpu_count() -> int | None:
    for key in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE"):
        val = os.environ.get(key)
        if not val:
            continue
        try:
            return max(1, int(str(val).split(",")[0]))
        except ValueError:
            continue
    return None


def _resolve_max_workers(requested: int | None) -> int:
    if requested and requested > 0:
        return requested
    slurm = _slurm_cpu_count()
    if slurm:
        return slurm
    cpu = os.cpu_count()
    return cpu if cpu else 8


class _ResidueRangeSelect(Select):
    def __init__(self, chain_id: str, resnums: set):
        self.chain_id = chain_id
        self.resnums = resnums

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, res):
        hetflag, resseq, _ = res.get_id()
        return hetflag == " " and resseq in self.resnums


# ── per-residue AF-Bind scoring ───────────────────────────────────────────────
def _score_residues(pdb_path: str, chain_id: str, resnum_range: tuple[int, int],
                    dssp_bin: str) -> list[dict]:
    """
    Score residues in [resnum_range[0], resnum_range[1]] from pdb_path.
    Runs DSSP on the full PDB (correct RSA context) and filters to the range.
    Returns list of dicts: resnum, plddt, rsa, rsa_smooth, af_bind_w.
    """
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    model = struct[0]

    plddt_map = {}
    chain = next((c for c in model.get_chains() if c.id == chain_id), None)
    if chain is None:
        return []
    for res in chain:
        hetflag, resseq, _ = res.get_id()
        if hetflag == " " and "CA" in res:
            plddt_map[resseq] = res["CA"].get_bfactor()

    tmp_pdb = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp_pdb = tmp.name
        io = PDBIO()
        io.set_structure(struct)
        io.save(tmp_pdb)
        dssp = DSSP(model, tmp_pdb, dssp=dssp_bin)
    finally:
        if tmp_pdb and os.path.exists(tmp_pdb):
            os.unlink(tmp_pdb)

    r_start, r_end = resnum_range
    records = []
    for (dc, res_id) in dssp.keys():
        if dc != chain_id:
            continue
        hetflag, resseq, _ = res_id
        if hetflag != " ":
            continue
        if not (r_start <= resseq <= r_end):
            continue
        entry = dssp[(dc, res_id)]
        rsa = entry[3]
        rsa = float(rsa) if rsa is not None else 0.0
        records.append({
            "resnum": resseq,
            "plddt": plddt_map.get(resseq, 50.0),
            "rsa": rsa,
        })

    if not records:
        return []

    # Sort by residue number (DSSP order should already be correct, but be safe)
    records.sort(key=lambda r: r["resnum"])

    # Smooth RSA with a 5-residue window
    rsa_arr = np.array([r["rsa"] for r in records])
    w = AF_BIND_RSA_SMOOTH_WIN
    rsa_smooth = np.convolve(rsa_arr, np.ones(w) / w, mode="same")

    for i, rec in enumerate(records):
        p = rec["plddt"]
        rs = float(rsa_smooth[i])
        in_morf = AF_BIND_PLDDT_LOW < p < AF_BIND_PLDDT_HIGH
        if in_morf and rs > AF_BIND_RSA_THRESHOLD * 0.7:
            rec["af_bind_w"] = (p / AF_BIND_PLDDT_HIGH) * (rs / AF_BIND_RSA_THRESHOLD)
        else:
            rec["af_bind_w"] = 0.0

    return records


# ── window scoring ────────────────────────────────────────────────────────────
def _score_windows(records: list[dict], window_sizes: list[int],
                   stride_frac: float) -> list[dict]:
    """Score all sliding windows; returns list sorted by score descending."""
    if not records:
        return []
    af_bind_w = np.array([r["af_bind_w"] for r in records])
    plddt_arr = np.array([r["plddt"] for r in records])
    n = len(records)
    windows = []

    for wlen in window_sizes:
        if wlen > n:
            continue
        stride = max(1, int(wlen * stride_frac))
        for start in range(0, n - wlen + 1, stride):
            end = start + wlen - 1
            w_plddt = plddt_arr[start:start + wlen]
            mean_plddt = float(w_plddt.mean())
            score = float(af_bind_w[start:start + wlen].mean())
            boundary_penalty = float(w_plddt.std()) / 100.0
            windows.append({
                "start_idx": start,
                "end_idx": end,
                "resnum_start": records[start]["resnum"],
                "resnum_end": records[end]["resnum"],
                "window_size": wlen,
                "mean_plddt": round(mean_plddt, 2),
                "mean_af_bind_w": round(score, 4),
                "af_bind_frac": round(
                    float(np.mean([1 if r["af_bind_w"] > 0 else 0
                                   for r in records[start:start + wlen]])), 4),
                "score": round(score - 0.5 * boundary_penalty, 4),
            })

    windows.sort(key=lambda x: -x["score"])
    return windows


# ── greedy crop selection ─────────────────────────────────────────────────────
def _select_crops(windows: list[dict], n_residues: int, min_spacing: int) -> list[dict]:
    """
    Greedy: pick highest-scoring window, mask ± min_spacing, repeat.
    Only selects windows with score > 0 (i.e. with non-zero AF-Bind content).
    """
    masked = np.zeros(n_residues, dtype=bool)
    selected = []
    for w in windows:
        if w["score"] <= 0.0:
            continue
        s, e = w["start_idx"], w["end_idx"]
        if masked[s:e + 1].any():
            continue
        selected.append(w)
        m_start = max(0, s - min_spacing)
        m_end = min(n_residues, e + min_spacing + 1)
        masked[m_start:m_end] = True
    return selected


# ── PDB writing ───────────────────────────────────────────────────────────────
def _write_crop_pdb(full_pdb: str, chain_id: str,
                    resnum_start: int, resnum_end: int, out_path: str) -> None:
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", full_pdb)
    resnums = set(range(resnum_start, resnum_end + 1))
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path, _ResidueRangeSelect(chain_id, resnums))


def _write_settings_json(crop_name: str, crop_pdb: str, design_root: str,
                         binder_length: int, n_final: int, out_path: str) -> None:
    settings = {
        "design_path": os.path.join(design_root, crop_name, ""),
        "binder_name": crop_name,
        "starting_pdb": str(Path(crop_pdb).resolve()),
        "chains": "A",
        "target_hotspot_residues": None,
        "lengths": [binder_length, binder_length],
        "number_of_final_designs": n_final,
    }
    with open(out_path, "w") as f:
        json.dump(settings, f, indent=4)


# ── per-domain processing ─────────────────────────────────────────────────────
@dataclass(frozen=True)
class IdpCropConfig:
    crops_dir: str
    settings_dir: str
    design_root: str
    advanced_settings: str
    dssp_binary_path: str
    window_sizes: tuple
    stride_frac: float
    min_spacing: int
    binder_length: int
    n_final_designs: int
    local_base: str
    remote_base: str

    def to_remote(self, path: str) -> str:
        if self.local_base and self.remote_base:
            return path.replace(self.local_base, self.remote_base)
        return path


def _process_idp_domain(domain_row: dict, config: IdpCropConfig) -> list[dict]:
    """
    Process one IDP domain: score windows → select crops → write files.
    Returns a list of crop dicts (empty if no suitable crops found).
    """
    full_pdb = domain_row.get("full_pdb_path_local", "")
    if not full_pdb or not os.path.exists(full_pdb):
        return []

    chain_id = str(domain_row.get("chain_id", "A"))
    chopping = str(domain_row.get("chopping", ""))

    # Parse domain residue range from chopping string (e.g. "100-250" or "100-150_200-250")
    try:
        parts = [tuple(map(int, seg.split("-"))) for seg in chopping.split("_")]
        r_start = min(s for s, e in parts)
        r_end = max(e for s, e in parts)
    except Exception:
        return []

    records = _score_residues(full_pdb, chain_id, (r_start, r_end), config.dssp_binary_path)
    if not records:
        return []

    windows = _score_windows(records, list(config.window_sizes), config.stride_frac)
    if not windows:
        return []

    selected = _select_crops(windows, len(records), config.min_spacing)
    if not selected:
        return []

    uniprot = domain_row["uniprot"]
    domain_name = domain_row["domain_name"]
    pae_json = domain_row.get("json_path_local", domain_row.get("json_path", ""))
    full_pdb_remote = domain_row.get("full_pdb_path", config.to_remote(str(Path(full_pdb).resolve())))

    crops = []
    for crop in selected:
        crop_name = f"{domain_name}_crop{crop['resnum_start']:04d}_{crop['resnum_end']:04d}"
        crop_pdb_local = str(Path(os.path.join(config.crops_dir, f"{crop_name}.pdb")).resolve())
        settings_json_local = str(Path(os.path.join(config.settings_dir, f"{crop_name}_idp.json")).resolve())

        if not os.path.exists(crop_pdb_local):
            _write_crop_pdb(full_pdb, chain_id,
                            crop["resnum_start"], crop["resnum_end"], crop_pdb_local)

        if not os.path.exists(settings_json_local):
            _write_settings_json(crop_name, crop_pdb_local, config.design_root,
                                 config.binder_length, config.n_final_designs, settings_json_local)

        crops.append({
            # Core columns matching domain_stats.csv schema (for step 05 compatibility)
            "uniprot": uniprot,
            "domain_name": crop_name,
            "domain_pdb_path_local": crop_pdb_local,
            "domain_pdb_path": config.to_remote(crop_pdb_local),
            "full_pdb_path_local": str(Path(full_pdb).resolve()),
            "full_pdb_path": full_pdb_remote,
            "json_path_local": pae_json,
            "json_path": config.to_remote(pae_json) if pae_json else "",
            "domain_length": crop["window_size"],
            "chain_id": chain_id,
            # IDP-specific crop columns
            "resnum_start": crop["resnum_start"],
            "resnum_end": crop["resnum_end"],
            "mean_plddt": crop["mean_plddt"],
            "mean_af_bind_w": crop["mean_af_bind_w"],
            "af_bind_frac": crop["af_bind_frac"],
            "window_score": crop["score"],
            "bindcraft_settings_json": settings_json_local,
            "bindcraft_advanced_settings": config.advanced_settings,
            # Parent domain metadata (for traceability)
            "parent_domain_name": domain_name,
            "parent_mean_plddt": domain_row.get("mean_plddt", float("nan")),
            "parent_coil_frac": domain_row.get("coil_frac", float("nan")),
        })

    return crops


# ── main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    args = _parse_args()
    print("Parsed arguments:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")

    stats_csv = Path(args.stats_csv)
    if not stats_csv.exists():
        raise FileNotFoundError(f"domain_stats.csv not found: {stats_csv}")

    df = pd.read_csv(stats_csv)
    print(f"\nLoaded {len(df)} domains from {stats_csv}")

    # Filter to IDP-like domains
    idp_mask = (df["mean_plddt"] < args.idp_max_plddt) & (df["coil_frac"] >= args.min_coil_frac)
    idp_df = df[idp_mask].copy()
    print(f"IDP-like domains (pLDDT < {args.idp_max_plddt} AND coil_frac >= {args.min_coil_frac}): {len(idp_df)}")

    if len(idp_df) == 0:
        print("No IDP domains found. Exiting.")
        return

    if args.debug:
        idp_df = idp_df.head(10)
        print(f"[debug] Limiting to {len(idp_df)} domains.")

    crops_dir = Path(args.crops_dir)
    settings_dir = Path(args.settings_dir)
    crops_dir.mkdir(parents=True, exist_ok=True)
    settings_dir.mkdir(parents=True, exist_ok=True)

    local_base = str(Path(args.local_data_base_dir).resolve())
    remote_base = str(Path(args.remote_data_base_dir).resolve())

    config = IdpCropConfig(
        crops_dir=str(crops_dir),
        settings_dir=str(settings_dir),
        design_root=args.design_root,
        advanced_settings=args.bindcraft_advanced_settings,
        dssp_binary_path=args.dssp_binary_path,
        window_sizes=tuple(args.window_sizes),
        stride_frac=args.stride_frac,
        min_spacing=args.min_spacing,
        binder_length=args.binder_length,
        n_final_designs=args.n_final_designs,
        local_base=local_base,
        remote_base=remote_base,
    )

    domain_records = idp_df.to_dict(orient="records")
    max_workers = _resolve_max_workers(args.max_workers)
    print(f"Processing {len(domain_records)} IDP domains with up to {max_workers} workers...")

    process_func = partial(_process_idp_domain, config=config)
    all_crops: list[dict] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        for crop_list in tqdm(
            ex.map(process_func, domain_records),
            total=len(domain_records),
            desc="Generating IDP crops",
        ):
            all_crops.extend(crop_list)

    print(f"\nTotal crops generated: {len(all_crops)}")

    # Write idp_crops_stats.csv
    crops_stats_path = Path(args.idp_crops_stats_path)
    crops_stats_path.parent.mkdir(parents=True, exist_ok=True)
    crops_df = pd.DataFrame(all_crops)
    crops_df.to_csv(crops_stats_path, index=False)
    print(f"Wrote IDP crop stats: {crops_stats_path}")

    # Write bindcraft jobs file (one settings JSON path per line)
    jobs_file = Path(args.bindcraft_jobs_file)
    jobs_file.parent.mkdir(parents=True, exist_ok=True)
    with jobs_file.open("w") as f:
        for crop in all_crops:
            f.write(crop["bindcraft_settings_json"] + "\n")
    print(f"Wrote BindCraft jobs file: {jobs_file}  ({len(all_crops)} entries)")

    # Summary
    domains_with_crops = len({c["parent_domain_name"] for c in all_crops})
    print(f"\nSummary:")
    print(f"  IDP domains processed: {len(domain_records)}")
    print(f"  Domains with ≥1 crop:  {domains_with_crops}")
    print(f"  Total crops:           {len(all_crops)}")
    if all_crops:
        af_vals = [c["mean_af_bind_w"] for c in all_crops]
        print(f"  AF-Bind w mean/max:    {np.mean(af_vals):.3f} / {np.max(af_vals):.3f}")


if __name__ == "__main__":
    main()
