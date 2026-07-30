import pandas as pd
import os
import tempfile
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.DSSP import DSSP
from tqdm.auto import tqdm
from pathlib import Path
from itertools import groupby
import argparse
import numpy as np
from scipy.spatial import cKDTree
import json
from Bio.PDB.Polypeptide import is_aa
import concurrent.futures
from dataclasses import dataclass
from functools import partial

# Inputs
DEFAULT_PDB_DIR = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "PDB"
DEFAULT_JSON_PAE_DIR = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "PAE"
DEFAULT_DOMAIN_CSV_DIR = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "CSV_Domains"
DEFAULT_UNIPROT_CSV = Path(__file__).resolve().parent.parent / "Data" / "UniProt" / "uniprot_human_reviewed.csv"
DEFAULT_LOCAL_DATA_BASE_DIR = Path(__file__).resolve().parent.parent / "Data"
DEFAULT_REMOTE_DATA_BASE_DIR = Path("/leonardo_work/AIFAC_S03_099/config/scripts/Bindome/Data")
DEFAULT_DSSP_BINARY_PATH = Path(__file__).resolve().parent.parent / "BindCraft" / "functions" / "dssp"

# Outputs
DEFAULT_PDB_DOMAIN_DIR = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "PDB_Domains"
DEFAULT_DOMAIN_STATS_PATH = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Statistics" / "domain_stats.csv"

# Config
DEFAULT_DEBUG = False
DEFAULT_BUFFER_RESIDUES = 2
DEFAULT_MERGE_GAPS_UP_TO_N_RESIDUES = 20

# AF-Bind scoring constants (Piovesan et al. 2022)
AF_BIND_RSA_THRESHOLD = 0.581
AF_BIND_PLDDT_LOW = 25.0
AF_BIND_PLDDT_HIGH = 68.8
AF_BIND_RSA_SMOOTH_WIN = 5

three_to_one_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute statistics for AFDB domains (length, pLDDT, DSSP).")
    parser.add_argument("--pdb-dir", type=str, default=str(DEFAULT_PDB_DIR), help="Path to PDB directory.")
    parser.add_argument("--domain-csv-dir", type=str, default=str(DEFAULT_DOMAIN_CSV_DIR), help="Path to domains CSV directory.")
    parser.add_argument("--json-pae-dir", type=str, default=str(DEFAULT_JSON_PAE_DIR), help="Path to PAE JSON directory.")
    parser.add_argument("--uniprot-csv", type=str, default=str(DEFAULT_UNIPROT_CSV), help="Path to UniProt CSV file.")
    parser.add_argument("--pdb-domain-dir", type=str, default=str(DEFAULT_PDB_DOMAIN_DIR),
                        help="(Kept for compatibility) Path to output PDB domains directory. Stats script does not write PDB domains.")
    parser.add_argument("--domain-stats-path", type=str, default=str(DEFAULT_DOMAIN_STATS_PATH),
                        help="Path to write the domain stats TSV.")
    parser.add_argument("--local-data-base-dir", type=str, default=str(DEFAULT_LOCAL_DATA_BASE_DIR),
                        help="Base directory for local data.")
    parser.add_argument("--remote-data-base-dir", type=str, default=str(DEFAULT_REMOTE_DATA_BASE_DIR),
                        help="Base directory for remote data.")
                        
    parser.add_argument("--buffer-residues", type=int, default=DEFAULT_BUFFER_RESIDUES,
                        help="Residues to extend domain on each side.")
    parser.add_argument("--merge-gaps-up-to-n-residues", type=int, default=DEFAULT_MERGE_GAPS_UP_TO_N_RESIDUES,
                        help="Merge gaps up to this many residues when extending domains. 0 = no merging.")
    parser.add_argument("--dssp-binary-path", type=str, default=str(DEFAULT_DSSP_BINARY_PATH),
                        help="Path to the DSSP binary.")
    parser.add_argument("--debug", action=argparse.BooleanOptionalAction, default=DEFAULT_DEBUG,
                        help="Enable debug mode.")
    parser.add_argument("--max-workers", type=int, default=None,
                        help="Maximum parallel workers for domain processing (default: SLURM, then CPU count, else 8).")
    return parser.parse_args()


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
    slurm_cpus = _slurm_cpu_count()
    if slurm_cpus:
        return slurm_cpus
    cpu_count = os.cpu_count()
    if cpu_count:
        return cpu_count
    return 8


def indices_to_ranges(indices):
    """Turn [5,6,7,10,11] → [(5,7),(10,11)]"""
    if not indices:
        return []
    vals = sorted(set(indices))
    ranges = []
    # group consecutive runs using the (value - index) trick
    for _, grp in groupby(enumerate(vals), key=lambda t: t[1] - t[0]):
        g = list(grp)
        start, end = g[0][1] + 1, g[-1][1] + 1  # convert to 1-based inclusive
        ranges.append((start, end))
    return ranges


def ranges_to_chopping(ranges):
    return "_".join(f"{s}-{e}" for s, e in ranges)


def residue_in_ranges(resseq, ranges):
    for s, e in ranges:
        if s <= resseq <= e:
            return True
    return False


def compute_domain_length(struct, chain_id, ranges):
    """
    Compute domain length as number of residues in `ranges` (inclusive, 1-based).
    """
    model = next(struct.get_models())
    chain = next((c for c in model.get_chains() if c.id == chain_id), None)
    if chain is None:
        return 0

    length = 0
    for residue in chain:
        hetflag, resseq, icode = residue.get_id()
        if hetflag != " ":
            continue
        if residue_in_ranges(resseq, ranges):
            length += 1

    return length


def compute_domain_plddt(struct, chain_id, ranges, prefer_ca=True):
    """
    Average pLDDT over residues in `ranges` (inclusive, 1-based).
    For AlphaFold PDBs, pLDDT is stored in B-factors (0..100).
    Returns (mean_plddt, n_residues_used).
    """
    model = next(struct.get_models())
    chain = next((c for c in model.get_chains() if c.id == chain_id), None)
    if chain is None:
        return float("nan"), 0

    plddts = []
    for residue in chain:
        hetflag, resseq, icode = residue.get_id()
        if hetflag != " ":
            continue
        if not residue_in_ranges(resseq, ranges):
            continue

        val = residue["CA"].get_bfactor()
        plddts.append(val)

    return float(sum(plddts) / len(plddts)), len(plddts)


def compute_domain_secondary_structure_content(struct, chain_id, ranges, dssp_path):
    """
    Compute secondary structure content for residues in `ranges` (inclusive, 1-based)
    using DSSP.

    Returns a dict that is convenient for filtering, e.g.:
    {
        "n": 123,
        "helix_frac": 0.41,
        "strand_frac": 0.22,
        "coil_frac": 0.37,
        "helix_segments": 3,
        "strand_segments": 2,
        "coil_segments": 4,
        "is_only_single_helix": False,
        "is_only_double_helix": False,
        "is_only_loopy": False,
    }

    Notes:
    - Requires DSSP installed on the system (mkdssp or dssp binary in PATH).
    - DSSP codes: H/G/I = helix; E/B = beta; everything else = coil/loop/turn.
    """
    model = next(struct.get_models())

    # DSSP
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        tmp_pdb = tmp.name
    io = PDBIO()
    io.set_structure(struct)
    io.save(tmp_pdb)
    dssp = DSSP(
        model,
        tmp_pdb,
        dssp=str(dssp_path)
    )
    helix_codes = {"H", "G", "I"}
    strand_codes = {"E", "B"}

    # Collect per-residue DSSP assignments within the domain ranges
    ss_by_resseq = []
    rsa_by_resseq = []  # (resseq, rsa) for AF-Bind scoring
    for (dssp_chain, res_id) in dssp.keys():
        if dssp_chain != chain_id:
            continue
        hetflag, resseq, icode = res_id
        if hetflag != " ":
            continue
        if not residue_in_ranges(resseq, ranges):
            continue

        dssp_entry = dssp[(dssp_chain, res_id)]
        ss = dssp_entry[2]  # DSSP secondary structure code
        # DSSP can return ' ' for unknown; normalize that to coil
        if ss is None or ss == " ":
            ss = "-"
        ss_by_resseq.append((resseq, ss))
        rsa = dssp_entry[3]  # relative ASA (0-1)
        rsa_by_resseq.append((resseq, float(rsa) if rsa is not None else 0.0))

    if not ss_by_resseq:
        raise ValueError("DSSP returned no secondary structure assignments for the specified domain.")

    # Classify into 3-state labels to make filtering easier
    labels = []
    for _, ss in ss_by_resseq:
        if ss in helix_codes:
            labels.append("H")
        elif ss in strand_codes:
            labels.append("E")
        else:
            labels.append("C")

    dssp_length = len(labels)
    n_helix = sum(1 for x in labels if x == "H")
    n_strand = sum(1 for x in labels if x == "E")
    n_coil = sum(1 for x in labels if x == "C")

    # Count contiguous segments in the 3-state label string
    def count_segments(lbls, target):
        segs = 0
        in_seg = False
        for x in lbls:
            if x == target:
                if not in_seg:
                    segs += 1
                    in_seg = True
            else:
                in_seg = False
        return segs

    n_helix_segments = count_segments(labels, "H")
    n_strand_segments = count_segments(labels, "E")
    n_coil_segments = count_segments(labels, "C")

    # Booleans for filtering:
    # - "only single helix": all residues are helix or loopy AND exactly one helix segment
    # - "only double helix": there is no strand segment AND exactly two helix segments, everything else must be loopy
    # - "only loopy": no helix and no strand (everything coil)
    is_only_single_helix = (n_strand_segments == 0) and (n_helix_segments == 1)
    is_only_double_helix = (n_strand_segments == 0) and (n_helix_segments == 2)
    is_only_loopy = (n_coil == dssp_length)

    return {
        "dssp_length": dssp_length,
        "helix_frac": n_helix / dssp_length,
        "strand_frac": n_strand / dssp_length,
        "coil_frac": n_coil / dssp_length,
        "n_helix_segments": n_helix_segments,
        "n_strand_segments": n_strand_segments,
        "n_coil_segments": n_coil_segments,
        "is_only_single_helix": is_only_single_helix,
        "is_only_double_helix": is_only_double_helix,
        "is_only_loopy": is_only_loopy,
        "_per_residue_rsa": rsa_by_resseq,  # internal; consumed by compute_af_bind_stats
    }

def compute_af_bind_stats(per_residue_rsa: list, struct, chain_id: str, ranges: list) -> dict:
    """
    Compute domain-level AF-Bind statistics (Piovesan et al. 2022).

    per_residue_rsa: list of (resseq, rsa) tuples from DSSP, already filtered to domain ranges.
    Returns {"mean_af_bind_w": float, "af_bind_frac": float}.
    """
    if not per_residue_rsa:
        return {"mean_af_bind_w": 0.0, "af_bind_frac": 0.0}

    model = next(struct.get_models())
    chain = next((c for c in model.get_chains() if c.id == chain_id), None)
    plddt_map = {}
    if chain is not None:
        for res in chain:
            hetflag, resseq, _ = res.get_id()
            if hetflag != " " or not residue_in_ranges(resseq, ranges):
                continue
            if "CA" in res:
                plddt_map[resseq] = res["CA"].get_bfactor()

    # Smooth RSA over a 5-residue window (matches Piovesan et al.)
    rsa_arr = np.array([rsa for _, rsa in per_residue_rsa])
    w = AF_BIND_RSA_SMOOTH_WIN
    rsa_smooth = np.convolve(rsa_arr, np.ones(w) / w, mode="same")

    af_bind_w_vals = []
    af_bind_flags = []
    for i, (resseq, rsa) in enumerate(per_residue_rsa):
        p = plddt_map.get(resseq, 50.0)
        rs = float(rsa_smooth[i])
        in_morf = AF_BIND_PLDDT_LOW < p < AF_BIND_PLDDT_HIGH
        above_rsa = rsa > AF_BIND_RSA_THRESHOLD
        af_bind_flags.append(1 if (in_morf and above_rsa) else 0)
        if in_morf and rs > AF_BIND_RSA_THRESHOLD * 0.7:
            af_bind_w_vals.append((p / AF_BIND_PLDDT_HIGH) * (rs / AF_BIND_RSA_THRESHOLD))
        else:
            af_bind_w_vals.append(0.0)

    return {
        "mean_af_bind_w": float(np.mean(af_bind_w_vals)),
        "af_bind_frac": float(np.mean(af_bind_flags)),
    }


def compute_rg(ca_coords: np.ndarray) -> float:
    """
    ca_coords: shape (N, 3)
    Rg = sqrt(mean(||r_i - r_cm||^2))
    """
    r_cm = ca_coords.mean(axis=0)
    dif = ca_coords - r_cm
    rg2 = np.mean(np.sum(dif * dif, axis=1))
    return float(np.sqrt(rg2))

def compute_radius_of_gyration(struct, chain_id, ranges):
    model = next(struct.get_models())
    chain = next((c for c in model.get_chains() if c.id == chain_id), None)
    if chain is None:
        return float("nan"), float("nan")

    residues = []
    for residue in chain:
        hetflag, resseq, _ = residue.get_id()
        if hetflag != " ":
            continue
        if not residue_in_ranges(resseq, ranges):
            continue
        if "CA" in residue:
            residues.append(residue)

    n = len(residues)
    if n == 0:
        return float("nan"), float("nan")

    coords = np.array([res["CA"].get_coord() for res in residues], dtype=float)
    rg = compute_rg(coords)
    rg_norm = rg / (n ** (1.0 / 3.0))
    return {"radius_of_gyration": rg, "radius_of_gyration_normalized": rg_norm}

def compute_num_ca_contacts(struct, chain_id, ranges, distance_cutoff=8.0):
    model = next(struct.get_models())
    chain = next((c for c in model.get_chains() if c.id == chain_id), None)

    residues = []
    for residue in chain:
        hetflag, resseq, _ = residue.get_id()
        if hetflag != " ":
            continue
        if not residue_in_ranges(resseq, ranges):
            continue
        if "CA" in residue:
            residues.append(residue)

    coords = np.array([res["CA"].get_coord() for res in residues], dtype=float)
    tree = cKDTree(coords)

    pairs = tree.query_pairs(r=distance_cutoff)
    pairs = {(i, j) for (i, j) in pairs if (j - i) > 2}
    return {"num_ca_contacts": len(pairs),
            "average_ca_contacts_per_residue": len(pairs) / len(residues) if len(residues) > 0 else 0}

def compute_membrane_annotations(uniprot_annotations, struct, chain_id, ranges):
    uniprot_sequence = (uniprot_annotations.get("sequence") or "").strip()

    def _chain_sequence(struct, chain_id: str) -> str:
        model = next(struct.get_models())
        chain = next((c for c in model.get_chains() if c.id == chain_id), None)
        if chain is None:
            return ""
        seq = []
        for residue in chain:
            hetflag, resseq, icode = residue.get_id()
            if hetflag != " ":
                continue
            if not is_aa(residue, standard=True):
                continue
            try:
                seq.append(three_to_one_map[residue.get_resname()])
            except KeyError:
                continue
        return "".join(seq)

    afdb_sequence = _chain_sequence(struct, chain_id)
    uniprot_afdb_sequences_match = uniprot_sequence == afdb_sequence

    is_membrane_protein = bool(uniprot_annotations.get("membrane_protein"))
    membrane_embedded_residue_ranges = json.loads(uniprot_annotations.get("membrane_embedded_ranges") or "[]")
    has_membrane_embedded_residue_annotation = len(membrane_embedded_residue_ranges) > 0

    def _ranges_overlap(a_start, a_end, b_start, b_end) -> bool:
        return a_start <= b_end and b_start <= a_end

    domain_overlaps_membrane_embedded_residues = False
    for d_start, d_end in ranges:
        for m_start, m_end in membrane_embedded_residue_ranges:
            if _ranges_overlap(d_start, d_end, m_start, m_end):
                domain_overlaps_membrane_embedded_residues = True
                break
        if domain_overlaps_membrane_embedded_residues:
            break

    return {
        "is_membrane_protein": is_membrane_protein,
        "has_membrane_embedded_residue_annotation": has_membrane_embedded_residue_annotation,
        "membrane_embedded_residue_ranges": membrane_embedded_residue_ranges,
        "domain_overlaps_membrane_embedded_residues": domain_overlaps_membrane_embedded_residues,
        "uniprot_afdb_sequences_match": uniprot_afdb_sequences_match,
    }

class DomainSelect(Select):
    """Select only residues in the specified numeric ranges for one chain."""
    def __init__(self, chain_id, ranges):
        """
        chain_id : str, e.g. 'A'
        ranges   : list of (start, end) tuples, inclusive
        """
        self.chain_id = chain_id
        self.ranges = ranges

    def accept_chain(self, chain):
        return (chain.id == self.chain_id)

    def accept_residue(self, residue):
        # skip HETATM or water
        hetflag, resseq, icode = residue.get_id()
        if hetflag != " ":
            return False
        # check if its residue number is in any of our ranges
        for start, end in self.ranges:
            if start <= resseq <= end:
                return True
        return False



@dataclass(frozen=True)
class DomainProcessingConfig:
    pdb_folder: str
    output_pdb_folder: str
    json_pae_dir: str
    num_buffer_residues: int
    merge_gaps_up_to_n_residues: int
    dssp_binary_path: str
    uniprot_dict: dict
    local_base: str
    remote_base: str

    def to_remote(self, path: str) -> str:
        if self.local_base and self.remote_base:
            return path.replace(self.local_base, self.remote_base)
        return path


def _process_domain(domain: dict, config: DomainProcessingConfig) -> dict | None:
    uniprot = domain["target_name"]
    uniprot_annotations = config.uniprot_dict.get(uniprot, {})
    domain_name = domain["domain_name"]
    chopping = domain["chopping"]
    pdb_fp = os.path.join(config.pdb_folder, uniprot + ".pdb")

    ranges = [tuple(map(int, part.split("-"))) for part in chopping.split("_")]
    ranges = [
        (max(1, start - config.num_buffer_residues), end + config.num_buffer_residues)
        for start, end in ranges
    ]
    if config.merge_gaps_up_to_n_residues > 0:
        ranges = sorted(ranges, key=lambda r: r[0])
        merged = []
        for start, end in ranges:
            if not merged:
                merged.append([start, end])
                continue
            prev_start, prev_end = merged[-1]
            gap = start - prev_end - 1
            if gap <= config.merge_gaps_up_to_n_residues:
                merged[-1][1] = max(prev_end, end)
            else:
                merged.append([start, end])
        ranges = [(s, e) for s, e in merged]

    local_parser = PDBParser(QUIET=True)
    local_io = PDBIO()
    struct = local_parser.get_structure(domain_name, pdb_fp)
    model = next(struct.get_models())
    chain = next(model.get_chains())
    chain_id = chain.id

    selector = DomainSelect(chain_id, ranges)
    domain_length = sum(1 for residue in chain if selector.accept_residue(residue))
    out_fp = os.path.join(config.output_pdb_folder, f"{domain_name}_l{domain_length}.pdb")
    local_io.set_structure(struct)
    local_io.save(out_fp, selector)

    domain_length = compute_domain_length(struct, chain_id, ranges)
    mean_plddt, n_used = compute_domain_plddt(struct, chain_id, ranges)
    secondary_structure = compute_domain_secondary_structure_content(
        struct, chain_id, ranges, config.dssp_binary_path
    )
    per_residue_rsa = secondary_structure.pop("_per_residue_rsa", [])
    af_bind_stats = compute_af_bind_stats(per_residue_rsa, struct, chain_id, ranges)
    rg_stats = compute_radius_of_gyration(struct, chain_id, ranges)
    num_ca_contacts = compute_num_ca_contacts(struct, chain_id, ranges)
    membrane_infos = compute_membrane_annotations(uniprot_annotations, struct, chain_id, ranges)

    domain_pdb_path_local = str(Path(out_fp).resolve())
    full_pdb_path_local = str(Path(pdb_fp).resolve())
    json_path_local = str((Path(config.json_pae_dir) / f"{uniprot}.pae.json").resolve())

    return {
        "uniprot": uniprot,
        "domain_name": domain_name,
        "domain_pdb_path_local": domain_pdb_path_local,
        "full_pdb_path_local": full_pdb_path_local,
        "json_path_local": json_path_local,
        "domain_pdb_path": config.to_remote(domain_pdb_path_local),
        "full_pdb_path": config.to_remote(full_pdb_path_local),
        "json_path": config.to_remote(json_path_local),
        "domain_length": domain_length,
        "design_time": 0,
        "design_runs": 0,
        "num_accepted_designs": 0,
        "num_trajectories": 0,
        "default_chopping": chopping,
        "chopping": ranges_to_chopping(ranges),
        "chain_id": chain_id,
        "buffer_residues": config.num_buffer_residues,
        "mean_plddt": mean_plddt,
        "n_plddt_used": n_used,
        **secondary_structure,
        **rg_stats,
        **num_ca_contacts,
        **membrane_infos,
        **af_bind_stats,
    }

def main() -> None:
    args = _parse_args()
    print("Parsed arguments:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")

    pdb_folder = Path(args.pdb_dir)
    domain_csv_fps = Path(args.domain_csv_dir)
    output_pdb_folder = Path(args.pdb_domain_dir)
    json_pae_dir = Path(args.json_pae_dir)

    local_data_base_dir = Path(args.local_data_base_dir)
    remote_data_base_dir = Path(args.remote_data_base_dir)

    dssp_binary_path = args.dssp_binary_path
    num_buffer_residues_around_domain = args.buffer_residues
    merge_gaps_up_to_n_residues = args.merge_gaps_up_to_n_residues
    DEBUG = args.debug
    max_workers = _resolve_max_workers(args.max_workers)
    print(f"Using max_workers={max_workers} for domain processing.")
    io = PDBIO()
    parser = PDBParser(QUIET=True)

    os.makedirs(output_pdb_folder, exist_ok=True)

    # Load domain CSVs
    csv_fps = sorted(domain_csv_fps.glob("**/*.domains.csv"))
    loaded = []
    for fp in tqdm(csv_fps):
        lines = Path(fp).read_text().splitlines()
        df = pd.DataFrame({"domain_residue_idx": lines})
        df["target_name"] = fp.name.split(".")[0]
        df["domain_number"] = range(1, len(df) + 1)
        df["domain_name"] = df["target_name"] + f"_DOMAIN" + df["domain_number"].astype(str)
        df["domain_ranges"] = df["domain_residue_idx"].apply(
            lambda s: indices_to_ranges([int(x) for x in s.strip().split(",") if x != ""])
        )
        df["chopping"] = df["domain_ranges"].apply(ranges_to_chopping)
        loaded.append(df)

    domains = pd.concat(loaded, ignore_index=True, sort=False)

    # Load UniProt CSV for membrane information
    uniprot_df = pd.read_csv(args.uniprot_csv)
    uniprot_dict = uniprot_df.set_index("uniprot_id").to_dict(orient="index")

    local_base = str(local_data_base_dir.resolve())
    remote_base = str(remote_data_base_dir.resolve())
    config = DomainProcessingConfig(
        pdb_folder=str(pdb_folder),
        output_pdb_folder=str(output_pdb_folder),
        json_pae_dir=str(json_pae_dir),
        num_buffer_residues=num_buffer_residues_around_domain,
        merge_gaps_up_to_n_residues=merge_gaps_up_to_n_residues,
        dssp_binary_path=str(dssp_binary_path),
        uniprot_dict=uniprot_dict,
        local_base=local_base,
        remote_base=remote_base,
    )

    # Process domains
    domains_records = domains.to_dict(orient="records")
    if DEBUG:
        domains_records = domains_records[:150]

    rows = []
    print(f"Processing {len(domains_records)} domains with up to {max_workers} parallel workers...")
    process_func = partial(_process_domain, config=config)
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        for row in tqdm(
            ex.map(process_func, domains_records),
            total=len(domains_records),
            desc="Computing domain stats"
        ):
            if row:
                rows.append(row)

    stats_path = Path(args.domain_stats_path)
    stats_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(rows)
    df.to_csv(stats_path, index=False)

    print(f"Wrote domain stats CSV: {stats_path}")
    print("Computing domain stats completed successfully.")

if __name__ == "__main__":
    main()
