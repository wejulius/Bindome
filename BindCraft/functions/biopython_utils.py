####################################
################ BioPython functions
####################################
### Import dependencies
import os
import math
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List, Dict
from collections import defaultdict
from scipy.spatial import cKDTree
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, Superimposer
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import is_aa
from .generic_utils import apply_advanced_rounding

three_to_one_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# analyze sequence composition of design
def validate_design_sequence(sequence, num_clashes, advanced_settings):
    note_array = []

    # Check if protein contains clashes after relaxation
    if num_clashes > 0:
        note_array.append('Relaxed structure contains clashes.')

    # Check if the sequence contains disallowed amino acids
    if advanced_settings["omit_AAs"]:
        restricted_AAs = advanced_settings["omit_AAs"].split(',')
        for restricted_AA in restricted_AAs:
            if restricted_AA in sequence:
                note_array.append('Contains: '+restricted_AA+'!')

    # Analyze the protein
    analysis = ProteinAnalysis(sequence)

    # Calculate the reduced extinction coefficient per 1% solution
    extinction_coefficient_reduced = analysis.molar_extinction_coefficient()[0]
    molecular_weight = apply_advanced_rounding(analysis.molecular_weight() / 1000, advanced_settings)
    extinction_coefficient_reduced_1 = apply_advanced_rounding(extinction_coefficient_reduced / molecular_weight * 0.01, advanced_settings)

    # Check if the absorption is high enough
    if extinction_coefficient_reduced_1 <= 2:
        note_array.append(f'Absorption value is {extinction_coefficient_reduced_1}, consider adding tryptophane to design.')

    # Join the notes into a single string
    notes = ' '.join(note_array)

    return notes

# temporary function, calculate RMSD of input PDB and trajectory target
def target_pdb_rmsd(trajectory_pdb, starting_pdb, chain_ids_string, advanced_settings):
    # Parse the PDB files
    parser = PDBParser(QUIET=True)
    structure_trajectory = parser.get_structure('trajectory', trajectory_pdb)
    structure_starting = parser.get_structure('starting', starting_pdb)
    
    # Extract chain A from trajectory_pdb
    chain_trajectory = structure_trajectory[0]['A']
    
    # Extract the specified chains from starting_pdb
    chain_ids = chain_ids_string.split(',')
    residues_starting = []
    for chain_id in chain_ids:
        chain_id = chain_id.strip()
        chain = structure_starting[0][chain_id]
        for residue in chain:
            if is_aa(residue, standard=True):
                residues_starting.append(residue)
    
    # Extract residues from chain A in trajectory_pdb
    residues_trajectory = [residue for residue in chain_trajectory if is_aa(residue, standard=True)]
    
    # Ensure that both structures have the same number of residues
    min_length = min(len(residues_starting), len(residues_trajectory))
    residues_starting = residues_starting[:min_length]
    residues_trajectory = residues_trajectory[:min_length]
    
    # Collect CA atoms from the two sets of residues
    atoms_starting = [residue['CA'] for residue in residues_starting if 'CA' in residue]
    atoms_trajectory = [residue['CA'] for residue in residues_trajectory if 'CA' in residue]
    
    # Calculate RMSD using structural alignment
    sup = Superimposer()
    sup.set_atoms(atoms_starting, atoms_trajectory)
    rmsd = sup.rms
    
    return apply_advanced_rounding(rmsd, advanced_settings)

# detect C alpha clashes for deformed trajectories
def calculate_clash_score(pdb_file, threshold=2.4, only_ca=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    atoms = []
    atom_info = []  # Detailed atom info for debugging and processing

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'H':  # Skip hydrogen atoms
                        continue
                    if only_ca and atom.get_name() != 'CA':
                        continue
                    atoms.append(atom.coord)
                    atom_info.append((chain.id, residue.id[1], atom.get_name(), atom.coord))

    tree = cKDTree(atoms)
    pairs = tree.query_pairs(threshold)

    valid_pairs = set()
    for (i, j) in pairs:
        chain_i, res_i, name_i, coord_i = atom_info[i]
        chain_j, res_j, name_j, coord_j = atom_info[j]

        # Exclude clashes within the same residue
        if chain_i == chain_j and res_i == res_j:
            continue

        # Exclude directly sequential residues in the same chain for all atoms
        if chain_i == chain_j and abs(res_i - res_j) == 1:
            continue

        # If calculating sidechain clashes, only consider clashes between different chains
        if not only_ca and chain_i == chain_j:
            continue

        valid_pairs.add((i, j))

    return len(valid_pairs)

def _residues(chain):
    return [r for r in Selection.unfold_entities(chain, "R") if is_aa(r, standard=True)]

def _ca_for_resnums(chain, resnums_set):
    ca = []
    for r in _residues(chain):
        if r.id[1] in resnums_set and "CA" in r:
            ca.append(r["CA"])
    return ca

def _chain_ca_atoms(structure, chain_id):
    out = []
    for a in Selection.unfold_entities(structure, "A"):
        if a.get_parent().get_parent().id == chain_id and a.id == "CA":
            out.append(a)
    return out

def require_float(cfg, key):
    try:
        return float(cfg[key])
    except KeyError:
        raise ValueError(f"Missing required setting: {key}") from None

def save_superposition(fixed_struct, moved_struct, out_path, save_full_structure=False):
    builder = StructureBuilder()
    builder.init_structure("superposition")

    builder.init_model(0)
    for ch in moved_struct[0]:
        builder.structure[0].add(ch.copy())

    if save_full_structure:
        builder.init_model(1)
        for ch in fixed_struct[0]:
            builder.structure[1].add(ch.copy())

    io = PDBIO()
    io.set_structure(builder.get_structure())
    io.save(str(out_path))

def calculate_clashes_after_superpostion_on_full_target(
    model_pdb_path: str,
    target_full_structure: Structure,
    target_full_pae: np.ndarray,
    target_crop_indices: List[int],     # PDB residue numbers (1-based)
    superimposed_filepath: Optional[str] = None,
    clash_radius: float = 3.0,             # Clash distance cutoff (Å)
    target_chain_id: str = "A",                         # Target chain ID
    binder_chain_id: str = "B",
    clash_plddt_threshold: float = 70.0,
    clash_pae_threshold: float = 15.0,
) -> Dict[str, float]:
    """
    Superimpose a cropped target structure onto the full target,
    then detect CA-based clashes between the binder and the
    full target outside the crop region.

    Parameters
    ----------
    model_pdb_path : str
        Path to PDB containing cropped target + binder.
    target_full_structure : Bio.PDB.Structure
        Full target structure.
    target_full_pae : np.ndarray
        Full-target PAE matrix (NxN, zero-based indexing).
    target_crop_indices : List[int]
        PDB residue numbers defining the crop region.
    target_chain_id : str
        Chain ID of the target.
    advanced_settings : dict
        Settings containing clash radius.
    clash_radius : float
        Distance cutoff for CA-CA clashes (Å).
    binder_chain_id : str
        Chain ID of the binder.

    Returns
    -------
    Dict[str, float]
        Dictionary containing clash counts and confidence metrics.
    """

    crop_resnums = set(target_crop_indices)

    parser = PDBParser(QUIET=True)
    target_crop_binder = parser.get_structure("protein", model_pdb_path)

    target_full_chain = target_full_structure[0][target_chain_id]
    target_crop_chain = target_crop_binder[0][target_chain_id]

    # 1) Superimpose crop target onto full target (CA atoms on crop resnums)
    fixed_ca = _ca_for_resnums(target_full_chain, crop_resnums)
    moving_ca = [r["CA"] for r in _residues(target_crop_chain) if "CA" in r]

    if len(fixed_ca) != len(moving_ca):
        raise ValueError(f"Mismatch in number of CA atoms for alignment: {model_pdb_path}")

    sup = Superimposer()
    sup.set_atoms(fixed=fixed_ca, moving=moving_ca)
    sup.apply(target_crop_binder.get_atoms())

    # Save combined PDB (optional but useful)
    if superimposed_filepath:
        Path(superimposed_filepath).parent.mkdir(parents=True, exist_ok=True)
        save_superposition(target_full_structure, target_crop_binder, superimposed_filepath)

    # 2) Build CA lists for full target excluding crop, and binder CA
    full_res = _residues(target_full_chain)

    crop_zero_idx = np.array([i for i, r in enumerate(full_res) if r.id[1] in crop_resnums], dtype=int)

    full_ex_crop_zero_idx = []
    full_ex_crop_ca_atoms = []
    for i, r in enumerate(full_res):
        if r.id[1] in crop_resnums:
            continue
        if "CA" in r:
            full_ex_crop_zero_idx.append(i)     # maps KD-tree index -> residue index in full_res
            full_ex_crop_ca_atoms.append(r["CA"])

    binder_ca_atoms = _chain_ca_atoms(target_crop_binder, binder_chain_id)

    coords_target = np.asarray([a.get_coord() for a in full_ex_crop_ca_atoms], float)
    coords_binder = np.asarray([a.get_coord() for a in binder_ca_atoms], float)

    if coords_target.size == 0:
        # The crop covers the entire target
        return {"num_superposition_clashes": 0, 
                "avg_clash_plddt": float("nan"), 
                "max_clash_plddt": float("nan"),
                "avg_clash_pae": float("nan"), 
                "min_clash_pae": float("nan"),
                "clash_pass": True}

    # 3) KD-tree query: any target CA within threshold of any binder CA
    tree = cKDTree(coords_target)
    neighbor_lists = tree.query_ball_point(coords_binder, r=clash_radius, workers=4)

    if not any(neighbor_lists):
        return {"num_superposition_clashes": 0, 
                "avg_clash_plddt": float("nan"), 
                "max_clash_plddt": float("nan"),
                "avg_clash_pae": float("nan"), 
                "min_clash_pae": float("nan"),
                "clash_pass": True}
    
    hit_idx = np.unique(np.concatenate([np.asarray(x, dtype=int) for x in neighbor_lists if x]))
    clash_res_indices = np.unique(np.asarray(full_ex_crop_zero_idx, dtype=int)[hit_idx])  # indices into full_res, zero-indexed

    # 4) Clash confidence metrics: pLDDT and PAE
    plddts = [full_res[i]["CA"].get_bfactor() for i in clash_res_indices if "CA" in full_res[i]]
    avg_plddt = float(np.mean(plddts)) if plddts else float("nan")
    max_plddt = float(np.max(plddts)) if plddts else float("nan")

    sub1 = target_full_pae[np.ix_(clash_res_indices, crop_zero_idx)]
    sub2 = target_full_pae[np.ix_(crop_zero_idx, clash_res_indices)].T
    avg_pae = float(((sub1 + sub2) * 0.5).mean())
    min_pae = float(((sub1 + sub2) * 0.5).min())

    # Pass if clashing residues are "not confident" i.e. low pLDDT OR high PAE
    clash_pass = bool(max_plddt <= clash_plddt_threshold or min_pae >= clash_pae_threshold)

    return {
        "num_superposition_clashes": int(clash_res_indices.size),
        "avg_clash_plddt": avg_plddt,
        "max_clash_plddt": max_plddt,
        "avg_clash_pae": avg_pae,
        "min_clash_pae": min_pae,
        "clash_pass": clash_pass,
    }

def superimpose_on_full_target_and_save_pdb(best_model_pdb, target_full_structure, target_crop_indices, target_settings, mpnn_design_name, best_model_number, design_paths):
    crop_resnums = set(target_crop_indices)

    parser = PDBParser(QUIET=True)
    target_crop_binder = parser.get_structure("protein", best_model_pdb)

    target_full_chain = target_full_structure[0][target_settings["chains"]]
    target_crop_chain = target_crop_binder[0][target_settings["chains"]]

    fixed_ca = _ca_for_resnums(target_full_chain, crop_resnums)
    moving_ca = [r["CA"] for r in _residues(target_crop_chain) if "CA" in r]

    if len(fixed_ca) != len(moving_ca):
        raise ValueError(f"Mismatch in number of CA atoms for alignment: {best_model_pdb}")

    sup = Superimposer()
    sup.set_atoms(fixed=fixed_ca, moving=moving_ca)
    sup.apply(target_crop_binder.get_atoms())

    superimposed_filepath = os.path.join(design_paths["Accepted/Superimposed"], f"{mpnn_design_name}_model{best_model_number}.pdb")
    
    Path(superimposed_filepath).parent.mkdir(parents=True, exist_ok=True)
    save_superposition(target_full_structure, target_crop_binder, superimposed_filepath)

# identify interacting residues at the binder interface
def hotspot_residues(trajectory_pdb, binder_chain="B", atom_distance_cutoff=4.0):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", trajectory_pdb)

    # Get the specified chain
    binder_atoms = Selection.unfold_entities(structure[0][binder_chain], 'A')
    binder_coords = np.array([atom.coord for atom in binder_atoms])

    # Get atoms and coords for the target chain
    target_atoms = Selection.unfold_entities(structure[0]['A'], 'A')
    target_coords = np.array([atom.coord for atom in target_atoms])

    # Build KD trees for both chains
    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    # Prepare to collect interacting residues
    interacting_residues = {}

    # Query the tree for pairs of atoms within the distance cutoff
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

    # Process each binder atom's interactions
    for binder_idx, close_indices in enumerate(pairs):
        binder_residue = binder_atoms[binder_idx].get_parent()
        binder_resname = binder_residue.get_resname()

        # Convert three-letter code to single-letter code using the manual dictionary
        if binder_resname in three_to_one_map:
            aa_single_letter = three_to_one_map[binder_resname]
            for close_idx in close_indices:
                target_residue = target_atoms[close_idx].get_parent()
                interacting_residues[binder_residue.id[1]] = aa_single_letter

    return interacting_residues

# calculate secondary structure percentage of design
def calc_ss_percentage(pdb_file, advanced_settings, chain_id="B", atom_distance_cutoff=4.0):
    # Parse the structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    model = structure[0]  # Consider only the first model in the structure

    # Calculate DSSP for the model
    dssp = DSSP(model, pdb_file, dssp=advanced_settings["dssp_path"])

    # Prepare to count residues
    ss_counts = defaultdict(int)
    ss_interface_counts = defaultdict(int)
    plddts_interface = []
    plddts_ss = []

    # Get chain and interacting residues once
    chain = model[chain_id]
    interacting_residues = set(hotspot_residues(pdb_file, chain_id, atom_distance_cutoff).keys())

    for residue in chain:
        residue_id = residue.id[1]
        if (chain_id, residue_id) in dssp:
            ss = dssp[(chain_id, residue_id)][2]  # Get the secondary structure
            ss_type = 'loop'
            if ss in ['H', 'G', 'I']:
                ss_type = 'helix'
            elif ss == 'E':
                ss_type = 'sheet'

            ss_counts[ss_type] += 1

            if ss_type != 'loop':
                # calculate secondary structure normalised pLDDT
                avg_plddt_ss = sum(atom.bfactor for atom in residue) / len(residue)
                plddts_ss.append(avg_plddt_ss)

            if residue_id in interacting_residues:
                ss_interface_counts[ss_type] += 1

                # calculate interface pLDDT
                avg_plddt_residue = sum(atom.bfactor for atom in residue) / len(residue)
                plddts_interface.append(avg_plddt_residue)

    # Calculate percentages
    total_residues = sum(ss_counts.values())
    total_interface_residues = sum(ss_interface_counts.values())

    percentages = calculate_percentages(total_residues, ss_counts['helix'], ss_counts['sheet'], advanced_settings)
    interface_percentages = calculate_percentages(total_interface_residues, ss_interface_counts['helix'], ss_interface_counts['sheet'], advanced_settings)

    i_plddt = apply_advanced_rounding(sum(plddts_interface) / len(plddts_interface) / 100, advanced_settings) if plddts_interface else 0
    ss_plddt = apply_advanced_rounding(sum(plddts_ss) / len(plddts_ss) / 100, advanced_settings) if plddts_ss else 0

    return (*percentages, *interface_percentages, i_plddt, ss_plddt)

def calculate_percentages(total, helix, sheet, advanced_settings):
    helix_percentage = apply_advanced_rounding((helix / total) * 100, advanced_settings) if total > 0 else 0
    sheet_percentage = apply_advanced_rounding((sheet / total) * 100, advanced_settings) if total > 0 else 0
    loop_percentage = apply_advanced_rounding(((total - helix - sheet) / total) * 100, advanced_settings) if total > 0 else 0

    return helix_percentage, sheet_percentage, loop_percentage