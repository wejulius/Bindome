#!/bin/bash

############################################
# Configuration (user variables)
############################################
sweep_name='Human_Sweep_batch2-14'
number_of_final_designs=50

############################################
# Resolve script directory
############################################
SCRIPT_PATH="$(realpath "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

############################################
# Defaults (relative to this script; override via CLI)
############################################
sweep_dir_default="${SCRIPT_DIR}/../Sweeps"
domain_index_file_default="${SCRIPT_DIR}/../Data/Filtered/sweep_index_batch2.tsv"

sweep_dir="${sweep_dir_default}"
domain_index_file="${domain_index_file_default}"

# External options:
#   -s <sweep_dir>
#   -i <domain_index_file>
#   -n <sweep_name>
#   -f <number_of_final_designs>
while getopts ":s:i:n:f:" opt; do
  case "${opt}" in
    s) sweep_dir="${OPTARG}" ;;
    i) domain_index_file="${OPTARG}" ;;
    n) sweep_name="${OPTARG}" ;;
    f) number_of_final_designs="${OPTARG}" ;;
    *) echo "Usage: $0 [-s sweep_dir] [-i domain_index_file] [-n sweep_name] [-f number_of_final_designs]" >&2; exit 1 ;;
  esac
done

############################################
# Creation of sweep directories
############################################
mkdir -p "${sweep_dir}/${sweep_name}/targets"
cp "${domain_index_file}" "${sweep_dir}/${sweep_name}/sweep_index.tsv"

# (TSV parsing: skip comments/empty lines, preserve file order)
while IFS=$'\t' read -r uniprot_name target_name domain_pdb_path full_pdb_path json_path length design_time design_runs num_accepted_designs num_trajectories ; do
  # skip comments (starting with #) and empty lines
  [[ -z "${target_name}" || "${target_name}" == \#* ]] && continue

  echo "Creating folder for target: ${target_name} with length ${length}"

  prefix="${target_name:0:2}"
  mkdir -p "${sweep_dir}/${sweep_name}/targets/${prefix}"
  target_folder="${sweep_dir}/${sweep_name}/targets/${prefix}/${target_name}"
  cp -r "${SCRIPT_DIR}/Template" "${target_folder}" 

  # Rename the template files
  mv "${target_folder}/target/__name__.json" "${target_folder}/target/${target_name}.json"

  # Local target paths (sweep-local copies)
  target_pdb_local="${target_folder}/target/${target_name}.pdb"
  target_full_local="${target_folder}/target/${uniprot_name}.pdb"
  target_pae_local="${target_folder}/target/${uniprot_name}.pae.json"

  # Copy the all files to the target folder (consistent order)
  cp "${domain_pdb_path}" "${target_pdb_local}"
  cp "${full_pdb_path}" "${target_full_local}"
  cp "${json_path}" "${target_pae_local}"

  # within the json file, replace the placeholder with the actual target name
  # __name__ -> target_name
  # __pdbfilepath__ -> target_pdb_local
  # __number_of_final_designs__ --> number_of_final_designs
  # __fulltargetstructure__ --> target_full_local
  # __fulltargetpaejson__ --> target_pae_local
  sed -i "s/__name__/${target_name}/g" "${target_folder}/target/${target_name}.json"
  sed -i "s|__pdbfilepath__|${target_pdb_local}|g" "${target_folder}/target/${target_name}.json"
  sed -i "s|__fulltargetstructure__|${target_full_local}|g" "${target_folder}/target/${target_name}.json"
  sed -i "s|__fulltargetpaejson__|${target_pae_local}|g" "${target_folder}/target/${target_name}.json"
  sed -i "s/\"__number_of_final_designs__\"/${number_of_final_designs}/g" \
      "${target_folder}/target/${target_name}.json"

done < "${domain_index_file}"

echo "Sweep directories created successfully at: ${sweep_dir}/${sweep_name}"