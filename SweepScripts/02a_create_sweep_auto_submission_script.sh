#!/usr/bin/env bash
set -euo pipefail

############################################
# Resolve script directory
############################################
SCRIPT_PATH="$(realpath "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

############################################
# Configuration (defaults; override via CLI)
############################################
SWEEP_DIRECTORY_DEFAULT="${SCRIPT_DIR}/../Sweeps/Human_Sweep_batch11-14"   # ADAPT THIS FOLDER
SWEEP_DOMAIN_INDEX_FILE_DEFAULT="${SWEEP_DIRECTORY_DEFAULT}/sweep_index.tsv"

SWEEP_DIRECTORY="${SWEEP_DIRECTORY_DEFAULT}"
SWEEP_DOMAIN_INDEX_FILE="${SWEEP_DOMAIN_INDEX_FILE_DEFAULT}"

SBATCH_JOB_NAME='sweep_batch11-14' # ADAPT THIS JOB NAME
SBATCH_TIME='12:00:00'
SBATCH_NODES='1'
SBATCH_GPUS_PER_NODE='1'
SBATCH_CPUS_PER_GPU='8'
SBATCH_PARTITION='boost_usr_prod'
#SBATCH_MEM='46gb'

PYTHON_PATH='/leonardo_work/AIFAC_S03_099/config/envs/BindCraft/bin/python'
BINDCRAFT_DIR="${SCRIPT_DIR}/../BindCraft"

# Filters for target selection (based on TSV columns)
MIN_ACCEPTED_DESIGNS=0
MAX_ACCEPTED_DESIGNS=48
MIN_DESIGN_TIME=0
MAX_DESIGN_TIME=6048000  # 10 weeks
AUTORESUBMIT=1

# External options:
# -d <sweep_directory>  -i <sweep_index.tsv>
# -j <job_name>         -t <time>          -N <nodes>
# -g <gpus_per_node>    -c <cpus_per_gpu>  -p <partition>
# -m <mem>
# -a <min_accepted>     -A <max_accepted>  -x <min_time> -X <max_time>
# -r (disable autoresubmit)
while getopts ":d:i:j:t:N:g:c:p:m:a:A:x:X:r" opt; do
  case "${opt}" in
    d) SWEEP_DIRECTORY="${OPTARG}" ;;
    i) SWEEP_DOMAIN_INDEX_FILE="${OPTARG}" ;;
    j) SBATCH_JOB_NAME="${OPTARG}" ;;
    t) SBATCH_TIME="${OPTARG}" ;;
    N) SBATCH_NODES="${OPTARG}" ;;
    g) SBATCH_GPUS_PER_NODE="${OPTARG}" ;;
    c) SBATCH_CPUS_PER_GPU="${OPTARG}" ;;
    p) SBATCH_PARTITION="${OPTARG}" ;;
    m) SBATCH_MEM="${OPTARG}" ;;
    a) MIN_ACCEPTED_DESIGNS="${OPTARG}" ;;
    A) MAX_ACCEPTED_DESIGNS="${OPTARG}" ;;
    x) MIN_DESIGN_TIME="${OPTARG}" ;;
    X) MAX_DESIGN_TIME="${OPTARG}" ;;
    r) AUTORESUBMIT=0 ;;
    *) echo "Usage: $0 [-d sweep_dir] [-i domain_index.tsv] [-j job] [-t time] [-N nodes] [-g gpus] [-c cpus] [-p part] [-m mem] [-a min_acc] [-A max_acc] [-x min_time] [-X max_time] [-r (disable autoresubmit)]" >&2; exit 1 ;;
  esac
done

############################################
# Build ordered target list (preserve TSV order)
# TSV columns assumed:
# uniprot_name  target_name  domain_pdb_path  full_pdb_path  json_path  length design_time design_runs num_accepted_designs num_trajectories
############################################
TARGET_LIST_FILE="${SWEEP_DIRECTORY}/targets.list"  # generated alongside this script
: > "${TARGET_LIST_FILE}"

while IFS=$'\t' read -r uniprot_name target_name domain_pdb_path full_pdb_path json_path length design_time design_runs num_accepted_designs num_trajectories; do
  [[ -z "${target_name}" || "${target_name}" == \#* ]] && continue
  # Coerce to integers (strip decimals)
  design_time_int="${design_time%%.*}"
  num_accepted_int="${num_accepted_designs%%.*}"
  [[ -z "${design_time_int}" ]] && design_time_int=0
  [[ -z "${num_accepted_int}" ]] && num_accepted_int=0

  # Apply filters
  if (( num_accepted_int < MIN_ACCEPTED_DESIGNS || num_accepted_int > MAX_ACCEPTED_DESIGNS )); then continue; fi
  if (( design_time_int < MIN_DESIGN_TIME || design_time_int > MAX_DESIGN_TIME )); then continue; fi
  # Remaining per-target budget (seconds)
  remaining_budget=$(( MAX_DESIGN_TIME - design_time_int ))
  if (( remaining_budget < 0 )); then remaining_budget=0; fi
  printf "%s\t%s\t%d\t%d\n" \
    "${uniprot_name}" "${target_name}" "${design_time_int}" "${remaining_budget}" \
  >> "${TARGET_LIST_FILE}"
done < "${SWEEP_DOMAIN_INDEX_FILE}"

NUM_TASKS="$(wc -l < "${TARGET_LIST_FILE}" | tr -d ' ')"
[[ "${NUM_TASKS}" -gt 0 ]] || { echo "ERROR: No targets found in ${SWEEP_DOMAIN_INDEX_FILE}" >&2; exit 1; }

# Slurm arrays are typically 0-based or 1-based; here we use 0-based
ARRAY_RANGE="0-$((NUM_TASKS-1))"

############################################
# Generate sbatch array script
############################################
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
SBATCH_SUBMIT_SCRIPT="02_run_sweep_array_${TIMESTAMP}.sh"
SBATCH_LOG_REL_PATH='logs/%A/02_run_sweep_%A_%a.log'

SBATCH_PATH="${SWEEP_DIRECTORY}/${SBATCH_SUBMIT_SCRIPT}"

SBATCH_SIGNAL_LINE=''
if [[ "${AUTORESUBMIT}" -eq 1 ]]; then
  SBATCH_SIGNAL_LINE='#SBATCH --signal=B:USR1@120'
fi

cat > "${SBATCH_PATH}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${SBATCH_JOB_NAME}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --nodes=${SBATCH_NODES}
#SBATCH --ntasks=1
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --qos=qos_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --array=${ARRAY_RANGE}
${SBATCH_SIGNAL_LINE}
#SBATCH --output=${SWEEP_DIRECTORY}/${SBATCH_LOG_REL_PATH}
#SBATCH --account=aifac_s03_099
#SBATCH --exclude=lrdn0980,lrdn2997,lrdn0614

set -euo pipefail
echo "Running on node: \$(hostname)"
#module load gcc/13.2.0
#module load cuda/12.4.1

SWEEP_DIRECTORY="${SWEEP_DIRECTORY}"
TARGET_LIST_FILE="${TARGET_LIST_FILE}"
AUTORESUBMIT="${AUTORESUBMIT}"
THIS_SBATCH_SCRIPT="${SBATCH_PATH}"

maybe_autoresubmit() {
  [[ "${AUTORESUBMIT}" == "1" ]] || return 0
  local reason="\${1:-end}"
  local array_log_dir="\${SWEEP_DIRECTORY}/logs/\${SLURM_ARRAY_JOB_ID}"
  local lock_file="\${array_log_dir}/.autoresubmit.lock"
  local max_active_task

  mkdir -p "\${array_log_dir}"

  (
    flock -x 200

    [[ -e "\${lock_file}.done" ]] && exit 0

    max_active_task="\$(
      squeue -h -r -j "\${SLURM_ARRAY_JOB_ID}" -t PD,R,CG,CF -o "%i" |
      awk -F'_' -v ajid="\${SLURM_ARRAY_JOB_ID}" 'NF==2 && \$1==ajid && \$2 ~ /^[0-9]+$/ { if (\$2+0 > max) max=\$2+0 } END { if (max == "") print -1; else print max }'
    )"

    if [[ "\${SLURM_ARRAY_TASK_ID}" -eq "\${max_active_task}" ]]; then
      echo "Autoresubmit (\${reason}) from highest active task \${SLURM_ARRAY_TASK_ID}"
      sbatch --dependency=afterany:\${SLURM_ARRAY_JOB_ID} "\${THIS_SBATCH_SCRIPT}"
      touch "\${lock_file}.done"
    else
      echo "Not resubmitting (\${reason}): task \${SLURM_ARRAY_TASK_ID} is not highest active task (max=\${max_active_task})"
    fi
  ) 200>"\${lock_file}"
}

on_usr1() {
  echo "[\$(date -Is)] Received USR1 (2 min before walltime)."
  maybe_autoresubmit "USR1"
}

if [[ "${AUTORESUBMIT}" == "1" ]]; then
  trap 'on_usr1' USR1
fi

# Get line (0-based index) for this array task
LINE="\$(sed -n "\$((SLURM_ARRAY_TASK_ID+1))p" "\${TARGET_LIST_FILE}")"
UNIPROT_NAME="\$(echo "\${LINE}" | cut -f1)"
TARGET_NAME="\$(echo "\${LINE}" | cut -f2)"

PREFIX="\${TARGET_NAME:0:2}"
TARGET_DIR="\${SWEEP_DIRECTORY}/targets/\${PREFIX}/\${TARGET_NAME}"

echo "[\$(date -Is)] SLURM_JOB_ID=\${SLURM_JOB_ID} TASK_ID=\${SLURM_ARRAY_TASK_ID}"
echo "Target: \${TARGET_NAME} (uniprot: \${UNIPROT_NAME})"
echo "Working dir: \${TARGET_DIR}"

if [[ ! -d "\${TARGET_DIR}" ]]; then
  echo "ERROR: Target directory not found: \${TARGET_DIR}" >&2
  exit 2
fi

if [[ ! -x "\${TARGET_DIR}/bindcraft.sh" ]]; then
  echo "ERROR: bindcraft.sh not found or not executable in: \${TARGET_DIR}" >&2
  echo "Tip: chmod +x \${TARGET_DIR}/bindcraft.sh" >&2
  exit 3
fi

cd "\${TARGET_DIR}"

echo "Slurm will enforce hard stop via USR1@T-120s"

# Run bindcraft in background and wait, so the shell can handle signals/traps cleanly.
./bindcraft.sh --bindcraft-dir "${BINDCRAFT_DIR}" --python-dir "${PYTHON_PATH}" &
main_pid=\$!

wait_rc=0
wait "\${main_pid}" || wait_rc=\$?
RC="\${wait_rc}"

if [[ "\${RC}" -eq 124 || "\${RC}" -eq 137 || "\${RC}" -eq 143 ]]; then
  echo "bindcraft terminated due to timeout/signal (rc=\${RC}). Treating as budget stop."
  exit 0
fi

echo "bindcraft completed for target \${TARGET_NAME}"
exit "\${RC}"
EOF

chmod +x "${SBATCH_PATH}"

############################################
# Submit
############################################
echo "Prepared target list: ${TARGET_LIST_FILE} (${NUM_TASKS} tasks)"
echo "Wrote sbatch script:  ${SBATCH_PATH}"
echo "Logs will go to:      ${SWEEP_DIRECTORY}/${SBATCH_LOG_REL_PATH}"
echo

sbatch "${SBATCH_PATH}"
echo "Submitted sbatch job array."