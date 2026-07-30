#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Resolve script directory
###############################################################################
SCRIPT_PATH="$(realpath "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

###############################################################################
# Configuration (defaults; override via CLI)
###############################################################################
SWEEP_DIRECTORY_DEFAULT="${SCRIPT_DIR}/../Sweeps/Human_Sweep_batch14-14"
SWEEP_DOMAIN_INDEX_FILE_DEFAULT="${SWEEP_DIRECTORY_DEFAULT}/sweep_index.tsv"
TARGETS_DIR_DEFAULT="${SWEEP_DIRECTORY_DEFAULT}/targets"

SWEEP_DIRECTORY="${SWEEP_DIRECTORY_DEFAULT}"
SWEEP_DOMAIN_INDEX_FILE="${SWEEP_DOMAIN_INDEX_FILE_DEFAULT}"
TARGETS_DIR="${TARGETS_DIR_DEFAULT}"

BACKUP="${BACKUP:-0}"

while getopts ":d:i:t:b:" opt; do
  case "${opt}" in
    d) SWEEP_DIRECTORY="${OPTARG}" ;;
    i) SWEEP_DOMAIN_INDEX_FILE="${OPTARG}" ;;
    t) TARGETS_DIR="${OPTARG}" ;;
    b) BACKUP="${OPTARG}" ;;
    *) echo "Usage: $0 [-d sweep_dir] [-i sweep_index.tsv] [-t targets_dir] [-b backup_flag]" >&2; exit 1 ;;
  esac
done

# Expand template once, then run per batch
if [[ "$SWEEP_DIRECTORY" == *"{i}"* ]]; then
  for i in $(seq 1 14); do
    sweep_i="${SWEEP_DIRECTORY//\{i\}/$i}"
    index_i="${sweep_i}/sweep_index.tsv"
    targets_i="${sweep_i}/targets"

    [[ -d "$sweep_i" ]] || { echo "Skipping missing: $sweep_i"; continue; }
    bash "$0" -d "$sweep_i" -i "$index_i" -t "$targets_i" -b "$BACKUP"
  done
  exit 0
fi

###############################################################################
# Validation
###############################################################################
if [[ ! -f "$SWEEP_DOMAIN_INDEX_FILE" ]]; then
  echo "ERROR: index file not found: $SWEEP_DOMAIN_INDEX_FILE" >&2
  exit 1
fi
if [[ ! -d "$TARGETS_DIR" ]]; then
  echo "ERROR: targets dir not found: $TARGETS_DIR" >&2
  exit 1
fi

###############################################################################
# Step 1 — Update the sweep index with design metrics
###############################################################################
echo "Updating sweep index..."

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

awk -v targets_dir="$TARGETS_DIR" '
BEGIN { FS="\t"; OFS="\t" }

# Preserve comments and blank lines verbatim
/^#/ || NF==0 { print; next }

{
  target = $2
  prefix = substr(target, 1, 2)
  rdir = targets_dir "/" prefix "/" target "/runs"

  design_runs = 0
  design_time = 0.0
  accepted    = 0
  n_trajectories = 0

  # 1) Enumerate run directories under runs/
  cmd = "ls -1d " rdir "/* 2>/dev/null"
  while ((cmd | getline path) > 0) {
    sub(/\n$/, "", path)  # strip trailing newline from ls output

    # 2) Count only directories
    testdir = sprintf("test -d \"%s\"", path)
    if (system(testdir) != 0) continue

    # 3) If runtime_log.csv exists, read its last data row
    log_fp = sprintf("%s/runtime_log.csv", path)
    testlog_cmd = sprintf("test -f \"%s\"", log_fp)
    if (system(testlog_cmd) != 0) continue

    design_runs++

    parse_cmd = sprintf("awk -F, '\''NR>1 && $1!=\"\" { last=$0 } END{ if(last!=\"\") print last }'\'' \"%s\"", log_fp)
    lastline = ""
    parse_cmd | getline lastline
    close(parse_cmd)

    # 4) Accumulate time, trajectories, and accepted designs (cols 1,3,4)
    if (lastline != "") {
      n = split(lastline, f, ",")
      if (n >= 4) {
        design_time    += (f[1] + 0)
        n_trajectories += (f[3] + 0)
        accepted       += (f[4] + 0)
      }
    }
  }
  close(cmd)

  # Write: first 6 columns as-is, then updated metrics
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\n", \
    $1, $2, $3, $4, $5, $6, design_time, design_runs, accepted, n_trajectories
}
' "$SWEEP_DOMAIN_INDEX_FILE" > "$tmp"

if [[ "$BACKUP" == "1" ]]; then
  cp -a "$SWEEP_DOMAIN_INDEX_FILE" "${SWEEP_DOMAIN_INDEX_FILE}.bak.$(date +%Y%m%d_%H%M%S)"
fi

mv "$tmp" "$SWEEP_DOMAIN_INDEX_FILE"
trap - EXIT

echo "Updated: $SWEEP_DOMAIN_INDEX_FILE"

###############################################################################
# Step 2 — Collect Accepted PDBs (parallelized per target)
###############################################################################
echo "Collecting accepted designs..."

command -v rsync >/dev/null 2>&1 || {
  echo "ERROR: rsync not found in PATH" >&2
  exit 1
}

command -v xargs >/dev/null 2>&1 || {
  echo "ERROR: xargs not found in PATH" >&2
  exit 1
}

process_target() {
  local target_dir="$1"
  local designs_dir="$target_dir/designs"

  echo "Processing target: $(basename "$target_dir")"

  [[ -d "$target_dir/runs" ]] || return 0
  mkdir -p "$designs_dir"

  # start fresh on each script run (avoid duplicate appends across reruns)
  rm -f \
    "$designs_dir/final_design_stats_overall.csv" \
    "$designs_dir/failure_csv_overall.csv" \
    "$designs_dir/mpnn_design_stats_overall.csv" \
    "$designs_dir/trajectory_stats_overall.csv"

  for run_dir in "$target_dir"/runs/*; do
    [[ -d "$run_dir" ]] || continue

    acc_dir="$run_dir/Accepted"
    acc_superimposed_dir="$acc_dir/Superimposed"
    acc_pae_json_binder="$acc_dir/Pae_json/Binder"
    acc_pae_json_complex="$acc_dir/Pae_json/Complex"

    [[ -d "$acc_dir" ]] || continue

    # ------------------------------------------------------------
    # Copy accepted PDBs (fast, no overwrite)
    # ------------------------------------------------------------
    rsync -a --ignore-existing \
      --include='*.pdb' --exclude='*' \
      "$acc_dir"/ "$designs_dir"/

    # Copy superimposed PDBs
    if [[ -d "$acc_superimposed_dir" ]]; then
      mkdir -p "$designs_dir/Superimposed"
      rsync -a --ignore-existing \
        --include='*.pdb' --exclude='*' \
        "$acc_superimposed_dir"/ "$designs_dir/Superimposed"/
    fi

    # Copy PAE JSONs of binder
    if [[ -d "$acc_pae_json_binder" ]]; then
      mkdir -p "$designs_dir/Pae_json/Binder"
      rsync -a --ignore-existing \
        --include='*.json' --exclude='*' \
        "$acc_pae_json_binder"/ "$designs_dir/Pae_json/Binder"/
    fi

    # Copy PAE JSONs of complex
    if [[ -d "$acc_pae_json_complex" ]]; then
      mkdir -p "$designs_dir/Pae_json/Complex"
      rsync -a --ignore-existing \
        --include='*.json' --exclude='*' \
        "$acc_pae_json_complex"/ "$designs_dir/Pae_json/Complex"/
    fi

    # ------------------------------------------------------------
    # Merge CSVs (append without duplicating headers)
    # ------------------------------------------------------------

    merge_csv () {
      local src="$1"
      local dst="$2"
      if [[ -f "$src" ]]; then
        if [[ ! -f "$dst" ]]; then
          cp "$src" "$dst"
        else
          tail -n +2 "$src" >> "$dst"
        fi
      fi
    }

    merge_csv "$run_dir/final_design_stats.csv" \
              "$designs_dir/final_design_stats_overall.csv"

    merge_csv "$run_dir/failure_csv.csv" \
              "$designs_dir/failure_csv_overall.csv"

    merge_csv "$run_dir/mpnn_design_stats.csv" \
              "$designs_dir/mpnn_design_stats_overall.csv"

    merge_csv "$run_dir/trajectory_stats.csv" \
              "$designs_dir/trajectory_stats_overall.csv"

  done

}

export -f process_target

# Collect target directories
mapfile -t TARGET_DIRS < <(
  find "$TARGETS_DIR" -mindepth 2 -maxdepth 2 -type d | sort
)

# Number of parallel jobs (default 8; override with JOBS=...)
JOBS="${JOBS:-8}"

printf '%s\0' "${TARGET_DIRS[@]}" \
  | xargs -0 -n 1 -P "$JOBS" bash -c 'process_target "$@"' _

# Merge per-target final design CSVs into one sweep-level CSV
merge_sweep_csv() {
  local rel_name="$1"   # e.g. final_design_stats_overall.csv
  local out_csv="$2"
  local first=1

  : > "$out_csv"

  while IFS= read -r fp; do
    if [[ $first -eq 1 ]]; then
      cat "$fp" > "$out_csv"
      first=0
    else
      tail -n +2 "$fp" >> "$out_csv"
    fi
  done < <(find "$TARGETS_DIR" -type f -path "*/*/designs/${rel_name}" | sort)

  # remove empty output if no input files were found
  [[ -s "$out_csv" ]] || rm -f "$out_csv"
}

echo "Merging final design stats into sweep-level CSV..."
merge_sweep_csv "final_design_stats_overall.csv" \
  "${SWEEP_DIRECTORY}/final_design_stats_overall_sweep.csv"

echo "Accepted designs collected."