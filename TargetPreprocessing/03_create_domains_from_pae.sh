# pip install python-igraph
#
# Usage:
#  ./02_create_domains_from_pae.sh [INPUT_PAE_DIR] [OUTPUT_DOMAIN_DIR]
#
# Arguments (optional):
#  INPUT_PAE_DIR       - directory containing *.pae.json (defaults to $DEFAULT_PAE_INPUT_DIR)
#  OUTPUT_DOMAIN_DIR   - directory to write *.domains.csv (defaults to $DEFAULT_DOMAIN_OUTPUT_DIR)
#
# Environment-overridable variables (can export before running):
#  PYTHON         - path to the python executable (default set below)
#  SCRIPT         - path to pae_to_domains.py (default set below)
#  RESOLUTION     - float resolution passed to pae_to_domains (default 1.0)
#  INPUT_PAE_DIR  - override input dir without positional args
#  OUTPUT_DOMAIN_DIR - override output dir without positional args
#
# If INPUT_PAE_DIR is not provided as the first positional arg, a sensible default
# is used. OUTPUT_DOMAIN_DIR defaults to the default domain output dir.

set -euo pipefail

# Defaults (override by exporting env vars or passing args)
SCRIPT_DIR=$(dirname "$0")
PYTHON="${PYTHON:-/work/upcorreia/users/wenckste/miniconda3/envs/BindCraft/bin/python}"
SCRIPT="${SCRIPT:-${SCRIPT_DIR}/pae_to_domains/pae_to_domains.py}"
DEFAULT_PAE_INPUT_DIR=${SCRIPT_DIR}/../Data/AFDB_v6/PAE
DEFAULT_DOMAIN_OUTPUT_DIR=${SCRIPT_DIR}/../Data/AFDB_v6/CSV_Domains

RESOLUTION="${RESOLUTION:-1.0}"
PAE_CUTOFF="${PAE_CUTOFF:-15.0}"

INPUT_PAE_DIR="${INPUT_PAE_DIR:-$DEFAULT_PAE_INPUT_DIR}"
OUTPUT_DOMAIN_DIR="${OUTPUT_DOMAIN_DIR:-$DEFAULT_DOMAIN_OUTPUT_DIR}"

# Positional overrides (if provided)
if [ $# -ge 1 ] && [ -n "${1:-}" ]; then
	INPUT_PAE_DIR="$1"
fi
if [ $# -ge 2 ] && [ -n "${2:-}" ]; then
	OUTPUT_DOMAIN_DIR="$2"
fi


echo "============================================================"
echo "        create_domains_from_pae.sh configuration"
echo "============================================================"
printf "%-22s : %s\n" "Python executable" "$PYTHON"
printf "%-22s : %s\n" "Python script" "$SCRIPT"
printf "%-22s : %s\n" "Input PAE directory" "$INPUT_PAE_DIR"
printf "%-22s : %s\n" "Output domain dir" "$OUTPUT_DOMAIN_DIR"
printf "%-22s : %s\n" "Resolution" "$RESOLUTION"
printf "%-22s : %s\n" "PAE cutoff" "$PAE_CUTOFF"
printf "%-22s : %s\n" "Script directory" "$SCRIPT_DIR"
echo "============================================================"
echo


mkdir -p "$OUTPUT_DOMAIN_DIR"


shopt -s nullglob
mapfile -d '' files < <(find "$INPUT_PAE_DIR" -maxdepth 1 -type f -name "*.pae.json" -print0)
if [ ${#files[@]} -eq 0 ]; then
	echo "No .pae.json files found in '$INPUT_PAE_DIR'" >&2
	exit 1
fi

total=${#files[@]}
idx=0
for f in "${files[@]}"; do
	idx=$((idx + 1))
	base=$(basename -- "$f")
	# remove the trailing .pae.json (exact match)
	name="${base%%.pae.json}"
	out="${OUTPUT_DOMAIN_DIR}/${name}.domains.csv"
	echo "[$idx/$total] Processing: $f -> $out"
	"$PYTHON" "$SCRIPT" "$f" --output_file "$out" --resolution "$RESOLUTION" --library igraph --pae_cutoff "$PAE_CUTOFF"
done

echo "All done. Outputs written to: $OUTPUT_DOMAIN_DIR"