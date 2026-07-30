#!/bin/bash

# Defaults (can be overridden by CLI args).
BINDCRAFT_DIR_DEFAULT="/mnt/lpdi/scratch/wenckste/code/binderdesign/BindCraft"
PYTHON_DIR_DEFAULT="/home/wenckste/bin/miniconda3/envs/BindCraft/bin/python"

# Parse external arguments.
# Usage: ./bindcraft.sh [--bindcraft-dir PATH] [--python-dir PATH]
BINDCRAFT_DIR="${BINDCRAFT_DIR_DEFAULT}"
PYTHON_DIR="${PYTHON_DIR_DEFAULT}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bindcraft-dir)
      BINDCRAFT_DIR="$2"
      shift 2
      ;;
    --python-dir)
      PYTHON_DIR="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 [--bindcraft-dir PATH] [--python-dir PATH]"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      echo "Usage: $0 [--bindcraft-dir PATH] [--python-dir PATH]"
      exit 1
      ;;
  esac
done

# Create a run directory with timestamp
SCRIPT_PATH="$(realpath "$0")"
TIMESTAMP=$(date +%Y%m%d_%H%M%S%3N)
RUN_PATH="$(dirname "$SCRIPT_PATH")/runs/$TIMESTAMP"
mkdir -p "$RUN_PATH"

# Modify the run didrectory in the json file and copy it to the run directory
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
TARGET_NAME="$(basename "$SCRIPT_DIR")"
cp "$SCRIPT_DIR/target/${TARGET_NAME}.json" "$RUN_PATH/${TARGET_NAME}.json"
JSON="$RUN_PATH/${TARGET_NAME}.json"
sed -i "s|__designpath__|${RUN_PATH}|" "$JSON"

# Run BindCraft with provided paths.
"$PYTHON_DIR" "$BINDCRAFT_DIR/bindcraft.py" \
  --settings "$RUN_PATH/${TARGET_NAME}.json" \
  --filters "$BINDCRAFT_DIR/settings_filters/default_filters.json" \
  --advanced "$BINDCRAFT_DIR/settings_advanced/default_4stage_multimer.json" \
  2>&1 | tee "$RUN_PATH/bindcraft_${TARGET_NAME}.log"

echo "Writing to $RUN_PATH/bindcraft_${TARGET_NAME}.log"

