# Sweep Scripts

Turns a `sweep_index.tsv` (from `../TargetPreprocessing`) into a running BindCraft slurm array sweep; will require some adaptation based on your available HPC cluster. 

1. **`01_create_sweep_dirs.sh -i <index.tsv> -n <sweep_name> [-s sweep_dir] [-f n_final_designs]`**
   Creates `<sweep_dir>/<sweep_name>/targets/<XX>/<target>/` per row (from `Template/`).

2. **`02a_create_sweep_auto_submission_script.sh -d <sweep_dir> [-a/-A min/max_accepted] [-x/-X min/max_design_time] [-j/-t/-N/-g/-c/-p/-m slurm opts]`**
   (or `02b_...` — no auto-resubmit, uses a timeout instead). Builds + submits a Slurm array job running `bindcraft.sh` per target. 02a self-resubmits near walltime (`-r` to disable).

3. **`03_update_sweep_index_and_collect_designs.sh -d <sweep_dir> [-b 1]`**
   Refreshes design-progress columns in `sweep_index.tsv` and copies accepted designs into each target's `designs/`. Rerun anytime to check progress. `-d` accepting `{i}` loops batches 1-14.

Rerun 2 and 3 as needed while the sweep is running.
