#!/usr/bin/env -S bash -l
#SBATCH --job-name=run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=12:00:00
#SBATCH --output=/work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/05_run_filter_domains_and_create_sweep_index.out



conda activate cluster_complement_4f
python /work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/05_filter_domains_and_create_sweep_index.py

echo "Done."