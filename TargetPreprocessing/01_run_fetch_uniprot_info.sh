#!/usr/bin/env -S bash -l
#SBATCH --job-name=run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=12:00:00
#SBATCH --output=/work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/01_run_fetch_uniprot_info.out



conda activate BindCraft
python /work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/01_fetch_uniprot_info.py

