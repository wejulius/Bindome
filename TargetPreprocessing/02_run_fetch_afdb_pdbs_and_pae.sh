#!/usr/bin/env -S bash -l
#SBATCH --job-name=run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=12:00:00
#SBATCH --output=/work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/02_run_fetch_afdb_pdbs_and_pae.out




conda activate cluster_complement_4f
python /work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/02_fetch_afdb_pdbs_and_pae.py

