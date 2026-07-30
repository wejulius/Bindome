#!/usr/bin/env bash
#SBATCH --job-name=run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=12:00:00
#SBATCH --output=/work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/03_run_create_domains_from_pae.out


bash /work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/03_create_domains_from_pae.sh

