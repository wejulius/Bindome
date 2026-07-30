#!/usr/bin/env -S bash -l
#SBATCH --job-name=run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=72gb
#SBATCH --time=12:00:00
#SBATCH --output=/work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/04_segment_and_compute_domain_stats.out

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
lscpu | grep '^CPU(s)'

conda activate cluster_complement_4f
python /work/upcorreia/users/wenckste/Bindome/TargetPreprocessing/04_segment_and_compute_domain_stats.py

