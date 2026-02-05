#!/bin/bash -l

# output and error
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
# Job Name:
#SBATCH -D /home/light/analysispy_initial_qc/analysispy
#SBATCH -J analysispy_initial_qc
#SBATCH --nodes=1
#SBATCH -c 64
#SBATCH --time=48:00:00
#SBATCH --mem=100000


conda activate local_analysispy

python ~/local_analysis_initial_qc/il_ancient_snakemake_base/local_analysis_initial_qc/local_analysis_initial_qc_main.py -p test_local.json
