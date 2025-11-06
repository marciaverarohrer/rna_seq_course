#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/mrohrer/rna_seq_course/results/slurm_output/fastqc_%j.out
#SBATCH --partition=pibu_el8
#SBATCH --error=/data/users/mrohrer/rna_seq_course/results/slurm_output/fastqc_%j.err
#SBATCH --mail-user=marcia.rohrer@students.unibe.ch
#SBATCH --mail-type=end,error

# Load the Apptainer module to use it
module load apptainer

# Run FastQC
apptainer exec $1 fastqc -t 1 -o $2 $3