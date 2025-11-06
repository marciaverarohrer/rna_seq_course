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

# directories for shortcut and to make the code more overseeable
WORKDIR=/data/users/mrohrer/rna_seq_course/fastqc
READS=/data/users/mrohrer/rna_seq_course/fastqc/reads_Lung
CONTAINER=/containers/apptainer/fastqc-0.12.1.sif
OUTDIR=$WORKDIR/results

#if directory is non-existent and to direct into the working directory
mkdir -p $OUTDIR
cd $WORKDIR

for file in `ls -1 $READS/*.fastq.gz`;
do  
    sbatch /data/users/mrohrer/rna_seq_course/scripts/01_slurm_script.sh $CONTAINER $OUTDIR $file;
done