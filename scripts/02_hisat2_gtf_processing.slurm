#!/usr/bin/env bash

#SBATCH --mail-type=end,error
#SBATCH --mail-user=marcia.rohrer@students.unibe.ch
#SBATCH --output=/data/users/mrohrer/rna_seq_course/results/slurm_output/slurm_hisat2_gtf_processing/slurm-alignement-gtf-processing-%j.out
#SBATCH --partition=pibu_el8
#SBATCH --job-name=rna_seq_alignment_gtf_processing
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --time=01:00:00

#important again to have the unzipped gtf annotation (.gtf , not .gtf.gz)!

CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

GTF_FILE=/data/users/mrohrer/rna_seq_course/data/Mus_musculus.GRCm39.115.gtf
OUTPUT_DIRECTORY=/data/users/mrohrer/rna_seq_course/results/alignment

# this will give information about the splice sites and exons from the gtf file
apptainer exec -B /data/users/mrohrer $CONTAINER hisat2_extract_splice_sites.py $GTF_FILE > $OUTPUT_DIRECTORY/splice_sites.txt
apptainer exec -B /data/users/mrohrer $CONTAINER hisat2_extract_exons.py  $GTF_FILE > $OUTPUT_DIRECTORY/exons.txt