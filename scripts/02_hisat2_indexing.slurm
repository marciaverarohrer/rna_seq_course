#!/usr/bin/env bash

#SBATCH --job-name=hisat2_index
#SBATCH --mail-type=end,error
#SBATCH --mail-user=marcia.rohrer@students.unibe.ch
#SBATCH --output=/data/users/mrohrer/rna_seq_course/results/slurm_output/hisat2_indexing/alignment_index_%j.out
#SBATCH --error=/data/users/mrohrer/rna_seq_course/results/slurm_output/hisat2_indexing/alignment_index_%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --partition=pibu_el8

# loading the path of the container to use hisat2 (provided)
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

# it is important to unzip the fasta file to use it!
# Input and output paths
#REF_GENOME_FASTA=/data/users/mrohrer/rna_seq_course/data/mouse_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
REF_GENOME_FASTA_UNZIPPED=/data/users/mrohrer/rna_seq_course/data/Mus_musculus.GRCm39.dna.primary_assembly.fa
OUTPUT_DIRECTORY=/data/users/mrohrer/rna_seq_course/results/alignment/hisat2_index

# create the output directory
mkdir -p $OUTPUT_DIRECTORY

# to unzip the reference fasta genome (apparently Hisat2 cannot read it zipped)
# gunzip -c $REF_GENOME_FASTAGENOME_FASTA > $REF_GENOME_FASTA_UNZIPPED
# since this didn't work I unzipped it separately in the terminal and used the link above.

# Run HISAT2 inside the provided container
apptainer exec $CONTAINER hisat2-build $REF_GENOME_FASTA_UNZIPPED $OUTPUT_DIRECTORY/Mus_musculus.GRCm39 -p ${SLURM_CPUS_PER_TASK}
