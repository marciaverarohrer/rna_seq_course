#!/usr/bin/env bash

# filepaths
SAM_DIR=/data/users/mrohrer/rna_seq_course/results/alignment/hisat2_mapping/sam
OUT_DIR=/data/users/mrohrer/rna_seq_course/results/alignment/hisat2_mapping/bam
REF_FILE=/data/users/mrohrer/rna_seq_course/data/Mus_musculus.GRCm39.dna.primary_assembly.fa

# each fastq file will be mapped to the reference genome
for file in `ls -1 $SAM_DIR/*.sam`; do
    sbatch 02_samtools.slurm $file $REF_FILE ${OUT_DIR}/$(basename ${file%.sam}).bam
done