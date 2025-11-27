#!/usr/bin/env bash

BAM_DIRECTORY=/data/users/mrohrer/rna_seq_course/results/alignment/hisat2_mapping/bam
REFERENCE_GTF=/data/users/mrohrer/rna_seq_course/data/Mus_musculus.GRCm39.115.gtf
OUTPUT_FILE=/data/users/mrohrer/rna_seq_course/results/03_counts/featureCounts_counts.txt

sbatch 03_count_reads_per_gene.slurm $BAM_DIRECTORY $REFERENCE_GTF $OUTPUT_FILE