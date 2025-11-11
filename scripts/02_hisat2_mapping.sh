#!/usr/bin/env bash

# filepaths
FASTQ_DIR=/data/users/mrohrer/rna_seq_course/fastqc/reads_Lung
OUT_DIR=/data/users/mrohrer/rna_seq_course/results/alignment/hisat2_mapping
INDEXED_DIR=/data/users/mrohrer/rna_seq_course/results/alignment/hisat2_index
SPLICE_SITE_FILE=/data/users/mrohrer/rna_seq_course/results/alignment/hisat2_gtf_processing/splice_sites.txt

# each fastq file from the data (1 and 2 of one sample are treated together) is mapped against the reference
for file_1 in `ls -1 $FASTQ_DIR/*_1.fastq.gz`; do
    file_2=${file_1%_1.fastq.gz}_2.fastq.gz
    sample_name=$(basename ${file_1%_1.fastq.gz})
    sbatch 02_hisat2_mapping.slurm $file_1 $file_2 $sample_name $INDEXED_DIR $OUT_DIR $SPLICE_SITE_FILE
done

# for this sbatch command to work the bash script has to be executed from inside the script directory!

# the command basename strips the path from the beginning of the string, where the other part strips the ending (_1fasta.gz )
# this line inside the for loop will create a sample name with only the sample title, without ending or directory.
