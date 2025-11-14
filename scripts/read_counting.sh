#!/usr/bin/env bash

#this script should effectively caclulate the read counts for each sample.

#four lines describe one read, where the third line always start with an @
# by doing grep "@" we should obtain the number of reads


OUTPUT_DIRECTORY=/data/users/mrohrer/rna_seq_course/results
READS_DIRECTORY=/data/users/mrohrer/rna_seq_course/fastqc/reads_Lung

cd $OUTPUT_DIRECTORY
touch readcount.txt

#for filename in `ls -1 $READS_DIRECTORY/*.fastq.gz`; do
#    read_count=zcat $filename | grep "@" | wc -l
#    echo "${filename%.fastq.gz} : $read_count" >> readcount.txt
#done

for filename in $READS_DIRECTORY/*.fastq.gz; do
    read_count=$(zcat "$filename" | grep -c "^@")
    echo "$(basename ${filename%.fastq.gz}) : $read_count" >> readcount.txt
done