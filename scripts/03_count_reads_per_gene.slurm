#!/usr/bin/env bash

#SBATCH --job-name=featureCounts
#SBATCH --output=/data/users/mrohrer/rna_seq_course/results/slurm_output/slurm_counts/slurm_featurecounts_%j.out
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --mail-type=end
#SBATCH --mail-user=marcia.rohrer@students.unibe.ch

BAM_DIRECTORY=$1
REFERENCE_GTF=$2
OUTPUT_FILE=$3

CONTAINER=/containers/apptainer/subread_2.0.6.sif

# again bind the directory so the container can access the files
apptainer exec --bind /data/users/mrohrer:/data/users/mrohrer $CONTAINER \
featureCounts \
  -p --countReadPairs \
  -t exon \
  -s 2 \
  -T ${SLURM_CPUS_PER_TASK} \
  -a $REFERENCE_GTF \
  -o $OUTPUT_FILE \
  $BAM_DIRECTORY/*.bam

# flags explained:
# -p: paired-end data
# --countReadPairs: counts the fragments
# -t exon: this counts only genes that overlap an exon region
# -g gene_id: indicates the count per gene, this should be default.
# -s 2: 2 means reverse stranded library
# -T: number of threads
# -a: GTF annotation of reference genome
# -o: output file
# input BAMs: all files in the directory



