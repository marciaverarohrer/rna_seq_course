#!/usr/bin/env bash

#SBATCH --mail-type=end
#SBATCH --mail-user=marcia.rohrer@students.unibe.ch
#SBATCH --output=/data/users/mrohrer/rna_seq_course/results/slurm_output/slurm_hisat2_mapping/slurm-alignement-mapping-%j.out
#SBATCH --job-name=rna_seq_alignment_mapping
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8

CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

FILE_1=$1
FILE_2=$2
SAMPLE_NAME=$3
INDEXED_DIR=$4
OUT_DIR=$5
SPLICE_SITE_FILE=$6

# Change the host paths to container paths with the /mnt added
# !!! for indexed directory hisat2 needs the filename also not just the directory!
# !!! /Mus_musculus.GRCm39 added to the indexed directory for inside the container.

FILE_1_IN_CONTAINER=/mnt/${FILE_1#*/data/users/mrohrer/}
FILE_2_IN_CONTAINER=/mnt/${FILE_2#*/data/users/mrohrer/}
INDEXED_DIR_IN_CONTAINER=/mnt/${INDEXED_DIR#*/data/users/mrohrer/}/Mus_musculus.GRCm39
OUT_DIR_IN_CONTAINER=/mnt/${OUT_DIR#*/data/users/mrohrer/}
SPLICE_SITE_FILE_IN_CONTAINER=/mnt/${SPLICE_SITE_FILE#*/data/users/mrohrer/}

#to make the filepaths accessible inside the container, /mnt was added to distinguish.
#above all filepaths were adjusted and /data/users/mrohrer was replaced with /mnt

apptainer exec -B /data/users/mrohrer/:/mnt $CONTAINER hisat2 -p ${SLURM_CPUS_PER_TASK} \
-x $INDEXED_DIR_IN_CONTAINER \
-1 $FILE_1_IN_CONTAINER -2 $FILE_2_IN_CONTAINER \
-S ${OUT_DIR_IN_CONTAINER}/${SAMPLE_NAME}.sam \
--summary-file ${OUT_DIR_IN_CONTAINER}/summary/${SAMPLE_NAME}_summary.txt \
--rna-strandness RF \
--known-splicesite-infile $SPLICE_SITE_FILE_IN_CONTAINER

