#!/usr/bin/env bash

#SBATCH --mail-type=end
#SBATCH --mail-user=marcia.rohrer@students.unibe.ch
#SBATCH --output=/data/users/mrohrer/rna_seq_course/results/slurm_output/slurm_samtools/slurm-samtools-%j.out
#SBATCH --job-name=rna_seq_alignment_samtools
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8

set -e #so the jobs stop when facing an error

SAM_FILE=$1
REF_FILE=$2
OUT_FILE=$3
UNSORTED_FILE=${OUT_FILE%.bam}.unsorted.bam
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

# 'faidx' : Creates an index to the reference sequence in the FASTA format
# makes it quicker for samtools to map reads to positions during later steps
apptainer exec --bind /data/users/mrohrer:/data/users/mrohrer $CONTAINER samtools faidx $REF_FILE

# 'view' converts sam to bam formats
#apptainer exec $CONTAINER samtools view -hBs -@ ${SLURM_CPUS_PER_TASK} -t ${REF_FILE}.fai ${SAM_FILE} > ${UNSORTED_FILE}

apptainer exec --bind /data/users/mrohrer:/data/users/mrohrer \
  $CONTAINER \
  samtools view -@ ${SLURM_CPUS_PER_TASK} -b -h -o ${UNSORTED_FILE} ${SAM_FILE}

# 'sort' 
apptainer exec --bind /data/users/mrohrer:/data/users/mrohrer $CONTAINER samtools sort -@ ${SLURM_CPUS_PER_TASK} -o $OUT_FILE -m 4G ${UNSORTED_FILE}
rm ${UNSORTED_FILE}

#index bam (make the bai)
apptainer exec --bind /data/users/mrohrer:/data/users/mrohrer $CONTAINER samtools index -@ ${SLURM_CPUS_PER_TASK} $OUT_FILE