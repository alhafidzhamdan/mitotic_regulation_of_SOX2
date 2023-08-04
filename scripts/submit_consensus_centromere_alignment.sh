#!/bin/bash

# To run this script, do 
# qsub -t n submit_consensus_centromere_alignment.sh

#$ -N consensus_centr
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -pe sharedmem 16
#$ -l h_rt=40:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh


CONFIG=$1
IDS=$2
READ_DIR=$3
REPEAT_ALIGMENT_DIR=$4
FASTA_LIST=$5

source $CONFIG

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1`

cd $READ_DIR

while read REF
do
    BASE=`basename $REF`
    SUFFIX=${BASE%.fa}
    if [[ ! -f $REPEAT_ALIGMENT_DIR/${SAMPLE_ID}.${SUFFIX}.sorted.bam.stats.out ]]; then
        echo "Running for $SUFFIX"
        bwa mem -M -t 20 $REF ${SAMPLE_ID}_R1.fastq.gz ${SAMPLE_ID}_R2.fastq.gz | samtools sort -@ 10 - -T ${SAMPLE_ID} -o $REPEAT_ALIGMENT_DIR/${SAMPLE_ID}.${SUFFIX}.sorted.bam
        samtools index $REPEAT_ALIGMENT_DIR/${SAMPLE_ID}.${SUFFIX}.sorted.bam
        samtools flagstat $REPEAT_ALIGMENT_DIR/${SAMPLE_ID}.${SUFFIX}.sorted.bam > $REPEAT_ALIGMENT_DIR/${SAMPLE_ID}.${SUFFIX}.sorted.bam.stats.out
        rm $REPEAT_ALIGMENT_DIR/${SAMPLE_ID}.${SUFFIX}.sorted.bam
        rm $REPEAT_ALIGMENT_DIR/${SAMPLE_ID}.${SUFFIX}.sorted.bam.bai
    fi
done < $FASTA_LIST

















