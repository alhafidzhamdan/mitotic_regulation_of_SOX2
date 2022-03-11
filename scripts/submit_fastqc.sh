#!/bin/bash

# To run this script, do 
# qsub -t 1-n -tc N submit_fastqc.sh <DIR> <IDS>
# IDS point to lane fastq files
#
#$ -N fastqc
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -pe sharedmem 8G
#$ -l h_rt=4:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh

module load igmm/apps/FastQC/0.11.9

CONFIG=$1
IDS=$2
DIR=$3

FASTQ=`head -n $SGE_TASK_ID $IDS | tail -n 1`

cd $DIR

fastqc -t 20 $FASTQ

