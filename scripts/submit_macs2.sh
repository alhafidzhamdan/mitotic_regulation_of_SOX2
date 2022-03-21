!#/bin/bash

#$ -N macs2
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=2:00:00
#$ -pe sharedmem 8

CONFIG=$1
IDS=$2
BAM_DIR=$3
OUTPUT_DIR=$4

source $CONFIG

unset MODULEPATH
. /etc/profile.d/modules.sh

module load igmm/apps/MACS2/2.1.1  

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1`
CHIP_BAM=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 2`
INPUT_BAM=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 3`

cd $OUTPUT_DIR

macs2 callpeak \
    -t $BAM_DIR/$CHIP_BAM \
    -c $BAM_DIR/$INPUT_BAM \
    -n ${SAMPLE_ID} \
    -p 1e-3 \
    -g mm
    
    
    
    
    
