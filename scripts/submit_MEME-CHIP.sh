#!/bin/bash

# To run this script, do 
# qsub -t n submit_MEME-CHIP.sh
#
#$ -N MEME-CHIP
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=12G
#$ -pe sharedmem 12
#$ -l h_rt=48:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh

BED_LIST=$1
MEME_CHIP_DIR=$2
FASTA_DIR=$3

####
SAMPLE_ID=`head -n $SGE_TASK_ID $BED_LIST | cut -f 1 | tail -n 1`
BED_FILE=`head -n $SGE_TASK_ID $BED_LIST | cut -f 2 | tail -n 1`
MEME_CHIP_DIR_SAMPLE=$MEME_CHIP_DIR/${SAMPLE_ID}
FASTA_FILE=$FASTA_DIR/${SAMPLE_ID}.fasta

## For TomTom and Centrimo (not for de novo discovery)
DB_DIR=/exports/igmm/eddie/Glioblastoma-WGS/resources/JASPAR2022_CORE_non-redundant.meme

BEDTOOLS=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/snakemake/bin/bedtools
MEME_CHIP=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/snakemake/bin/meme-chip
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/scripts/bin:/exports/igmm/eddie/Glioblastoma-WGS/scripts/libexec/meme-5.3.0:$PATH
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/meme-chip/bin:$PATH

#### MEME_CHIP 
#### Based on https://www.biostars.org/p/401084/

### $BEDTOOLS getfasta -fi $MM10_REF_FASTA -bed $BED_FILE > $FASTA_FILE

$MEME_CHIP -oc $MEME_CHIP_DIR_SAMPLE -db $DB_DIR $FASTA_FILE -meme-p 12 -meme-searchsize 100000 -streme-pvt 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 -spamo-skip
    



    
    
    
    
    
    
    
