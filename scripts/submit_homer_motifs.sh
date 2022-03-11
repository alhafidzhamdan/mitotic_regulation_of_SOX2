#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_homer_motif.sh
#
#$ -N homer_motif
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -l h_rt=10:00:00

## BED_LIST contains full paths to bed files 

CONFIG=$1
BED=$2
OUTDIR=$3

source $CONFIG

#SAMPLE=`head -n $SGE_TASK_ID $BED_LIST | cut -f 1| tail -n 1`
#BED=`head -n $SGE_TASK_ID $BED_LIST | cut -f 2| tail -n 1`
##BED_FILE=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/SOX_TF/downstream/beds/$BED
##OUTDIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/SOX_TF/downstream/homer/motifs/${BED%.bed}.motif_output

## Get enriched motifs
findMotifsGenome.pl \
     $BED \
     hg38 \
     $OUTDIR \
     -size given
     
###scanMotifGenomeWide.pl $MOTIF $INPUT_DIR -keepAll -homer2 -bed > ${SAMPLE_ID}_homer_called_motif.bed
### Use scrambleFasta.pl to generate background sequence

