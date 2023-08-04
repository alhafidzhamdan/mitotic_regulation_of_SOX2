#!/bin/bash

#$ -N diffbind_between_variants
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=350G
#$ -l h_rt=20:00:00

CONFIG=$1
PHASE=$2

source $CONFIG

DIFFBIND_SCRIPT=/exports/igmm/eddie/Glioblastoma-WGS/scripts/fame_gbm/R/diffbind_between_variants.R

Rscript $DIFFBIND_SCRIPT $PHASE



