#!/bin/bash

#$ -N diffbind_between_phases
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=350G
#$ -l h_rt=20:00:00

CONFIG=$1
VARIANT=$2

source $CONFIG

DIFFBIND_SCRIPT=/exports/igmm/eddie/Glioblastoma-WGS/scripts/fame_gbm/R/diffbind_between_phases.R

Rscript $DIFFBIND_SCRIPT $VARIANT



