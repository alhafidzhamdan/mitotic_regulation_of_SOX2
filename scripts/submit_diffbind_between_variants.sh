#!/bin/bash

#$ -N diffbind_between_variants
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=350G
#$ -l h_rt=20:00:00

export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/Diffbind/bin:$PATH

PHASE=$1

DIFFBIND=/exports/igmm/eddie/Glioblastoma-WGS/scripts/fame_gbm/R/diffbind_between_variants.R

Rscript $DIFFBIND $PHASE



