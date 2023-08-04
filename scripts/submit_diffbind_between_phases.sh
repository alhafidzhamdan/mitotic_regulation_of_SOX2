#!/bin/bash

#$ -N diffbind_between_phases
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=350G
#$ -l h_rt=20:00:00

export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/Diffbind/bin:$PATH

VARIANT=$1

DIFFBIND=/exports/igmm/eddie/Glioblastoma-WGS/scripts/fame_gbm/R/diffbind_between_phases.R

Rscript $DIFFBIND $VARIANT



