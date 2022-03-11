#!/bin/bash

# To run this script, do 
# qsub submit_Diffbind.sh
#
#$ -N diffbind
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=48:00:00

export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/Diffbind/bin:$PATH

DIFFBIND=/exports/igmm/eddie/Glioblastoma-WGS/scripts/fame_gbm/R/Diffbind.R

Rscript $DIFFBIND 

