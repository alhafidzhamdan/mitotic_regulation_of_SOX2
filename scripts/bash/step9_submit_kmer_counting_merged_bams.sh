#!/bin/bash

## Count kmers at merged bam level
#$ -N merged_bams_kmercounting
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=150G
####$ -pe sharedmem 4
#$ -l h_rt=24:00:00

CONFIG=$1
IDS=$2

unset MODULEPATH
. /etc/profile.d/modules.sh

source $CONFIG
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py2/bin:$PATH

FACTOR=`head -n $SGE_TASK_ID $IDS | tail -n 1`
ALL_READS_BAM=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/merged/sorted/${FACTOR}.merged.sorted.bam
ALL_READS_FASTQ=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/merged/sorted/all/${FACTOR}.merged.sorted.fq
ALL_READS_FASTA=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/merged/sorted/all/${FACTOR}.merged.sorted.fa
JELLYFISH_RESULTS_ALL=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/repeats/jellyfish/all

## Jellyfish installed using conda https://anaconda.org/bioconda/jellyfish
## Guides: https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md#Counting-high-frequency-k-mers

echo "PROCESSING ALL READS..."
if [ ! -f $ALL_READS_FASTA ]; then
    ## Convert to fastq file
    bamtofastq filename=$ALL_READS_BAM > $ALL_READS_FASTQ
    awk '{ if ('NR%4==1' || 'NR%4==2'){ gsub("@",">",$1); print }}' $ALL_READS_FASTQ > $ALL_READS_FASTA
fi

## Count k-mers and dump in csv file:
## 12-mer represents tandem repeats in telomeric ends
## 17-mer represents the 17 nucleotides in the centromeric protein CENBP

cd $JELLYFISH_RESULTS_ALL
for kmer_length in 12 17; do
    if [ ! -f ${FACTOR}_all_${kmer_length}-mer_counts.csv ]; then
         jellyfish count -m ${kmer_length} -s 100M -t 10 -C --min-quality 10 $ALL_READS_FASTA -o ${FACTOR}_all_${kmer_length}-mer_counts.jf 
         jellyfish dump -c ${FACTOR}_all_${kmer_length}-mer_counts.jf -o ${FACTOR}_all_${kmer_length}-mer_counts.csv
    fi
done











