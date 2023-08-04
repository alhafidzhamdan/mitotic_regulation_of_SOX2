#!/bin/bash

## Count kmers at replicate bam level
#$ -N replicate_kmercounting
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=150G
####$ -pe sharedmem 4
#$ -l h_rt=48:00:00

CONFIG=$1
IDS=$2

unset MODULEPATH
. /etc/profile.d/modules.sh

source $CONFIG
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py2/bin:$PATH

FACTOR=`head -n $SGE_TASK_ID $IDS | tail -n 1` ## This points to sample ID of replicate bams

ALL_READS_BAM=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/${FACTOR}.sorted.bam
ALL_READS_FASTQ=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/all/${FACTOR}.sorted.fq
ALL_READS_FASTA=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/all/${FACTOR}.sorted.fa
JELLYFISH_RESULTS_ALL=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/repeats/jellyfish/all/replicates

UNMAPPED_READS_BAM=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/unmapped/${FACTOR}.sorted.unmapped.bam
UNMAPPED_READS_FASTQ=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/unmapped/${FACTOR}.sorted.unmapped.fq
UNMAPPED_READS_FASTA=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/unmapped/${FACTOR}.sorted.unmapped.fa
JELLYFISH_RESULTS_UNMAPPED=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/repeats/jellyfish/unmapped/replicates

MAPPED_READS_BAM=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/mapped/${FACTOR}.sorted.mapped.bam
MAPPED_READS_FASTQ=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/mapped/${FACTOR}.sorted.mapped.fq
MAPPED_READS_FASTA=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/bams/replicates/sorted/mapped/${FACTOR}.sorted.mapped.fa
JELLYFISH_RESULTS_MAPPED=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/repeats/jellyfish/mapped/replicates

## Jellyfish installed using conda https://anaconda.org/bioconda/jellyfish
## Guides: https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md#Counting-high-frequency-k-mers

########## ALL READS #########
echo "PROCESSING ALL READS..."
if [ ! -f $ALL_READS_FASTA ]; then
    ## Convert to fastq file
    bamtofastq filename=$ALL_READS_BAM > $ALL_READS_FASTQ
    awk '{ if ('NR%4==1' || 'NR%4==2'){ gsub("@",">",$1); print }}' $ALL_READS_FASTQ > $ALL_READS_FASTA
fi

## Count N k-mers within all (both mapped and unmapped) genomic regions, and dump in csv file:
cd $JELLYFISH_RESULTS_ALL
for kmer_length in 12 17; do
    if [ ! -f ${FACTOR}_all_${kmer_length}-mer_counts.csv ]; then
         jellyfish count -m ${kmer_length} -s 100M -t 10 -C --min-quality 10 $ALL_READS_FASTA -o ${FACTOR}_all_${kmer_length}-mer_counts.jf 
         jellyfish dump -c ${FACTOR}_all_${kmer_length}-mer_counts.jf -o ${FACTOR}_all_${kmer_length}-mer_counts.csv
    fi
done

########## UNMAPPED READS #########
echo "PROCESSING UNMAPPED READS..."
if [ ! -f $UNMAPPED_READS_FASTA ]; then
    ## Extract unmapped reads from bam file
    samtools view -b -f 4 $ALL_READS_BAM > $UNMAPPED_READS_BAM

    ## Convert to fastq file
    bamtofastq filename=$UNMAPPED_READS_BAM > $UNMAPPED_READS_FASTQ
    awk '{ if ('NR%4==1' || 'NR%4==2'){ gsub("@",">",$1); print }}' $UNMAPPED_READS_FASTQ > $UNMAPPED_READS_FASTA
fi

## Count k-mers within the unmapped regions, and dump in csv file:
cd $JELLYFISH_RESULTS_UNMAPPED
for kmer_length in 17; do
    if [ ! -f ${FACTOR}_unmapped_${kmer_length}-mer_counts.csv ]; then
         jellyfish count -m ${kmer_length} -s 100M -t 10 -C --min-quality 10 $UNMAPPED_READS_FASTA -o ${FACTOR}_unmapped_${kmer_length}-mer_counts.jf 
         jellyfish dump -c ${FACTOR}_unmapped_${kmer_length}-mer_counts.jf -o ${FACTOR}_unmapped_${kmer_length}-mer_counts.csv
    fi
done

########## MAPPED READS #########
echo "PROCESSING MAPPED READS..."
if [ ! -f $MAPPED_READS_FASTA ]; then
    ## Extract mapped reads from bam file
    samtools view -h -F 4 -b $ALL_READS_BAM > $MAPPED_READS_BAM

    ## Convert to fastq file
    bamtofastq filename=$MAPPED_READS_BAM > $MAPPED_READS_FASTQ
    awk '{ if ('NR%4==1' || 'NR%4==2'){ gsub("@",">",$1); print }}' $MAPPED_READS_FASTQ > $MAPPED_READS_FASTA
fi

## Count k-mers within the mapped regions, and dump in csv file:
cd $JELLYFISH_RESULTS_MAPPED
for kmer_length in 17; do
    if [ ! -f ${FACTOR}_mapped_${kmer_length}-mer_counts.csv ]; then
        jellyfish count -m ${kmer_length} -s 100M -t 10 -C --min-quality 10 $MAPPED_READS_FASTA -o ${FACTOR}_mapped_${kmer_length}-mer_counts.jf 
        jellyfish dump -c ${FACTOR}_mapped_${kmer_length}-mer_counts.jf -o ${FACTOR}_mapped_${kmer_length}-mer_counts.csv
    fi
done










