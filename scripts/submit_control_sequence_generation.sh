#!/bin/bash

#$ -N control_seq
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=100G
#$ -l h_rt=2:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh

CONFIG=$1
REPEAT_FASTA_DIR=$2

source $CONFIG


## Generate 100 random sequences to generate NULL distribution
## For major satellite sequence:
### Generate random coordinates, excluding blacklisted regions:
cd ${REPEAT_FASTA_DIR}/major_sat_control

bedtools random -g $MM10_CHROM_SIZE_NO_ALT  -l 702 -n 1000 | \
    bedtools intersect -v -a - -b $BLACKLISTED_PEAKS_MM10 | \
    shuf -n 100 | awk '{print $1 "\t" $2 "\t" $3 "\t" $1 "_" $2 "_" $3}' > major_sat_control_100.bed

### Get nucleotide sequence from these coordinates
bedtools getfasta -fi $MM10_REF_FASTA -bed major_sat_control_100.bed -fo major_sat_control_100.fa

### Split them into individual fasta files:
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.major_sat_control.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < major_sat_control_100.fa

for file in *.fa; do samtools faidx $file;done
for file in *.fa; do bwa index $file;done

## For major satellite sequence:
### Generate random coordinates, excluding blacklisted regions:
cd ${REPEAT_FASTA_DIR}/minor_sat_control

bedtools random -g $MM10_CHROM_SIZE_NO_ALT  -l 360 -n 1000 | \
    bedtools intersect -v -a - -b $BLACKLISTED_PEAKS_MM10 | \
    shuf -n 100 | awk '{print $1 "\t" $2 "\t" $3 "\t" $1 "_" $2 "_" $3}' > minor_sat_control_100.bed

### Get nucleotide sequence from these coordinates
bedtools getfasta -fi $MM10_REF_FASTA -bed minor_sat_control_100.bed -fo minor_sat_control_100.fa

### Split them into individual fasta files:
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.minor_sat_control.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < minor_sat_control_100.fa

for file in *.fa; do samtools faidx $file;done
for file in *.fa; do bwa index $file;done




