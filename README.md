This repository contains the scripts used to analyse the ChIP-seq data in the paper **Mitotic Regulation of SOX2** (preprint ...)

All the bash/R scripts can be found in the folder **scripts**

#### Step 1: QC fastq files
`qsub -t 1-N -tc N submit_fastqc.sh <DIR> <IDS>`
#### Step 2: Trim adaptors +/- demultiplexing + repeat QC
`qsub -t 1-N submit_fastp.sh <IDS>`
#### Step 3: Align to a reference genome
`qsub -t 1-N -tc N submit_bwa.sh <CONFIG> <IDS>`
#### Step 4: QC alignment
`qsub submit_plotCorrelation.sh <CONFIG> <BATCH>`
#### Step 5: Generate coverage tracks
`qsub -t 1-N -tc N submit_bamCoverage.sh <CONFIG> <IDS>`
#### Step 6: Filter duplicated and blacklisted regions
`qsub -t 1-N -tc N submit_filter_bams.sh <CONFIG> <IDS>`
#### Step 7: Call peaks +/- merge replicates + QC +/- IDR
`qsub -t 1-N -tc N submit_macs2IDR.sh <CONFIG> <IDS> <CHIP> <BATCH>`
#### Step 8: Downstream analyses



