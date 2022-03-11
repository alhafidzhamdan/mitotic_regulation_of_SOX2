This repository contains the scripts used to analyse the ChIP-seq data from the paper **Mitotic Regulation of SOX2** (preprint ...)

#### Step 1: QC fastq files
`qsub -t 1-N -tc N submit_fastqc.sh <DIR> <IDS>`
#### Step 2: Trim adaptors +/- demultiplexing + repeat QC
`qsub -t 1-N submit_fastp.sh <IDS>`
#### Step 3: Align to a reference genome
`qsub -t 1-N -tc N submit_bwa.sh <CONFIG> <IDS>`
#### Step 4: QC alignment
`qsub submit_plotCorrelation.sh <CONFIG> <BATCH>`
#### Step 5: Filter duplicated and blacklisted regions
`qsub -t 1-N -tc N submit_filter_bams.sh <CONFIG> <IDS>`
#### Step 6: Generate coverage tracks
#### Step 7: Call peaks +/- merge replicates + QC +/- IDR
#### Step 8: Downstream analyses



