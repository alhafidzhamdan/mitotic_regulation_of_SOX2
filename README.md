This repository contains the scripts used to analyse the ChIP-seq data from the paper **Mitotic Regulation of SOX2** (preprint ...)

#### Step 1: QC fastq files
`qsub -t 1-n -tc N submit_fastqc.sh <DIR> <IDS>`
#### Step 2: Trim adaptors +/- demultiplexing + repeat QC
#### Step 3: Align to a reference genome
#### Step 4: QC alignment
#### Step 5: Filter duplicated and blacklisted regions
#### Step 6: Generate coverage tracks
#### Step 7: Call peaks +/- merge replicates + QC +/- IDR
#### Step 8: Downstream analyses



