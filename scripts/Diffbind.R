#!/usr/bin/env Rscript

## Set wd:
setwd("/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/params")

## Load libraries:
library(DiffBind, quietly = TRUE)

## Load args:
## args <- commandArgs(trailingOnly = TRUE)

## Load sample sheet:
samples <- read.delim("diffbind_metadata_v4.tsv")
print(samples)

# Create a diffbind object
db_config = data.frame(th = 0.01, RunParallel = FALSE)
db = dba(sampleSheet = samples, config = db_config)

# Create consensus peak sets for each of GSC and normal tissues
# Consensus peakset for different cell types - present in at least half of samples of the same cell type
# Based on https://github.com/MarioniLab/Spermatogenesis2018/blob/master/Analysis/CutAndRun/DiffBind_CUTnRUN_H3K4me3.Rmd
# and https://www.biostars.org/p/464713/

message(" =========== get consensus peak sets ============")
Consensus <- dba.peakset(db, consensus = DBA_CONDITION, minOverlap = 0.5)
Consensus <- dba(Consensus, mask=Consensus$masks$Consensus, minOverlap=1)
Consensus_peaks <- dba.peakset(Consensus, bRetrieve = TRUE)

# Calculate affinity binding matrix at each consensus site
message(" ========== calculate binding matrix ========== ")

# calculate a binding matrix with scores based on read counts for every sample (affinity scores),
# rather than confidence scores for only those peaks called in a specific sample (occupancy scores)
Consensus_counts = dba.count(db, peaks = Consensus_peaks, score = DBA_SCORE_TMM_READS_FULL, bUseSummarizeOverlaps=TRUE, bRemoveDuplicates = TRUE)

# save
outDir = "/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/"
dbCountRData = paste0(outDir, "db.count.RData")
save(Consensus_counts, file = dbCountRData)

# end

