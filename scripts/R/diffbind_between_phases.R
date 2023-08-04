#!/usr/bin/env Rscript

## Set wd:
setwd("/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/diffbind/extended_2q1vsq7")

## Load libraries:
.p <- c("tidyverse", "DiffBind")
suppressPackageStartupMessages(lapply(.p, require, character.only=T))

## Load args:
args <- commandArgs(trailingOnly = TRUE)
variant <- as.character(args[1])

## Load sample sheet:
sample_sheet <- read.delim("/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/diffbind/sheets/diffbind_sheet_replicates_normalised_2q1vsq7_extended.txt")
print(sample_sheet)

# Create a diffbind object
db_config = data.frame(th = 0.01, bUsePval=FALSE, RunParallel = FALSE)
db = dba(sampleSheet = sample_sheet, config = db_config)
# print(db)
# 
# ## Plot heatmap of peak occupancies 
# db$config$condition <- "Phase"
# db$config$tissue <- "Variant"
# pdf("Diffbind_replicate_no_bl6_heatmap.pdf")
# dba.plotHeatmap(db, colScheme="Blues")
# dev.off()
# 
# ## Plot replicate PCA, showing phases
# pdf("Diffbind_replicate_no_bl6_pca.pdf")
# dba.plotPCA(db, dotSize=1, attributes = DBA_FACTOR) 
# dev.off()

message("### FINDING DIFFERENTIALLY BOUND REGIONS ###\n\n")
## Perform differential binding analyse
## Inspired by https://github.com/j-andrews7/OneLinerOmics/blob/master/R/ChIPseq.R
## 1) Phases (Interphase vs Mitosis) (DBA_FACTOR in Diffbind term)
n.consensus <-2 ## This will count binding affinity that exist within peaks that occur in at least 2 replicates
rlevel <- DBA_FACTOR

sample_sheet_comp <- sample_sheet %>% filter(Tissue %in% variant)
samps <- dba(sampleSheet = sample_sheet_comp, minOverlap = n.consensus)

count <- dba.count(samps, minOverlap = n.consensus) 
saveRDS(count, paste0(variant, "_interphase_vs_mitosis_DBA_count.rds"))

cont <- dba.contrast(count, categories = rlevel)
saveRDS(cont, paste0(variant, "_interphase_vs_mitosis_DBA_cont.rds"))

results <- dba.analyze(cont, method = DBA_DESEQ2)
saveRDS(results, paste0(variant, "_interphase_vs_mitosis_DBA_results.rds"))

report <- dba.report(results, th = 1, bCalled = TRUE, bCounts = TRUE, method = DBA_DESEQ2, bCalledDetail = TRUE)
saveRDS(report, paste0(variant, "_interphase_vs_mitosis_DBA_report.rds"))
