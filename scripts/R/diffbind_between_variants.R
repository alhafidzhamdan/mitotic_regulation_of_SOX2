#!/usr/bin/env Rscript

## Set wd:
setwd("/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/diffbind/extended_2q1vsq7")

## Load libraries:
.p <- c("tidyverse", "DiffBind")
suppressPackageStartupMessages(lapply(.p, require, character.only=T))

## Load args:
args <- commandArgs(trailingOnly = TRUE)
phase <- as.character(args[1])

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
n.consensus <-2 ## This will count binding affinity that exist within peaks that occur in at least 2 samples
## Perform differential binding analyse
## Inspired by https://github.com/j-andrews7/OneLinerOmics/blob/master/R/ChIPseq.R

rlevel <- DBA_TISSUE

## For WT vs 3A 
message(paste0("Running comparisons for WT vs 3A at ", phase))
sample_sheet_comp <- sample_sheet %>% filter(Condition %in% phase) %>% filter(Tissue %in% c("WT", "3A"))
samps <- dba(sampleSheet = sample_sheet_comp, minOverlap = n.consensus)

count <- dba.count(samps, minOverlap = n.consensus) 
saveRDS(count, paste0("WT_vs_3A_", phase, "_DBA_count.rds"))

cont <- dba.contrast(count, categories = rlevel)
saveRDS(cont, paste0("WT_vs_3A_", phase, "_DBA_cont.rds"))

results <- dba.analyze(cont, method = DBA_DESEQ2)
saveRDS(results, paste0("WT_vs_3A_", phase, "_DBA_results.rds"))

report <- dba.report(results, th = 1, bCalled = TRUE, bCounts = TRUE, method = DBA_DESEQ2, bCalledDetail = TRUE)
saveRDS(report, paste0("WT_vs_3A_", phase, "_DBA_report.rds"))

## For WT vs 5MK
message(paste0("Running comparisons for WT vs 5MK at ", phase))
sample_sheet_comp <- sample_sheet %>% filter(Condition %in% phase) %>% filter(Tissue %in% c("WT", "5MK"))
samps <- dba(sampleSheet = sample_sheet_comp, minOverlap = n.consensus)

count <- dba.count(samps, minOverlap = n.consensus)
saveRDS(count, paste0("WT_vs_5MK_", phase, "_DBA_count.rds"))

cont <- dba.contrast(count, categories = rlevel)
saveRDS(cont, paste0("WT_vs_5MK_", phase, "_DBA_cont.rds"))

results <- dba.analyze(cont, method = DBA_DESEQ2)
saveRDS(results, paste0("WT_vs_5MK_", phase, "_DBA_results.rds"))

report <- dba.report(results, th = 1, bCalled = TRUE, bCounts = TRUE, method = DBA_DESEQ2, bCalledDetail = TRUE)
saveRDS(report, paste0("WT_vs_5MK_", phase, "_DBA_report.rds"))

## For 3A vs 5MK
message(paste0("Running comparisons for 3A vs 5MK at ", phase))
sample_sheet_comp <- sample_sheet %>% filter(Condition %in% phase) %>% filter(Tissue %in% c("3A", "5MK"))
samps <- dba(sampleSheet = sample_sheet_comp, minOverlap = n.consensus)

count <- dba.count(samps, minOverlap = n.consensus)
  saveRDS(count, paste0("3A_vs_5MK_", phase, "_DBA_count.rds"))

cont <- dba.contrast(count, categories = rlevel)
saveRDS(cont, paste0("3A_vs_5MK_", phase, "_DBA_cont.rds"))

results <- dba.analyze(cont, method = DBA_DESEQ2)
saveRDS(results, paste0("3A_vs_5MK_", phase, "_DBA_results.rds"))

report <- dba.report(results, th = 1, bCalled = TRUE, bCounts = TRUE, method = DBA_DESEQ2, bCalledDetail = TRUE)
saveRDS(report, paste0("3A_vs_5MK_", phase, "_DBA_report.rds"))

