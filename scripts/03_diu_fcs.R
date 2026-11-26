# Task 3: Differential isoform usage (DIU) + FCS------
library(tidyverse)
library(DESeq2)

# --- Load data  ------
d <- "/home/projects/exam_2026_22140_22150/group_4/data/"
iso_prkn   <- readRDS(paste0(d, "GSE90514_iso_counts.rds"))
meta_prkn  <- readRDS(paste0(d, "GSE90514_meta.rds"))
iso_lrrk2  <- readRDS(paste0(d, "GSE167529_25mM_iso_counts.rds"))
meta_lrrk2 <- readRDS(paste0(d, "GSE167529_25mM_meta.rds"))

isoInfo <- read_rds("/home/projects/22140/exam_data/Homo_sapiens.GRCh38.107.chr_patch_hapl_scaff_tx_info.Rds")

# -- Step 1: filter isoforms -------
# keep isoforms with >= 10 reads in >= (smallest group) samples,
# and only genes with 2+ isoforms
filter_iso <- function(iso, min_samples, min_count = 10) {
  keep <- rowSums(iso >= min_count) >= min_samples
  iso  <- iso[keep, ]
  info <- isoInfo[match(rownames(iso), isoInfo$transcript_id),
                  c("transcript_id", "gene_id", "gene_name")]
  found <- !is.na(info$gene_id)
  iso  <- iso[found, ];  info <- info[found, ]
  multi <- info$gene_id %in% info$gene_id[duplicated(info$gene_id)]
  list(counts = iso[multi, ], info = info[multi, ])
}
prkn_f  <- filter_iso(iso_prkn,  min_samples = 4)   # 74,905 isoforms, 13,134 genes
lrrk2_f <- filter_iso(iso_lrrk2, min_samples = 9)   # 74,318 isoforms, 13,101 genes

# ----- Step 2: PCA (quality check), dataset 1: PRKN ------
meta_prkn$condition <- factor(meta_prkn$condition, levels = c("Control", "PD"))
all(colnames(prkn_f$counts) == rownames(meta_prkn))   # must be TRUE

dds_prkn <- DESeqDataSetFromMatrix(countData = round(as.matrix(prkn_f$counts)),
                                   colData   = meta_prkn,
                                   design    = ~ condition)
rld_prkn <- rlog(dds_prkn)
plotPCA(rld_prkn, intgroup = "condition")
plotPCA(rld_prkn, intgroup = c("condition", "gender"))
