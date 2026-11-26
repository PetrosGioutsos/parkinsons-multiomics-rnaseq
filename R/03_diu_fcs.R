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

# ----- Step 3: Differential isoform usage (DEXSeq), dataset 1: PRKN ------
library(DEXSeq)

# isoform rows and info rows must be in the same order
all(rownames(prkn_f$counts) == prkn_f$info$transcript_id)   # must be TRUE

sample_prkn <- data.frame(row.names = colnames(prkn_f$counts),
                          condition = meta_prkn$condition)

dxd_prkn <- DEXSeqDataSet(countData  = round(as.matrix(prkn_f$counts)),
                          sampleData = sample_prkn,
                          design     = ~ sample + exon + condition:exon,
                          featureID  = prkn_f$info$transcript_id,   # isoform
                          groupID    = prkn_f$info$gene_id)         # gene

dxr_prkn <- DEXSeq(dxd_prkn,
                   reducedModel = ~ sample + exon,   # drops condition:exon = usage test
                   BPPARAM = BiocParallel::MulticoreParam(4))

saveRDS(dxr_prkn, "/home/projects/exam_2026_22140_22150/group_4/results/03_GSE90514_dexseq.rds")

# how many significant isoforms, and in how many genes?
sum(dxr_prkn$padj < 0.05, na.rm = TRUE)
length(unique(dxr_prkn$groupID[which(dxr_prkn$padj < 0.05)]))

# ----- Step 2: PCA (quality check), dataset 2: LRRK2 / sporadic ------
meta_lrrk2$condition <- factor(meta_lrrk2$condition, levels = c("Control", "PD"))
meta_lrrk2$batch     <- factor(meta_lrrk2$batch)
all(colnames(lrrk2_f$counts) == rownames(meta_lrrk2))   # must be TRUE

dds_lrrk2 <- DESeqDataSetFromMatrix(countData = round(as.matrix(lrrk2_f$counts)),
                                    colData   = meta_lrrk2,
                                    design    = ~ condition)
rld_lrrk2 <- rlog(dds_lrrk2)
plotPCA(rld_lrrk2, intgroup = "condition")
plotPCA(rld_lrrk2, intgroup = c("condition", "gender"))
plotPCA(rld_lrrk2, intgroup = c("condition", "batch"))

# ----- Step 3: Differential isoform usage (DEXSeq), dataset 2: LRRK2 / sporadic ------
all(rownames(lrrk2_f$counts) == lrrk2_f$info$transcript_id)   # must be TRUE

sample_lrrk2 <- data.frame(row.names = colnames(lrrk2_f$counts),
                           condition = meta_lrrk2$condition,
                           gender    = meta_lrrk2$gender,
                           batch     = meta_lrrk2$batch)

dxd_lrrk2 <- DEXSeqDataSet(countData  = round(as.matrix(lrrk2_f$counts)),
                           sampleData = sample_lrrk2,
                           design     = ~ sample + exon + batch:exon + gender:exon + condition:exon,
                           featureID  = lrrk2_f$info$transcript_id,   # isoform
                           groupID    = lrrk2_f$info$gene_id)         # gene

dxr_lrrk2 <- DEXSeq(dxd_lrrk2,
                    reducedModel = ~ sample + exon + batch:exon + gender:exon,   # usage test, corrected for batch + gender
                    BPPARAM = BiocParallel::MulticoreParam(4))

saveRDS(dxr_lrrk2, "/home/projects/exam_2026_22140_22150/group_4/results/03_GSE167529_dexseq.rds")

sum(dxr_lrrk2$padj < 0.05, na.rm = TRUE)
length(unique(dxr_lrrk2$groupID[which(dxr_lrrk2$padj < 0.05)]))
