# ==============================================================================
# Task 02: Quality Control & PCA (aligned with Week 5 TA solutions)
# Datasets: GSE90514 & GSE167529 (25 mM glucose)
# ==============================================================================

library(DESeq2)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Load cleaned data from group_4/data
# ------------------------------------------------------------------------------
data_dir    <- "/home/projects/exam_2026_22140_22150/group_4/data"
results_dir <- "/home/projects/exam_2026_22140_22150/group_4/results"

gene_90514       <- readRDS(file.path(data_dir, "GSE90514_gene_counts.rds"))
meta_90514       <- readRDS(file.path(data_dir, "GSE90514_meta.rds"))

gene_167529_25mM <- readRDS(file.path(data_dir, "GSE167529_25mM_gene_counts.rds"))
meta_167529_25mM <- readRDS(file.path(data_dir, "GSE167529_25mM_meta.rds"))
meta_167529      <- readRDS(file.path(data_dir, "GSE167529_all42_meta.rds"))

# Ensure age, age_scaled, and age_group are present in GSE167529 metadata
if (!"age" %in% colnames(meta_167529_25mM) || any(is.na(meta_167529_25mM$age))) {
  aad_lookup <- c(
    "SP-01" = 63, "SP-02" = 55, "SP-04" = 46, "SP-06" = 44, "SP-09" = 66,
    "SP-11" = 48, "SP-12" = 63, "SP-13" = 68, "FOX-073-001" = 66, "SP-16" = 51,
    "SP-17" = 52, "FOX-019-001" = 71, "FOX-075-001" = 51, "FOX-015-001" = 68,
    "SP-08" = 66, "IPD-008" = 65, "CTRL-001" = 68, "CTRL-002" = 73,
    "CTRL-003" = 55, "FOX-076-001" = 61, "FOX-069-001" = 64
  )
  extract_donor <- function(titles) {
    inside <- sub(".*\\(([^)]+)\\).*", "\\1", titles)
    trimws(gsub("\\s+(high|low)$", "", inside, ignore.case = TRUE))
  }
  
  meta_167529$donor_id        <- extract_donor(meta_167529$title)
  meta_167529$age             <- as.numeric(aad_lookup[meta_167529$donor_id])
  meta_167529_25mM$donor_id   <- extract_donor(meta_167529_25mM$title)
  meta_167529_25mM$age        <- as.numeric(aad_lookup[meta_167529_25mM$donor_id])
}

meta_167529$age_scaled      <- as.numeric(scale(meta_167529$age))
meta_167529$age_group       <- factor(ifelse(meta_167529$age >= 60, ">=60", "<60"))
meta_167529_25mM$age_scaled <- as.numeric(scale(meta_167529_25mM$age))
meta_167529_25mM$age_group  <- factor(ifelse(meta_167529_25mM$age >= 60, ">=60", "<60"))

saveRDS(meta_167529_25mM, file.path(data_dir, "GSE167529_25mM_meta.rds"))
saveRDS(meta_167529,      file.path(data_dir, "GSE167529_all42_meta.rds"))

has_age_90514 <- "age" %in% colnames(meta_90514) && !any(is.na(meta_90514$age))
if (has_age_90514) {
  meta_90514$age        <- as.numeric(meta_90514$age)
  meta_90514$age_scaled <- as.numeric(scale(meta_90514$age))
  meta_90514$age_group  <- factor(ifelse(meta_90514$age >= 60, ">=60", "<60"))
  saveRDS(meta_90514, file.path(data_dir, "GSE90514_meta.rds"))
}

stopifnot(all(rownames(meta_90514) == colnames(gene_90514)))
stopifnot(all(rownames(meta_167529_25mM) == colnames(gene_167529_25mM)))
stopifnot(!any(is.na(meta_167529_25mM$age)))

# ------------------------------------------------------------------------------
# 2. QC: Design Balance & Confounder Checks
# ------------------------------------------------------------------------------
cat("\n=== GSE90514 Design Table (Condition x Gender) ===\n")
print(table(meta_90514$condition, meta_90514$gender))
if (has_age_90514) {
  print(table(meta_90514$condition, meta_90514$age_group))
  print(tapply(meta_90514$age, meta_90514$condition, summary))
}

cat("\n=== GSE167529 Full Design Tables (Glucose x Condition / Batch / Run) ===\n")
print(table(meta_167529$glucose, meta_167529$condition))
print(table(meta_167529$glucose, meta_167529$batch))
print(table(meta_167529$glucose, meta_167529$run))

cat("\n=== GSE167529 (25 mM only) Confounder Balance ===\n")
print(table(meta_167529_25mM$condition, meta_167529_25mM$gender))
print(table(meta_167529_25mM$condition, meta_167529_25mM$batch))
print(table(meta_167529_25mM$condition, meta_167529_25mM$run))
print(table(meta_167529_25mM$condition, meta_167529_25mM$age_group))

cat("\n=== GSE167529 (25 mM only) Age Distribution by Condition ===\n")
print(tapply(meta_167529_25mM$age, meta_167529_25mM$condition, summary))

# ------------------------------------------------------------------------------
# 3. QC: Library Sizes (Total Read Counts per Sample, in millions)
# ------------------------------------------------------------------------------
cat("\n=== Library Sizes (in millions of reads) ===\n")
cat("GSE90514:\n")
print(round(colSums(gene_90514) / 1e6, 2))

cat("GSE167529 (25 mM) summary:\n")
print(summary(round(colSums(gene_167529_25mM) / 1e6, 2)))

# ------------------------------------------------------------------------------
# 4. Build PCA Plots — GSE90514 (Week 5 TA style: rlog + plotPCA)
# ------------------------------------------------------------------------------
dds_90514_base <- DESeqDataSetFromMatrix(
  countData = round(gene_90514),
  colData   = meta_90514,
  design    = ~ condition
)
rld_90514 <- rlog(dds_90514_base)

p1_90514 <- plotPCA(rld_90514, intgroup = c("condition")) +
  ggtitle("GSE90514 PCA — Condition only")

p2_90514 <- plotPCA(rld_90514, intgroup = c("gender", "condition")) +
  ggtitle("GSE90514 PCA — Gender + Condition")

plots_90514 <- list(p1_90514, p2_90514)

if (has_age_90514) {
  p3_90514 <- plotPCA(rld_90514, intgroup = "age") +
    ggtitle("GSE90514 PCA — Age (continuous)")
  p4_90514 <- plotPCA(rld_90514, intgroup = c("gender", "age_group", "condition")) +
    ggtitle("GSE90514 PCA — Gender + Age + Condition")
  plots_90514 <- c(plots_90514, list(p3_90514, p4_90514))
}

# ------------------------------------------------------------------------------
# 5. Build PCA Plots — GSE167529 (25 mM) (Week 5 TA style: rlog + plotPCA)
# ------------------------------------------------------------------------------
dds_167529_base <- DESeqDataSetFromMatrix(
  countData = round(gene_167529_25mM),
  colData   = meta_167529_25mM,
  design    = ~ condition
)
rld_167529 <- rlog(dds_167529_base)

p1_167529 <- plotPCA(rld_167529, intgroup = c("condition")) +
  ggtitle("GSE167529 (25 mM) PCA — Condition only")

p2_167529 <- plotPCA(rld_167529, intgroup = c("gender", "condition")) +
  ggtitle("GSE167529 (25 mM) PCA — Gender + Condition")

p3_167529 <- plotPCA(rld_167529, intgroup = c("condition", "batch")) +
  ggtitle("GSE167529 (25 mM) PCA — Condition + Batch")

p4_167529 <- plotPCA(rld_167529, intgroup = c("condition", "run")) +
  ggtitle("GSE167529 (25 mM) PCA — Condition + Sequencing Run")

p5_167529 <- plotPCA(rld_167529, intgroup = "age") +
  ggtitle("GSE167529 (25 mM) PCA — Age (continuous)")

p6_167529 <- plotPCA(rld_167529, intgroup = c("condition", "age_group")) +
  ggtitle("GSE167529 (25 mM) PCA — Condition + Age Group")

p7_167529 <- plotPCA(rld_167529, intgroup = c("gender", "age_group", "condition")) +
  ggtitle("GSE167529 (25 mM) PCA — Gender + Age + Condition")

plots_167529 <- list(p1_167529, p2_167529, p3_167529, p4_167529, p5_167529, p6_167529, p7_167529)

# ------------------------------------------------------------------------------
# 6. Display all plots in RStudio & save PDF to group_4/results
# ------------------------------------------------------------------------------
while (!is.null(dev.list())) dev.off()

plots <- c(plots_90514, plots_167529)
for (p in plots) print(p)

pdf_out <- file.path(results_dir, "02_qc_pca_plots.pdf")
pdf(pdf_out, width = 8, height = 5.5)
for (p in plots) print(p)
dev.off()

# Copy PDF to personal working directory for easy browser download
file.copy(pdf_out, "02_qc_pca_plots.pdf", overwrite = TRUE)