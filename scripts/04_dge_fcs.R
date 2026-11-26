# ==============================================================================
# Task 04: Differential Gene Expression (DESeq2) & Functional Class Scoring (FCS)
# Datasets: GSE90514 & GSE167529 (25 mM glucose)
# ==============================================================================

library(DESeq2)
library(ggplot2)
library(fgsea)
library(msigdbr)

# ------------------------------------------------------------------------------
# 1. Load cleaned data from group_4/data
# ------------------------------------------------------------------------------
data_dir    <- "/home/projects/exam_2026_22140_22150/group_4/data"
results_dir <- "/home/projects/exam_2026_22140_22150/group_4/results"

gene_90514       <- readRDS(file.path(data_dir, "GSE90514_gene_counts.rds"))
meta_90514       <- readRDS(file.path(data_dir, "GSE90514_meta.rds"))

gene_167529_25mM <- readRDS(file.path(data_dir, "GSE167529_25mM_gene_counts.rds"))
meta_167529_25mM <- readRDS(file.path(data_dir, "GSE167529_25mM_meta.rds"))

# Ensure reference levels and scaled age are properly set
meta_90514$condition        <- relevel(factor(meta_90514$condition), ref = "Control")
meta_167529_25mM$condition  <- relevel(factor(meta_167529_25mM$condition), ref = "Control")
meta_167529_25mM$age_scaled <- as.numeric(scale(meta_167529_25mM$age))

stopifnot(all(rownames(meta_90514) == colnames(gene_90514)))
stopifnot(all(rownames(meta_167529_25mM) == colnames(gene_167529_25mM)))

# ------------------------------------------------------------------------------
# 2. Load MSigDB gene sets & build Ensembl -> Symbol lookup
# ------------------------------------------------------------------------------
msig_h     <- msigdbr(species = "Homo sapiens", category = "H")
pathways_h <- split(msig_h$gene_symbol, msig_h$gs_name)

msig_r     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
pathways_r <- split(msig_r$gene_symbol, msig_r$gs_name)

msig_go    <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

id_map_df  <- unique(rbind(
  msig_h[, c("ensembl_gene", "gene_symbol")],
  msig_r[, c("ensembl_gene", "gene_symbol")],
  msig_go[, c("ensembl_gene", "gene_symbol")]
))
id_map_df  <- id_map_df[!is.na(id_map_df$ensembl_gene) & id_map_df$ensembl_gene != "" &
                          !is.na(id_map_df$gene_symbol)  & id_map_df$gene_symbol != "", ]
id_map_df  <- id_map_df[!duplicated(id_map_df$ensembl_gene), ]
ens2sym    <- setNames(id_map_df$gene_symbol, id_map_df$ensembl_gene)

map_to_symbol <- function(ids) {
  clean_ids <- sub("\\..*$", "", ids)
  if (mean(grepl("^ENSG", clean_ids)) > 0.5) {
    mapped <- unname(ens2sym[clean_ids])
    return(ifelse(!is.na(mapped) & mapped != "", mapped, clean_ids))
  }
  return(ids)
}

# ------------------------------------------------------------------------------
# 3. Low-count filtering (keep genes with >= 10 counts in >= smallest group size)
# ------------------------------------------------------------------------------
keep_90514      <- rowSums(gene_90514 >= 10) >= 4
gene_90514_filt <- round(gene_90514[keep_90514, ])

keep_167529      <- rowSums(gene_167529_25mM >= 10) >= 9
gene_167529_filt <- round(gene_167529_25mM[keep_167529, ])

cat("\n=== Retained Genes After Low-Count Filtering ===\n")
cat("GSE90514:       ", nrow(gene_90514_filt), "of", nrow(gene_90514), "\n")
cat("GSE167529_25mM: ", nrow(gene_167529_filt), "of", nrow(gene_167529_25mM), "\n")

# ------------------------------------------------------------------------------
# 4. DESeq2 — GSE90514 (Unadjusted vs. Adjusted for Gender)
# ------------------------------------------------------------------------------
dds_90514_unadj <- DESeqDataSetFromMatrix(gene_90514_filt, meta_90514, ~ condition)
dds_90514_unadj <- DESeq(dds_90514_unadj)
res_90514_unadj <- as.data.frame(
  results(dds_90514_unadj, contrast = c("condition", "PD", "Control"), alpha = 0.05)
)
res_90514_unadj$gene_id <- rownames(res_90514_unadj)
res_90514_unadj$symbol  <- map_to_symbol(res_90514_unadj$gene_id)

dds_90514_adj <- DESeqDataSetFromMatrix(gene_90514_filt, meta_90514, ~ gender + condition)
dds_90514_adj <- DESeq(dds_90514_adj)
res_90514_adj <- as.data.frame(
  results(dds_90514_adj, contrast = c("condition", "PD", "Control"), alpha = 0.05)
)
res_90514_adj$gene_id <- rownames(res_90514_adj)
res_90514_adj$symbol  <- map_to_symbol(res_90514_adj$gene_id)
res_90514_adj <- res_90514_adj[order(res_90514_adj$padj), ]

# ------------------------------------------------------------------------------
# 5. DESeq2 — GSE167529 (25 mM) (Unadjusted vs. Adjusted for Batch + Gender + Age)
# ------------------------------------------------------------------------------
dds_167529_unadj <- DESeqDataSetFromMatrix(gene_167529_filt, meta_167529_25mM, ~ condition)
dds_167529_unadj <- DESeq(dds_167529_unadj)
res_167529_unadj <- as.data.frame(
  results(dds_167529_unadj, contrast = c("condition", "PD", "Control"), alpha = 0.05)
)
res_167529_unadj$gene_id <- rownames(res_167529_unadj)
res_167529_unadj$symbol  <- map_to_symbol(res_167529_unadj$gene_id)

dds_167529_adj <- DESeqDataSetFromMatrix(
  gene_167529_filt, meta_167529_25mM, ~ batch + gender + age_scaled + condition
)
dds_167529_adj <- DESeq(dds_167529_adj)
res_167529_adj <- as.data.frame(
  results(dds_167529_adj, contrast = c("condition", "PD", "Control"), alpha = 0.05)
)
res_167529_adj$gene_id <- rownames(res_167529_adj)
res_167529_adj$symbol  <- map_to_symbol(res_167529_adj$gene_id)
res_167529_adj <- res_167529_adj[order(res_167529_adj$padj), ]

# ------------------------------------------------------------------------------
# 6. Summarize DEG Counts & Cross-Study Overlap
# ------------------------------------------------------------------------------
summarize_deg <- function(df, label) {
  sig05  <- sum(!is.na(df$padj) & df$padj < 0.05)
  up05   <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange > 0)
  down05 <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange < 0)
  sig_fc <- sum(!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1)
  data.frame(
    Model = label,
    Total_DEGs_padj05 = sig05,
    Up = up05,
    Down = down05,
    DEGs_padj05_absLFC1 = sig_fc
  )
}

deg_summary <- rbind(
  summarize_deg(res_90514_unadj,  "GSE90514 Unadjusted (~ condition)"),
  summarize_deg(res_90514_adj,    "GSE90514 Adjusted (~ gender + condition)"),
  summarize_deg(res_167529_unadj, "GSE167529 Unadjusted (~ condition)"),
  summarize_deg(res_167529_adj,   "GSE167529 Adjusted (~ batch + gender + age + condition)")
)

cat("\n=== Summary of Differentially Expressed Genes (PD vs Control) ===\n")
print(deg_summary, row.names = FALSE)

cat("\n=== Top 10 DEGs — GSE90514 (Adjusted) ===\n")
print(head(res_90514_adj[, c("gene_id", "symbol", "baseMean", "log2FoldChange", "pvalue", "padj")], 10), row.names = FALSE)

cat("\n=== Top 10 DEGs — GSE167529 (Adjusted) ===\n")
print(head(res_167529_adj[, c("gene_id", "symbol", "baseMean", "log2FoldChange", "pvalue", "padj")], 10), row.names = FALSE)

# Cross-study DEG overlap (adjusted models, padj < 0.05)
deg_90514_ids  <- res_90514_adj$gene_id[!is.na(res_90514_adj$padj) & res_90514_adj$padj < 0.05]
deg_167529_ids <- res_167529_adj$gene_id[!is.na(res_167529_adj$padj) & res_167529_adj$padj < 0.05]
overlap_ids    <- intersect(deg_90514_ids, deg_167529_ids)

cat("\n=== Overlapping Significant DEGs (padj < 0.05) Across Both Studies ===\n")
cat("Number of overlapping DEGs:", length(overlap_ids), "\n")
if (length(overlap_ids) > 0) {
  ov_df <- merge(
    res_90514_adj[res_90514_adj$gene_id %in% overlap_ids, c("gene_id", "symbol", "log2FoldChange", "padj")],
    res_167529_adj[res_167529_adj$gene_id %in% overlap_ids, c("gene_id", "log2FoldChange", "padj")],
    by = "gene_id", suffixes = c("_GSE90514", "_GSE167529")
  )
  print(ov_df, row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 7. Functional Class Scoring (FCS / GSEA with fgsea + MSigDB Hallmark & Reactome)
# ------------------------------------------------------------------------------
run_fgsea_clean <- function(res_df, pathways) {
  df <- res_df[!is.na(res_df$symbol) & res_df$symbol != "" & !is.na(res_df$stat), ]
  df <- df[order(abs(df$stat), decreasing = TRUE), ]
  df <- df[!duplicated(df$symbol), ]
  ranks <- setNames(df$stat, df$symbol)
  ranks <- sort(ranks, decreasing = TRUE)
  
  set.seed(42)
  fg <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
  fg <- as.data.frame(fg[order(fg$padj), ])
  fg$leadingEdge <- sapply(fg$leadingEdge, paste, collapse = ";")
  fg
}

fcs_90514_h      <- run_fgsea_clean(res_90514_adj, pathways_h)
fcs_167529_h     <- run_fgsea_clean(res_167529_adj, pathways_h)

fcs_90514_react  <- run_fgsea_clean(res_90514_adj, pathways_r)
fcs_167529_react <- run_fgsea_clean(res_167529_adj, pathways_r)

cat("\n=== Top 10 Hallmark Pathways — GSE90514 (Adjusted) ===\n")
print(head(fcs_90514_h[, c("pathway", "NES", "pval", "padj", "size")], 10), row.names = FALSE)

cat("\n=== Top 10 Hallmark Pathways — GSE167529 (Adjusted) ===\n")
print(head(fcs_167529_h[, c("pathway", "NES", "pval", "padj", "size")], 10), row.names = FALSE)

hallmark_comp <- merge(
  fcs_90514_h[, c("pathway", "NES", "padj")],
  fcs_167529_h[, c("pathway", "NES", "padj")],
  by = "pathway", suffixes = c("_GSE90514", "_GSE167529")
)
hallmark_comp <- hallmark_comp[order(pmin(hallmark_comp$padj_GSE90514, hallmark_comp$padj_GSE167529)), ]

cat("\n=== Shared / Top Hallmark Pathways Across Both Studies ===\n")
print(head(hallmark_comp, 15), row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. Save Result Tables to group_4/results
# ------------------------------------------------------------------------------
write.csv(deg_summary,      file.path(results_dir, "04_deg_summary.csv"), row.names = FALSE)
write.csv(res_90514_adj,    file.path(results_dir, "04_GSE90514_dge_adjusted.csv"), row.names = FALSE)
write.csv(res_167529_adj,   file.path(results_dir, "04_GSE167529_dge_adjusted.csv"), row.names = FALSE)
write.csv(fcs_90514_h,      file.path(results_dir, "04_GSE90514_fcs_hallmark.csv"), row.names = FALSE)
write.csv(fcs_167529_h,     file.path(results_dir, "04_GSE167529_fcs_hallmark.csv"), row.names = FALSE)
write.csv(fcs_90514_react,  file.path(results_dir, "04_GSE90514_fcs_reactome.csv"), row.names = FALSE)
write.csv(fcs_167529_react, file.path(results_dir, "04_GSE167529_fcs_reactome.csv"), row.names = FALSE)
write.csv(hallmark_comp,    file.path(results_dir, "04_hallmark_cross_study.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. Generate Summary Plots (Volcano + Cross-Study Hallmark NES)
# ------------------------------------------------------------------------------
make_volcano <- function(df, title_txt) {
  df$sig <- ifelse(!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1, "padj < 0.05 & |LFC| > 1",
                   ifelse(!is.na(df$padj) & df$padj < 0.05, "padj < 0.05", "NS"))
  ggplot(df[!is.na(df$padj), ], aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_color_manual(values = c("padj < 0.05 & |LFC| > 1" = "#D55E00",
                                  "padj < 0.05" = "#0072B2", "NS" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    theme_bw() +
    ggtitle(title_txt) +
    labs(x = "log2 Fold Change (PD vs Control)", y = "-log10(p-value)", color = "Significance")
}

v1 <- make_volcano(res_90514_adj,  "GSE90514 Volcano — Adjusted (~ gender + condition)")
v2 <- make_volcano(res_167529_adj, "GSE167529 (25 mM) Volcano — Adjusted (~ batch + gender + age + condition)")

top_h_names <- head(hallmark_comp$pathway, 15)
h_plot_df <- rbind(
  data.frame(pathway = fcs_90514_h$pathway, NES = fcs_90514_h$NES, padj = fcs_90514_h$padj, Dataset = "GSE90514"),
  data.frame(pathway = fcs_167529_h$pathway, NES = fcs_167529_h$NES, padj = fcs_167529_h$padj, Dataset = "GSE167529 (25 mM)")
)
h_plot_df <- h_plot_df[h_plot_df$pathway %in% top_h_names, ]
h_plot_df$pathway <- gsub("^HALLMARK_", "", h_plot_df$pathway)
h_plot_df$sig <- ifelse(h_plot_df$padj < 0.05, "FDR < 0.05", "FDR >= 0.05")

p_hallmark <- ggplot(h_plot_df, aes(x = Dataset, y = reorder(pathway, NES), color = NES, size = -log10(padj))) +
  geom_point(aes(shape = sig)) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1)) +
  theme_bw() +
  ggtitle("Cross-Study Hallmark FCS Comparison (PD vs Control)") +
  labs(x = "", y = "", size = "-log10(FDR)", shape = "Significance")

while (!is.null(dev.list())) dev.off()
print(v1)
print(v2)
print(p_hallmark)

pdf_out <- file.path(results_dir, "04_dge_fcs_plots.pdf")
pdf(pdf_out, width = 8, height = 5.5)
print(v1)
print(v2)
print(p_hallmark)
dev.off()

# Copy PDF to personal working directory for easy browser download
file.copy(pdf_out, "04_dge_fcs_plots.pdf", overwrite = TRUE)