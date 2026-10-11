# 02_qc.R
# Library-level QC, low-count filtering, variance-stabilised PCA and
# sample-to-sample distances across all 16 samples.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})
source(file.path("R", "utils.R"))

counts  <- readRDS(file.path("data", "processed", "counts.rds"))
samples <- readRDS(file.path("data", "processed", "samples.rds"))

# ---- Library size and detected genes ----------------------------------------
lib <- data.frame(
  sample     = colnames(counts),
  title      = samples$title,
  cell_line  = samples$cell_line,
  treatment  = samples$treatment,
  total_reads_M  = colSums(counts) / 1e6,
  genes_detected = colSums(counts > 0),
  genes_ge10     = colSums(counts >= 10)
)
write.csv(lib, file.path("results", "qc_library_stats.csv"), row.names = FALSE)

p_lib <- ggplot(lib, aes(x = reorder(title, total_reads_M), y = total_reads_M, fill = treatment)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = treatment_cols) +
  labs(x = NULL, y = "Assigned reads (millions)", fill = "Treatment",
       title = "Library size per sample") +
  theme_portfolio()
save_fig(p_lib, "01_library_sizes.png", width = 7, height = 5)

# ---- Filtering ---------------------------------------------------------------
# Keep genes with >= 10 counts in at least 4 samples (the size of the smallest group).
keep <- rowSums(counts >= 10) >= 4
message(sprintf("Genes retained after filtering: %d of %d", sum(keep), length(keep)))
counts_f <- counts[keep, ]

# ---- Count distributions before/after normalisation ---------------------------
dds_all <- DESeqDataSetFromMatrix(counts_f, samples, design = ~ cell_line + treatment)
dds_all <- estimateSizeFactors(dds_all)
sf <- data.frame(title = samples$title, size_factor = sizeFactors(dds_all))
write.csv(sf, file.path("results", "qc_size_factors.csv"), row.names = FALSE)

long_counts <- function(m, label) {
  data.frame(title = rep(samples$title, each = nrow(m)),
             treatment = rep(samples$treatment, each = nrow(m)),
             value = as.vector(log2(m + 1)), stage = label)
}
dist_df <- rbind(long_counts(counts(dds_all, normalized = FALSE), "Raw counts"),
                 long_counts(counts(dds_all, normalized = TRUE),  "Size-factor normalised"))
dist_df$stage <- factor(dist_df$stage, levels = c("Raw counts", "Size-factor normalised"))
p_box <- ggplot(dist_df, aes(x = title, y = value, fill = treatment)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  facet_wrap(~ stage, ncol = 1) +
  scale_fill_manual(values = treatment_cols) +
  labs(x = NULL, y = "log2(count + 1)", fill = "Treatment",
       title = "Per-sample count distributions (filtered genes)") +
  theme_portfolio() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_fig(p_box, "02_count_distributions.png", width = 9, height = 7)

# ---- Variance-stabilising transform, PCA ------------------------------------
vsd_all <- vst(dds_all, blind = TRUE)
pca <- prcomp(t(assay(vsd_all)[top_var_rows(assay(vsd_all), 500), ]))
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- cbind(samples, pca$x[, 1:4])
write.csv(data.frame(PC = paste0("PC", 1:length(pve)), percent_variance = pve),
          file.path("results", "qc_pca_variance.csv"), row.names = FALSE)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = treatment, shape = cell_line)) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = treatment_cols) +
  labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("PC2 (%.1f%%)", pve[2]),
       colour = "Treatment", shape = "Donor cell line",
       title = "PCA of all 16 samples (VST, top 500 variable genes)") +
  theme_portfolio()
save_fig(p_pca, "03_pca_all_samples.png", width = 7.5, height = 5.5)

# ---- Sample-to-sample distance heatmap --------------------------------------
d <- as.matrix(dist(t(assay(vsd_all))))
rownames(d) <- colnames(d) <- samples$title
ann <- data.frame(Treatment = samples$treatment, Donor = samples$cell_line,
                  row.names = samples$title)
png(file.path("figures", "04_sample_distances.png"), width = 8, height = 7, units = "in", res = 150)
pheatmap(d, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         annotation_row = ann, annotation_col = ann,
         annotation_colors = list(Treatment = treatment_cols, Donor = donor_cols),
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         main = "Euclidean distances between samples (VST)",
         fontsize = 8, show_colnames = FALSE)
invisible(dev.off())


# ---- Metadata sanity check: glucocorticoid-response score per sample --------
# Mean VST expression (centred per gene) of canonical glucocorticoid targets.
# Every sample annotated as receiving dexamethasone should score high.
genes <- readRDS(file.path("data", "processed", "genes.rds"))
gc_markers <- c("FKBP5", "TSC22D3", "PER1", "KLF15", "ZBTB16", "DUSP1", "CRISPLD2")
gc_ids <- rownames(vsd_all)[genes[rownames(vsd_all), "symbol"] %in% gc_markers]
gm <- assay(vsd_all)[gc_ids, ]
gc_score <- colMeans(gm - rowMeans(gm))
gc_df <- data.frame(title = samples$title, cell_line = samples$cell_line,
                    treatment = samples$treatment, gc_score = round(gc_score, 2))
write.csv(gc_df, file.path("results", "qc_glucocorticoid_score.csv"), row.names = FALSE)
p_gc <- ggplot(gc_df, aes(treatment, gc_score, colour = treatment)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = cell_line), size = 3, colour = "grey30",
                           nudge_x = 0.3, direction = "y", min.segment.length = 0.2) +
  scale_colour_manual(values = treatment_cols, guide = "none") +
  labs(x = "Annotated treatment (GEO metadata)", y = "Glucocorticoid-response score",
       title = "Metadata check: glucocorticoid-response score per sample",
       subtitle = paste("Mean centred VST expression of", paste(gc_markers, collapse = ", "))) +
  theme_portfolio()
save_fig(p_gc, "04b_glucocorticoid_score_qc.png", width = 8.5, height = 5)

saveRDS(keep, file.path("data", "processed", "keep_genes.rds"))
message("QC complete.")
