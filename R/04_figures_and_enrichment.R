# 04_figures_and_enrichment.R
# Visualise the Dex vs untreated results and run GO enrichment (limma::goana).

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(limma)
})
source(file.path("R", "utils.R"))

dds   <- readRDS(file.path("results", "dds_dex_vs_untreated.rds"))
de    <- read.csv(file.path("results", "deseq2_dex_vs_untreated_all_genes.csv"),
                  colClasses = c(entrez = "character"))
calls <- readRDS(file.path("results", "method_concordance.rds"))
alpha <- 0.05; lfc_th <- 1

# Canonical glucocorticoid-response genes used as positive controls
gc_markers <- c("FKBP5", "TSC22D3", "PER1", "KLF15", "ZBTB16", "DUSP1", "CRISPLD2")

# ---- PCA of the 8 analysed samples --------------------------------------
vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)[top_var_rows(assay(vsd), 500), ]))
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pdf8 <- cbind(as.data.frame(colData(dds)), pca$x[, 1:2])
p <- ggplot(pdf8, aes(PC1, PC2, colour = treatment, shape = cell_line)) +
  geom_point(size = 4) +
  scale_colour_manual(values = treatment_cols) +
  labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("PC2 (%.1f%%)", pve[2]),
       colour = "Treatment", shape = "Donor cell line",
       title = "PCA: dexamethasone vs untreated (8 samples)") +
  theme_portfolio()
save_fig(p, "05_pca_dex_vs_untreated.png", width = 7, height = 5)

# ---- Volcano plot -----------------------------------------------------------
vd <- subset(de, !is.na(padj))
vd$class <- "Not significant"
vd$class[vd$padj < alpha & vd$log2FC >  lfc_th] <- "Up in Dex"
vd$class[vd$padj < alpha & vd$log2FC < -lfc_th] <- "Down in Dex"
vd$class[vd$padj < alpha & abs(vd$log2FC) <= lfc_th] <- "FDR < 0.05, |log2FC| <= 1"
vd$class <- factor(vd$class, levels = c("Up in Dex", "Down in Dex",
                                        "FDR < 0.05, |log2FC| <= 1", "Not significant"))
vd$neglog10 <- -log10(pmax(vd$padj, 1e-300))
lab <- unique(c(head(vd$symbol[order(vd$padj)], 12),
                head(vd$symbol[vd$padj < alpha][order(-abs(vd$log2FC[vd$padj < alpha]))], 6),
                gc_markers))
vd$label <- ifelse(vd$symbol %in% lab, vd$symbol, NA)
n_up <- sum(vd$class == "Up in Dex"); n_dn <- sum(vd$class == "Down in Dex")
p <- ggplot(vd, aes(log2FC, neglog10, colour = class)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_th, lfc_th), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", linewidth = 0.3) +
  geom_text_repel(aes(label = label), colour = "black", size = 3, max.overlaps = 40,
                  min.segment.length = 0, na.rm = TRUE) +
  scale_colour_manual(values = c("Up in Dex" = "#d62728", "Down in Dex" = "#1f77b4",
                                 "FDR < 0.05, |log2FC| <= 1" = "#bdbdbd",
                                 "Not significant" = "#e5e5e5")) +
  labs(x = "log2 fold change (Dex / untreated)", y = "-log10 adjusted p-value",
       colour = NULL,
       title = "Volcano plot: dexamethasone vs untreated",
       subtitle = sprintf("DESeq2, paired by donor; FDR < %.2f & |log2FC| > %g: %d up, %d down",
                          alpha, lfc_th, n_up, n_dn)) +
  theme_portfolio() + theme(legend.position = "bottom")
save_fig(p, "06_volcano.png", width = 8, height = 6.5)

# ---- MA plot ---------------------------------------------------------------
p <- ggplot(vd, aes(log10(baseMean), log2FC_shrunk, colour = padj < alpha)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_colour_manual(values = c(`FALSE` = "#bdbdbd", `TRUE` = "#d62728"),
                      labels = c("FDR >= 0.05", "FDR < 0.05")) +
  labs(x = "log10 mean normalised count", y = "Shrunken log2 fold change",
       colour = NULL, title = "MA plot (normal-prior LFC shrinkage)") +
  theme_portfolio()
save_fig(p, "07_ma_plot.png", width = 7, height = 5)

# ---- Heatmap of top 40 DE genes ---------------------------------------------
top <- head(subset(de, !is.na(padj) & padj < alpha & abs(log2FC) > lfc_th), 40)
mat <- assay(vsd)[top$entrez, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- ifelse(is.na(top$symbol), top$entrez, top$symbol)
cd <- as.data.frame(colData(dds))
colnames(mat) <- cd$title
ann <- data.frame(Treatment = cd$treatment, Donor = cd$cell_line, row.names = cd$title)
png(file.path("figures", "08_heatmap_top40.png"), width = 7, height = 9, units = "in", res = 150)
pheatmap(mat, annotation_col = ann,
         annotation_colors = list(Treatment = treatment_cols[levels(cd$treatment)],
                                  Donor = donor_cols),
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
         breaks = seq(-max(abs(mat)), max(abs(mat)), length.out = 102),
         main = "Top 40 DE genes (VST, centred per gene)", fontsize_row = 8, fontsize_col = 8)
invisible(dev.off())

# ---- Positive-control genes, per donor ---------------------------------------
nc <- counts(dds, normalized = TRUE)
ids <- de$entrez[match(gc_markers, de$symbol)]
mk <- do.call(rbind, lapply(seq_along(ids), function(i)
  data.frame(gene = gc_markers[i], count = nc[ids[i], ], treatment = cd$treatment,
             donor = cd$cell_line)))
mk$gene <- factor(mk$gene, levels = gc_markers)
p <- ggplot(mk, aes(treatment, count, group = donor, colour = donor)) +
  geom_line(alpha = 0.7) + geom_point(size = 2) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  scale_colour_manual(values = donor_cols) +
  scale_x_discrete(labels = c("Untreated" = "Untr.", "Dexamethasone" = "Dex")) +
  labs(x = NULL, y = "Normalised count (log scale)", colour = "Donor",
       title = "Known glucocorticoid-response genes, paired by donor") +
  theme_portfolio()
save_fig(p, "09_glucocorticoid_markers.png", width = 10, height = 5.5)

# ---- Method concordance -------------------------------------------------------
p <- ggplot(calls, aes(lfc_DESeq2, lfc_edgeR)) +
  geom_point(size = 0.5, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "#d62728", linewidth = 0.4) +
  labs(x = "log2FC (DESeq2)", y = "log2FC (edgeR QL)",
       title = "DESeq2 vs edgeR fold-change estimates",
       subtitle = sprintf("Spearman r = %.3f; significant in all three methods: %d genes",
                          cor(calls$lfc_DESeq2, calls$lfc_edgeR, method = "spearman"),
                          sum(calls$DESeq2 & calls$edgeR & calls$limma_voom))) +
  theme_portfolio()
save_fig(p, "10_method_concordance.png", width = 6, height = 5.5)

# ---- GO enrichment (biological process), up- and down-regulated separately ---
universe <- de$entrez
up <- de$entrez[!is.na(de$padj) & de$padj < alpha & de$log2FC >  lfc_th]
dn <- de$entrez[!is.na(de$padj) & de$padj < alpha & de$log2FC < -lfc_th]
go <- goana(list(Up = up, Down = dn), universe = universe, species = "Hs")
go_bp <- go[go$Ont == "BP", ]
go_bp$FDR.Up   <- p.adjust(go_bp$P.Up,   method = "BH")
go_bp$FDR.Down <- p.adjust(go_bp$P.Down, method = "BH")
go_bp$GO_ID <- rownames(go_bp)
write.csv(go_bp[order(go_bp$P.Up), ], file.path("results", "go_bp_enrichment.csv"), row.names = FALSE)

top_go <- function(dir, n = 10) {
  pcol <- paste0("P.", dir); fcol <- paste0("FDR.", dir); dcol <- dir
  g <- go_bp[go_bp$N >= 15 & go_bp$N <= 500, ]
  g <- head(g[order(g[[pcol]]), ], n)
  data.frame(term = g$Term, direction = dir, DE = g[[dcol]], N = g$N,
             FDR = g[[fcol]], ratio = g[[dcol]] / g$N)
}
tg <- rbind(top_go("Up"), top_go("Down"))
write.csv(tg, file.path("results", "go_bp_top_terms.csv"), row.names = FALSE)
tg$term <- factor(paste0(tg$term, " (", tg$direction, ")"),
                  levels = rev(paste0(tg$term, " (", tg$direction, ")")))
p <- ggplot(tg, aes(-log10(FDR), term, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c(Up = "#d62728", Down = "#1f77b4")) +
  labs(x = "-log10 FDR (BH)", y = NULL, fill = NULL,
       title = "Top GO biological processes (gene sets of 15-500 genes)",
       subtitle = "limma::goana; DE genes with FDR < 0.05 and |log2FC| > 1 vs all tested genes") +
  theme_portfolio() + theme(axis.text.y = element_text(size = 8))
save_fig(p, "11_go_enrichment.png", width = 10, height = 6.5)

message("Figures and enrichment complete.")
