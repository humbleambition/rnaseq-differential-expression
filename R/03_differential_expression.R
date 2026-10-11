# 03_differential_expression.R
# Dexamethasone vs untreated, paired by donor cell line (8 samples).
# Primary analysis: DESeq2. Cross-checks: edgeR quasi-likelihood and limma-voom.

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
})

counts  <- readRDS(file.path("data", "processed", "counts.rds"))
samples <- readRDS(file.path("data", "processed", "samples.rds"))
genes   <- readRDS(file.path("data", "processed", "genes.rds"))

alpha  <- 0.05   # FDR threshold
lfc_th <- 1      # |log2 fold change| threshold used for "strong" DE calls

# ---- Subset to the contrast of interest -------------------------------------
sel <- samples$treatment %in% c("Untreated", "Dexamethasone")
s8  <- droplevels(samples[sel, ])
c8  <- counts[, sel]
keep <- rowSums(c8 >= 10) >= 4
c8  <- c8[keep, ]
message(sprintf("Genes tested: %d", nrow(c8)))

# ---- DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(c8, s8, design = ~ cell_line + treatment)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("treatment", "Dexamethasone", "Untreated"), alpha = alpha)
# Shrunken LFCs for ranking/visualisation ("normal" prior; apeglm/ashr not required)
res_shr <- lfcShrink(dds, contrast = c("treatment", "Dexamethasone", "Untreated"),
                     res = res, type = "normal", quiet = TRUE)

de <- data.frame(
  entrez         = rownames(res),
  symbol         = genes[rownames(res), "symbol"],
  name           = genes[rownames(res), "name"],
  biotype        = genes[rownames(res), "biotype"],
  baseMean       = res$baseMean,
  log2FC         = res$log2FoldChange,
  log2FC_shrunk  = res_shr$log2FoldChange,
  lfcSE          = res$lfcSE,
  stat           = res$stat,
  pvalue         = res$pvalue,
  padj           = res$padj,
  stringsAsFactors = FALSE
)
de <- de[order(de$padj, -abs(de$log2FC)), ]
write.csv(de, file.path("results", "deseq2_dex_vs_untreated_all_genes.csv"), row.names = FALSE)
sig <- subset(de, !is.na(padj) & padj < alpha)
write.csv(sig, file.path("results", "deseq2_dex_vs_untreated_significant.csv"), row.names = FALSE)

# ---- edgeR (quasi-likelihood) -----------------------------------------------
y <- DGEList(c8, group = s8$treatment)
y <- normLibSizes(y)
design <- model.matrix(~ cell_line + treatment, data = s8)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "treatmentDexamethasone")
tt_edger <- topTags(qlf, n = Inf)$table

# ---- limma-voom -------------------------------------------------------------
v <- voom(y, design)
vfit <- eBayes(lmFit(v, design))
tt_voom <- topTable(vfit, coef = "treatmentDexamethasone", number = Inf)

# ---- Concordance ------------------------------------------------------------
ids <- rownames(c8)
calls <- data.frame(
  entrez = ids,
  DESeq2 = !is.na(res[ids, "padj"]) & res[ids, "padj"] < alpha,
  edgeR  = tt_edger[ids, "FDR"] < alpha,
  limma_voom = tt_voom[ids, "adj.P.Val"] < alpha,
  lfc_DESeq2 = res[ids, "log2FoldChange"],
  lfc_edgeR  = tt_edger[ids, "logFC"],
  lfc_voom   = tt_voom[ids, "logFC"]
)
saveRDS(calls, file.path("results", "method_concordance.rds"))

summary_tbl <- data.frame(
  metric = c("Genes tested",
             sprintf("DESeq2 FDR < %.2f", alpha),
             sprintf("DESeq2 FDR < %.2f, up in Dex", alpha),
             sprintf("DESeq2 FDR < %.2f, down in Dex", alpha),
             sprintf("DESeq2 FDR < %.2f and |log2FC| > %g", alpha, lfc_th),
             sprintf("  of which up (log2FC > %g)", lfc_th),
             sprintf("  of which down (log2FC < -%g)", lfc_th),
             sprintf("edgeR QL FDR < %.2f", alpha),
             sprintf("limma-voom FDR < %.2f", alpha),
             "Significant in all three methods",
             "Spearman r of log2FC, DESeq2 vs edgeR",
             "Spearman r of log2FC, DESeq2 vs limma-voom"),
  value = c(nrow(c8),
            nrow(sig),
            sum(sig$log2FC > 0),
            sum(sig$log2FC < 0),
            sum(abs(sig$log2FC) > lfc_th),
            sum(sig$log2FC > lfc_th),
            sum(sig$log2FC < -lfc_th),
            sum(calls$edgeR),
            sum(calls$limma_voom),
            sum(calls$DESeq2 & calls$edgeR & calls$limma_voom),
            round(cor(calls$lfc_DESeq2, calls$lfc_edgeR, method = "spearman"), 3),
            round(cor(calls$lfc_DESeq2, calls$lfc_voom,  method = "spearman"), 3))
)
write.csv(summary_tbl, file.path("results", "de_summary.csv"), row.names = FALSE)
print(summary_tbl, row.names = FALSE)

saveRDS(dds, file.path("results", "dds_dex_vs_untreated.rds"))
message("Differential expression complete.")
