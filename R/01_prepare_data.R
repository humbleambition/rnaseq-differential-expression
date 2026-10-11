# 01_prepare_data.R
# Build the count matrix, sample table, and gene annotation.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

raw_dir  <- file.path("data", "raw")
proc_dir <- file.path("data", "processed")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Counts ---------------------------------------------------------------
counts <- read.delim(file.path(raw_dir, "GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
                     check.names = FALSE)
rownames(counts) <- as.character(counts$GeneID)
counts <- as.matrix(counts[, -1])
storage.mode(counts) <- "integer"

# ---- Sample metadata from the GEO series matrix -----------------------------
sm <- readLines(file.path(raw_dir, "GSE52778_series_matrix.txt.gz"))
get_field <- function(prefix) {
  lines <- grep(paste0("^", prefix, "\t"), sm, value = TRUE)
  lapply(lines, function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1]))
}
acc   <- get_field("!Sample_geo_accession")[[1]]
title <- get_field("!Sample_title")[[1]]
chars <- get_field("!Sample_characteristics_ch1")
char_df <- as.data.frame(lapply(chars, function(x) sub("^[^:]+:\\s*", "", x)))
names(char_df) <- vapply(chars, function(x) sub(":.*$", "", x[1]), character(1))

samples <- data.frame(
  sample    = acc,
  title     = title,
  cell_line = char_df[["cell line"]],
  treatment = char_df[["treatment"]],
  ercc_mix  = char_df[["ercc_mix"]],
  stringsAsFactors = FALSE
)
samples$treatment <- factor(samples$treatment,
  levels = c("Untreated", "Dexamethasone", "Albuterol", "Albuterol_Dexamethasone"))
samples$cell_line <- factor(samples$cell_line)
rownames(samples) <- samples$sample

stopifnot(identical(colnames(counts), rownames(samples)))

# ---- Gene annotation (Entrez ID -> symbol, name, biotype) -------------------
ids <- rownames(counts)
genes <- data.frame(
  entrez  = ids,
  symbol  = mapIds(org.Hs.eg.db, ids, "SYMBOL",   "ENTREZID", multiVals = "first"),
  name    = mapIds(org.Hs.eg.db, ids, "GENENAME", "ENTREZID", multiVals = "first"),
  biotype = mapIds(org.Hs.eg.db, ids, "GENETYPE", "ENTREZID", multiVals = "first"),
  stringsAsFactors = FALSE
)
rownames(genes) <- ids

saveRDS(counts,  file.path(proc_dir, "counts.rds"))
saveRDS(samples, file.path(proc_dir, "samples.rds"))
saveRDS(genes,   file.path(proc_dir, "genes.rds"))
write.csv(samples, file.path("results", "sample_table.csv"), row.names = FALSE)

message(sprintf("Count matrix: %d genes x %d samples", nrow(counts), ncol(counts)))
print(table(samples$cell_line, samples$treatment))
