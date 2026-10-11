# install_packages.R -- install the R packages used by this project.
#   Rscript R/install_packages.R
# Tested with R 4.5.0, DESeq2 1.46.0, edgeR 4.4.2, limma 3.62.2.

cran <- c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "matrixStats",
          "rmarkdown", "knitr")
bioc <- c("DESeq2", "edgeR", "limma", "AnnotationDbi", "org.Hs.eg.db", "GO.db")

repos <- "https://cloud.r-project.org"
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = repos)

missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}
message("All packages available.")
