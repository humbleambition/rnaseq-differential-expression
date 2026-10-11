# run_all.R -- run the full pipeline from the project root:
#   Rscript R/run_all.R
steps <- c("00_download_data.R", "01_prepare_data.R", "02_qc.R",
           "03_differential_expression.R", "04_figures_and_enrichment.R")
for (s in steps) {
  message("\n==== ", s, " ====")
  source(file.path("R", s), local = new.env())
}
message("\n==== Rendering report ====")
rmarkdown::render(file.path("report", "analysis.Rmd"), output_format = "all", quiet = TRUE)
writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
message("Pipeline complete.")
