# utils.R -- shared helpers for plotting and output

treatment_cols <- c(
  Untreated               = "#7f7f7f",
  Dexamethasone           = "#d62728",
  Albuterol               = "#1f77b4",
  Albuterol_Dexamethasone = "#9467bd"
)
donor_cols <- c(N052611 = "#1b9e77", N061011 = "#d95f02",
                N080611 = "#7570b3", N61311  = "#e7298a")

theme_portfolio <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
}

save_fig <- function(p, file, width = 7, height = 5, dpi = 150) {
  dir.create("figures", showWarnings = FALSE)
  ggplot2::ggsave(file.path("figures", file), p, width = width, height = height, dpi = dpi)
  invisible(file.path("figures", file))
}

top_var_rows <- function(m, n = 500) {
  v <- matrixStats::rowVars(m)
  order(v, decreasing = TRUE)[seq_len(min(n, length(v)))]
}
