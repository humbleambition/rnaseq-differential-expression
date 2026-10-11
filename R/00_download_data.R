# 00_download_data.R
# Download the public GEO data for GSE52778 (Himes et al. 2014).
# Raw data are not committed to the repository; run this script first.
#
#   - NCBI-generated raw gene counts (GRCh38.p13, NCBI/Entrez gene IDs)
#   - GEO series matrix (sample metadata)

raw_dir <- file.path("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

files <- list(
  counts = list(
    url  = paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts",
                  "&acc=GSE52778&format=file&file=GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
    dest = file.path(raw_dir, "GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  ),
  series_matrix = list(
    url  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/matrix/GSE52778_series_matrix.txt.gz",
    dest = file.path(raw_dir, "GSE52778_series_matrix.txt.gz")
  )
)

options(timeout = 300)
for (f in files) {
  if (file.exists(f$dest)) {
    message("Already present: ", f$dest)
    next
  }
  message("Downloading ", basename(f$dest))
  download.file(f$url, f$dest, mode = "wb", quiet = TRUE)
  # Guard against HTML error / captcha pages being saved instead of gzip data
  magic <- readBin(f$dest, "raw", n = 2)
  if (!identical(magic, as.raw(c(0x1f, 0x8b)))) {
    file.remove(f$dest)
    stop("Download of ", basename(f$dest), " did not return a gzip file. ",
         "Try again later or download manually from ", f$url)
  }
}
message("Done. Files in ", raw_dir, ":")
print(file.info(list.files(raw_dir, full.names = TRUE))[, "size", drop = FALSE])
