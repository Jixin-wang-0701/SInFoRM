args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript export_rebuilt_pfam.R rebuilt.rds rows.csv output_dir")
suppressPackageStartupMessages(library(Matrix))
pfam <- readRDS(args[1])
rows <- read.csv(args[2], stringsAsFactors = FALSE)
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(nrow(pfam) == nrow(rows), !anyDuplicated(pfam$bed1_id))
idx <- match(rows$bed_id, pfam$bed1_id)
stopifnot(!anyNA(idx), all(rows$ZSTAT == pfam$magma[idx]))
cols <- 5L:(ncol(pfam) - 1L)
i_list <- vector("list", length(cols))
nnz <- integer(length(cols))
for (j in seq_along(cols)) {
  value <- pfam[[cols[j]]][idx]
  stopifnot(!anyNA(value), all(value %in% 0:1))
  i_list[[j]] <- which(value == 1L)
  nnz[j] <- length(i_list[[j]])
}
X <- sparseMatrix(i = unlist(i_list, use.names = FALSE),
                  j = rep.int(seq_along(cols), nnz),
                  x = rep.int(1, sum(nnz)),
                  dims = c(nrow(rows), length(cols)))
writeMM(X, file.path(output_dir, "X_pfam_rebuilt.mtx"))
header <- strsplit(readLines("df_pfam.csv", n = 1L), ",", fixed = TRUE)[[1]]
old_raw <- header[5L:(length(header) - 1L)]
stopifnot(identical(make.unique(make.names(old_raw)), names(pfam)[cols[seq_along(old_raw)]]))
write.csv(data.frame(feature = names(pfam)[cols],
                     legacy_raw_name = c(old_raw, rep(NA_character_, length(cols) - length(old_raw)))),
          file.path(output_dir, "pfam_feature_map.csv"), row.names = FALSE)
cat("Exported", nrow(X), "x", ncol(X), "with", length(X@x), "nonzero values\n")
