args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript prepare_full_universe.R <source_dir> <output_dir>")
source_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(Matrix))

read_source <- function(name) readRDS(file.path(source_dir, name))
parse_key <- function(x) {
  parts <- strsplit(x, ":", fixed = TRUE)
  data.frame(
    chr = vapply(parts, `[`, character(1), 1L),
    start = as.integer(vapply(parts, `[`, character(1), 2L)),
    end = as.integer(vapply(parts, `[`, character(1), 3L)),
    stringsAsFactors = FALSE
  )
}
chr_num <- function(x) as.integer(sub("^chr", "", x))
compute_gene_covariates <- function(meta, genes, window_bp = 1000000) {
  out <- data.frame(
    overlap_gene_count = numeric(nrow(meta)),
    overlap_protein_coding_count = numeric(nrow(meta)),
    overlap_gene_size_mean = numeric(nrow(meta)),
    overlap_gene_size_max = numeric(nrow(meta)),
    local_gene_density_1mb = numeric(nrow(meta)),
    local_protein_coding_density_1mb = numeric(nrow(meta)),
    nearest_gene_distance = numeric(nrow(meta))
  )
  for (chr in unique(meta$Chromosome)) {
    ridx <- which(meta$Chromosome == chr)
    g <- genes[genes$chr == chr, , drop = FALSE]
    if (!nrow(g)) next
    starts <- sort(g$start)
    ends <- sort(g$end)
    pc <- g$gene_type == "protein_coding"
    pc_starts <- sort(g$start[pc])
    pc_ends <- sort(g$end[pc])
    region_start <- meta$j_start[ridx]
    region_end <- meta$j_end[ridx]
    window_start <- pmax(1, region_start - window_bp)
    window_end <- region_end + window_bp
    window_mb <- pmax(1, window_end - window_start + 1) / 1e6
    out$local_gene_density_1mb[ridx] <-
      (findInterval(window_end, starts) - findInterval(window_start - 1, ends)) / window_mb
    if (length(pc_starts)) {
      out$local_protein_coding_density_1mb[ridx] <-
        (findInterval(window_end, pc_starts) - findInterval(window_start - 1, pc_ends)) / window_mb
    }
    for (ii in seq_along(ridx)) {
      candidates <- which(g$start <= region_end[ii] & g$end >= region_start[ii])
      if (length(candidates)) {
        lens <- g$gene_length[candidates]
        out$overlap_gene_count[ridx[ii]] <- length(candidates)
        out$overlap_protein_coding_count[ridx[ii]] <-
          sum(g$gene_type[candidates] == "protein_coding", na.rm = TRUE)
        out$overlap_gene_size_mean[ridx[ii]] <- mean(lens)
        out$overlap_gene_size_max[ridx[ii]] <- max(lens)
      } else {
        distance <- pmin(abs(region_start[ii] - g$end), abs(g$start - region_end[ii]))
        out$nearest_gene_distance[ridx[ii]] <- min(distance, na.rm = TRUE)
      }
    }
    message("  gene covariates: ", chr, " (", length(ridx), " rows)")
  }
  out
}
write_block <- function(df, idx, cols, prefix, block) {
  n <- length(idx)
  i_list <- vector("list", length(cols))
  x_list <- vector("list", length(cols))
  nnz <- integer(length(cols))
  for (j in seq_along(cols)) {
    if (j %% 2000L == 0L) message(block, ": ", j, "/", length(cols))
    v <- suppressWarnings(as.numeric(df[[cols[j]]][idx]))
    nz <- which(is.finite(v) & v != 0)
    i_list[[j]] <- nz
    x_list[[j]] <- v[nz]
    nnz[j] <- length(nz)
  }
  X <- sparseMatrix(
    i = unlist(i_list, use.names = FALSE),
    j = rep.int(seq_along(cols), nnz),
    x = unlist(x_list, use.names = FALSE),
    dims = c(n, length(cols))
  )
  writeMM(X, file.path(output_dir, paste0("X_", block, ".mtx")))
  write.csv(data.frame(feature_index = seq_along(cols) - 1L,
                       feature = paste0(prefix, cols), block = block),
            file.path(output_dir, paste0("features_", block, ".csv")), row.names = FALSE)
  message(block, ": ", nrow(X), " x ", ncol(X), ", nnz=", length(X@x))
  rm(X, i_list, x_list)
  invisible(gc())
}

message("Reading original feature tables")
pfam <- read_source("df_pfam.rds")
go <- read_source("df_go.rds")
exp <- read_source("df_exp.rds")
baseline <- read_source("enhanced_genetic_baseline_covariates.rds")
stopifnot(nrow(pfam) == 95750L, !anyDuplicated(pfam$bed1_id),
          !anyDuplicated(rownames(go)), !anyDuplicated(baseline$bed_id))
exp_region <- parse_key(exp$region_key)
exp_key <- paste0(exp_region$chr, ":", exp_region$start - 1L, ":", exp_region$end + 1L)
stopifnot(!anyDuplicated(exp_key))

pfam_region <- parse_key(pfam$bed1_id)
ord <- order(chr_num(pfam_region$chr), pfam_region$start, pfam_region$end)
keys <- pfam$bed1_id[ord]
idx_go <- match(keys, rownames(go))
idx_exp <- match(keys, exp_key)
idx_base <- match(keys, baseline$bed_id)
stopifnot(!anyNA(idx_base), all(is.finite(pfam$magma[ord])))
go_overlap <- which(!is.na(idx_go))
exp_overlap <- which(!is.na(idx_exp))
stopifnot(isTRUE(all.equal(as.numeric(pfam$magma[ord][go_overlap]),
                           as.numeric(go$magma[idx_go[go_overlap]]), tolerance = 0)))
stopifnot(isTRUE(all.equal(as.numeric(pfam$magma[ord][exp_overlap]),
                           as.numeric(exp$ZSTAT[idx_exp[exp_overlap]]), tolerance = 0)))

rows <- data.frame(
  bed_id = keys,
  Chromosome = as.character(pfam$Chromosome[ord]),
  j_start = as.numeric(pfam$j_start[ord]),
  j_end = as.numeric(pfam$j_end[ord]),
  ZSTAT = as.numeric(pfam$magma[ord]),
  has_expression = !is.na(idx_exp),
  has_go = !is.na(idx_go),
  stringsAsFactors = FALSE
)
for (nm in setdiff(names(baseline), "bed_id")) rows[[nm]] <- baseline[[nm]][idx_base]
missing_gene_cov <- which(is.na(rows$overlap_gene_count))
if (length(missing_gene_cov)) {
  message("Computing missing GENCODE gene covariates for ", length(missing_gene_cov), " rows")
  genes <- read_source("gencode_v44_gene_annotation_cache.rds")
  generated <- compute_gene_covariates(rows[missing_gene_cov, , drop = FALSE], genes)
  for (nm in names(generated)) rows[[nm]][missing_gene_cov] <- generated[[nm]]
  rm(genes, generated)
  invisible(gc())
}
stopifnot(length(unique(rows$Chromosome)) == 22L,
          !anyDuplicated(rows$bed_id),
          sum(rows$ZSTAT > 1.96) == 9136L,
          all(vapply(rows, function(x) !anyNA(x), logical(1))))
write.csv(rows, file.path(output_dir, "rows.csv"), row.names = FALSE)
message("Rows: ", nrow(rows), "; expression present: ", sum(rows$has_expression),
        "; GO present: ", sum(rows$has_go))

write_block(exp, idx_exp,
            setdiff(names(exp), c("region_key", "ZSTAT")), "EXP__", "expression")
write_block(go, idx_go,
            setdiff(names(go), "magma"), "GO__", "go")
write_block(pfam, ord,
            setdiff(names(pfam), c("bed1_id", "Chromosome", "j_start", "j_end", "magma")),
            "PFAM__", "pfam")
message("Export complete")
