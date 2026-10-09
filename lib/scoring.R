# =============================================================================
# lib/scoring.R -- the four-step rank-mean scoring engine (Methods, per-sample scoring and dataset standardisation;
# Supplemental Method 4). This is the only code that turns expression values
# into module and axis scores.
#
#   Step 1  within-sample rank z: rank genes inside each sample, then centre and
#           scale each sample's ranks (independent of scale and of other samples)
#   Step 2  module score = mean rank-z of the module's measured genes; modules
#           with fewer than three measured genes are left unscored
#   Step 3  module scores standardised within the dataset; axis raw score =
#           equal-weight mean of the axis's module z scores
#   Step 4  axis raw scores standardised again within the dataset
# =============================================================================

zscore_rows <- function(x) {
  y <- t(scale(t(x)))
  y[is.na(y)] <- 0
  y
}

zscore_columns_by_group <- function(df, cols, group_col = "dataset") {
  out <- df
  for (idx in split(seq_len(nrow(df)), df[[group_col]])) {
    z <- scale(as.matrix(df[idx, cols, drop = FALSE]))
    z[is.na(z)] <- 0
    out[idx, cols] <- z
  }
  out
}

# gene sets and module -> axis assignment from the frozen framework
framework_sets <- function(fw) {
  gene_sets <- lapply(split(fw$gene_symbol, fw$module), clean_genes)
  module_axis <- tapply(fw$axis, fw$module, function(x) x[[1]])
  module_meta <- data.frame(module = names(module_axis), axis = unname(module_axis),
                            stringsAsFactors = FALSE)
  module_meta <- module_meta[order(match(module_meta$axis, AXIS_STORAGE_ORDER),
                                   module_meta$module), ]
  axis_modules <- split(module_meta$module, module_meta$axis)
  axis_modules <- axis_modules[AXIS_STORAGE_ORDER[AXIS_STORAGE_ORDER %in% names(axis_modules)]]
  list(gene_sets = gene_sets, module_axis = module_axis, axis_modules = axis_modules)
}

# Steps 1-2
rank_mean_modules <- function(expr_gene, gene_sets, min_genes = 3L) {
  rownames(expr_gene) <- toupper(trimws(rownames(expr_gene)))
  gs <- lapply(gene_sets, function(g) intersect(clean_genes(g), rownames(expr_gene)))
  gs <- gs[vapply(gs, length, integer(1)) >= min_genes]
  if (!length(gs)) stop("No gene set has at least ", min_genes, " measured genes", call. = FALSE)
  ranked <- apply(expr_gene, 2, rank, ties.method = "average", na.last = "keep")
  ranked <- scale(ranked)
  out <- do.call(rbind, lapply(gs, function(g) colMeans(ranked[g, , drop = FALSE], na.rm = TRUE)))
  rownames(out) <- names(gs)
  out
}

# Steps 1-4 for one dataset
score_dataset <- function(dataset_id, expr_gene, sets) {
  expr_gene <- as.matrix(expr_gene)
  rownames(expr_gene) <- toupper(trimws(rownames(expr_gene)))
  sc <- rank_mean_modules(expr_gene, sets$gene_sets)
  z <- zscore_rows(sc)
  module_z <- matrix(NA_real_, nrow = length(sets$gene_sets), ncol = ncol(z),
                     dimnames = list(names(sets$gene_sets), colnames(z)))
  module_z[rownames(z), colnames(z)] <- z
  axis_raw <- do.call(rbind, lapply(names(sets$axis_modules), function(a) {
    mods <- intersect(sets$axis_modules[[a]], rownames(module_z))
    colMeans(module_z[mods, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(axis_raw) <- names(sets$axis_modules)
  axis_z <- matrix(NA_real_, nrow = length(AXIS_STORAGE_ORDER), ncol = ncol(axis_raw),
                   dimnames = list(AXIS_STORAGE_ORDER, colnames(axis_raw)))
  axis_z[rownames(axis_raw), colnames(axis_raw)] <- axis_raw
  axis_z <- zscore_rows(axis_z)
  coverage <- do.call(rbind, lapply(names(sets$gene_sets), function(mod) {
    g <- sets$gene_sets[[mod]]
    present <- intersect(g, rownames(expr_gene))
    data.frame(dataset = dataset_id, axis = unname(sets$module_axis[[mod]]), module = mod,
               module_genes = length(g), present_genes = length(present),
               present_gene_symbols = paste(present, collapse = ";"),
               stringsAsFactors = FALSE)
  }))
  list(module_z = module_z, axis_z = axis_z, coverage = coverage)
}

# Steps 1-4 for all datasets; returns long per-sample tables
score_all <- function(expr_list, sets) {
  res <- lapply(names(expr_list), function(ds) {
    cat(sprintf("   scoring %-20s %6d genes x %4d samples\n", ds,
                nrow(expr_list[[ds]]), ncol(expr_list[[ds]])))
    score_dataset(ds, expr_list[[ds]], sets)
  })
  names(res) <- names(expr_list)
  module <- do.call(rbind, lapply(names(res), function(ds) {
    m <- t(res[[ds]]$module_z)
    data.frame(dataset = ds, sample_id = rownames(m), m, check.names = FALSE,
               stringsAsFactors = FALSE)
  }))
  axis <- do.call(rbind, lapply(names(res), function(ds) {
    m <- t(res[[ds]]$axis_z[AXIS_STORAGE_ORDER, , drop = FALSE])
    data.frame(dataset = ds, sample_id = rownames(m), m, check.names = FALSE,
               stringsAsFactors = FALSE)
  }))
  axis <- zscore_columns_by_group(axis, AXIS_STORAGE_ORDER, "dataset")
  coverage <- do.call(rbind, lapply(res, `[[`, "coverage"))
  rownames(module) <- rownames(axis) <- rownames(coverage) <- NULL
  list(module = module, axis = axis, coverage = coverage)
}
