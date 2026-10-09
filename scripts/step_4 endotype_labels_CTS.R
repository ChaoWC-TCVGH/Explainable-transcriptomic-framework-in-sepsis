library(sva)
library(randomForest)
library(AnnotationDbi)
library(hgu133plus2.db)
library(hgu219.db)

.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))
source(file.path(.script_dir, "lib", "expression_io.R"))
##自行安裝
#remotes::install_github("bpscicluna/ConsensusTranscriptomicSubtype")

ref <- new.env()
ref_dir <- Sys.getenv("SEPSIS_CTS_REFERENCE_DIR", unset = "")
if (nzchar(ref_dir)) {
  load(file.path(ref_dir, "exp_core_g.rda"), envir = ref)
  load(file.path(ref_dir, "core_samples.rda"), envir = ref)
} else {
  assert(requireNamespace("ConsensusTranscriptomicSubtype", quietly = TRUE),
         "install ConsensusTranscriptomicSubtype or set SEPSIS_CTS_REFERENCE_DIR")
  data("exp_core_g", package = "ConsensusTranscriptomicSubtype")
  assign("exp_core_g", exp_core_g, envir = ref)
  rm(exp_core_g)   # 清掉全域環境的暫存版本，避免混淆
  
  data("core_samples", package = "ConsensusTranscriptomicSubtype")
  assign("core_samples", core_samples, envir = ref)
  rm(core_samples)
}
cat(sprintf("   reference: %d genes x %d samples; randomForest %s; sva %s\n", nrow(ref$exp_core_g),
            ncol(ref$exp_core_g), utils::packageVersion("randomForest"), utils::packageVersion("sva")))

collapse_mean <- function(mat, id) {
  keep <- !is.na(id) & nzchar(id)
  n <- table(id[keep])
  s <- rowsum(mat[keep, , drop = FALSE], group = id[keep], reorder = FALSE)
  sweep(s, 1, as.numeric(n[rownames(s)]), FUN = "/")
}
ensembl_expression <- function(gse) {
  s <- read_series_matrix(file.path(RAW_DIR, "geo_matrix", paste0(gse, "_series_matrix.txt.gz")))
  e <- s$expr
  if (s$platform %in% c("GPL570", "GPL13667")) {
    db <- if (s$platform == "GPL570") hgu133plus2.db::hgu133plus2.db else hgu219.db::hgu219.db
    ens <- suppressMessages(AnnotationDbi::mapIds(db, keys = rownames(e), column = "ENSEMBL",
                                                  keytype = "PROBEID", multiVals = "first"))
    return(collapse_mean(e, ens[rownames(e)]))
  }
  # GPL17586: bundled probe -> symbol map, then symbol -> Ensembl for the 18 CTS genes
  map <- read_csv(file.path(PROBE_MAP_DIR, "platform_probe_gene_map_GPL17586.csv"))
  sym <- collapse_mean(e, setNames(map$gene_symbol, map$probe_id)[rownames(e)])
  lookup <- read_csv(file.path(CTS_LABEL_DIR, "cts18_symbol_to_ensembl.csv"))
  sym <- sym[rownames(sym) %in% lookup$symbol, , drop = FALSE]
  rownames(sym) <- setNames(lookup$ensembl, lookup$symbol)[rownames(sym)]
  sym
}
cts_classify <- function(x, seed) {
  rownames(x) <- sub("\\..*$", "", rownames(x))
  n <- table(rownames(x))
  x <- rowsum(x, group = rownames(x), reorder = FALSE)
  x <- sweep(x, 1, as.numeric(n[rownames(x)]), FUN = "/")
  common <- intersect(rownames(ref$exp_core_g), rownames(x))
  assert(length(common) / nrow(ref$exp_core_g) >= 0.8, "CTS reference gene overlap below 80%")
  rm <- ref$exp_core_g[common, , drop = FALSE]
  nm <- x[common, , drop = FALSE]
  lab <- ref$core_samples[colnames(rm), "CTS", drop = TRUE]
  corrected <- suppressMessages(sva::ComBat(dat = cbind(rm, nm), batch = c(rep("reference", ncol(rm)), rep("new", ncol(nm))),
                                            par.prior = TRUE, prior.plots = FALSE))
  set.seed(seed)
  rf <- randomForest::randomForest(x = t(corrected[, colnames(rm), drop = FALSE]), y = factor(lab),
                                   ntree = 500, importance = TRUE)
  prob <- predict(rf, newdata = t(corrected[, colnames(nm), drop = FALSE]), type = "prob")
  raw <- colnames(prob)[max.col(prob, ties.method = "first")]
  data.frame(Sample = rownames(prob), CTS = ifelse(grepl("^CTS", raw), raw, paste0("CTS", raw)),
             Max_probability = apply(prob, 1, max), stringsAsFactors = FALSE)
}

# seeds as used when the labels were generated
seeds <- c(GSE65682 = 20260702, GSE95233 = 20260702, GSE134347 = 20260706, GSE26440 = 20260706, GSE57065 = 20260706)
res <- do.call(rbind, lapply(names(seeds), function(co) {
  new <- cts_classify(ensembl_expression(co), seeds[[co]])
  old <- read_csv(file.path(CTS_LABEL_DIR, paste0("cts_labels_", co, ".csv")))
  j <- merge(old[, c("Sample", "CTS", "Max_probability")], new, by = "Sample", suffixes = c("_bundled", "_regenerated"))
  r <- data.frame(cohort = co, n = nrow(old), joined = nrow(j), identical_labels = sum(j$CTS_bundled == j$CTS_regenerated),
                  max_abs_probability_difference = max(abs(j$Max_probability_bundled - j$Max_probability_regenerated)))
  cat(sprintf("   %-10s %4d samples, %4d identical labels, max |dP| %.3g\n", co, r$n, r$identical_labels,
              r$max_abs_probability_difference))
  r
}))
write_csv(res, file.path(TABLE_DIR, "cts_label_check.csv"))
assert(all(res$joined == res$n & res$identical_labels == res$n), "regenerated CTS labels differ from the bundled labels")

