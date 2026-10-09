# =============================================================================
# lib/expression_io.R -- readers that turn the public raw files into gene-level
# expression matrices (Methods, gene expression preprocessing; Supplemental Method 2). Used by step 03.
#
#   Microarray (GEO series matrices): deposited values used as-is; probes mapped
#     to gene symbols with the bundled platform annotation and averaged.
#   RNA-seq counts (recount3, GSE208581, GSE185263, GSE216902): counts summed by
#     gene symbol, genes kept when CPM >= 1 in at least 10% of samples, then
#     edgeR log2-CPM (prior count 1).
#   RNA-seq FPKM (GSE205672): FPKM averaged by symbol, genes kept when > 0 in at
#     least 10% of samples, then log2(FPKM + 1).
# No cross-platform normalisation or batch correction is applied anywhere.
# =============================================================================

standardize_gene_symbols <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "---")] <- NA_character_
  x
}

ensembl_to_symbol <- function(ids) {
  ids <- sub("\\..*$", "", as.character(ids))
  out <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = unique(ids), keytype = "ENSEMBL",
    column = "SYMBOL", multiVals = "first"))
  standardize_gene_symbols(unname(out[ids]))
}

filter_low_counts <- function(counts, min_cpm = 1) {
  min_samples <- max(2, ceiling(ncol(counts) * 0.10))
  counts[rowSums(edgeR::cpm(counts) >= min_cpm) >= min_samples, , drop = FALSE]
}

# ---- GEO series matrix: expression + sample annotation --------------------------------
read_series_matrix <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  assert(length(begin) == 1L && length(end) == 1L && end > begin, "malformed series matrix: ", path)
  expr <- utils::read.delim(textConnection(lines[(begin + 1):(end - 1)]),
                            check.names = FALSE, stringsAsFactors = FALSE)
  rownames(expr) <- expr[[1]]
  expr <- as.matrix(expr[, -1, drop = FALSE])
  storage.mode(expr) <- "numeric"

  meta_lines <- lines[seq_len(begin - 1)]
  sample_meta <- list()
  for (ln in meta_lines[grepl("^!Sample_", meta_lines)]) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    key <- sub("^!Sample_", "", parts[1])
    base_key <- key
    dup_i <- 1
    while (key %in% names(sample_meta)) {
      dup_i <- dup_i + 1
      key <- paste0(base_key, "_", dup_i)
    }
    sample_meta[[key]] <- gsub('^"|"$', "", parts[-1])
  }
  n <- ncol(expr)
  pd <- data.frame(sample_index = seq_len(n), stringsAsFactors = FALSE)
  for (k in names(sample_meta)) {
    vals <- sample_meta[[k]]
    length(vals) <- n
    pd[[make.names(k)]] <- vals
  }
  if (!"geo_accession" %in% names(pd)) pd$geo_accession <- colnames(expr)
  pd <- characteristics_to_columns(pd)
  platform <- unique(gsub('^"|"$', "", unlist(strsplit(
    meta_lines[grepl("^!Series_platform_id", meta_lines)][1], "\t", fixed = TRUE))[-1]))[1]
  list(expr = expr, pheno = pd, platform = platform)
}

characteristics_to_columns <- function(pd) {
  char_cols <- grep("^characteristics", names(pd), ignore.case = TRUE, value = TRUE)
  out <- pd[, setdiff(names(pd), char_cols), drop = FALSE]
  for (cc in char_cols) {
    vals <- as.character(pd[[cc]])
    hit <- grepl(":", vals)
    keys <- trimws(sub(":.*$", "", vals[hit]))
    for (k in unique(keys)) {
      safe <- make.names(k)
      if (!safe %in% names(out)) out[[safe]] <- NA_character_
      idx <- hit & trimws(sub(":.*$", "", vals)) == k
      out[[safe]][idx] <- trimws(sub("^[^:]+:\\s*", "", vals[idx]))
    }
  }
  out
}

collapse_probes_to_genes <- function(expr, platform) {
  map <- read_csv(file.path(PROBE_MAP_DIR, paste0("platform_probe_gene_map_", platform, ".csv")))
  sym <- toupper(trimws(map$gene_symbol[match(rownames(expr), map$probe_id)]))
  keep <- !is.na(sym) & nzchar(sym) & sym != "---" & sym != "NA"
  out <- limma::avereps(expr[keep, , drop = FALSE], ID = sym[keep])
  rownames(out) <- toupper(trimws(rownames(out)))
  out
}

# the GEO sample fields the analyses use (MARS labels, case status, time point, outcome)
GEO_META_FIELDS <- c("dataset", "geo_accession", "title", "disease.state", "outcome",
                     "collection.time", "time.point", "icu_acquired_infection",
                     "mortality_event_28days", "time_to_event_28days",
                     "endotype_class", "survival")

geo_sample_metadata <- function(pheno, gse) {
  pheno$dataset <- gse
  for (f in setdiff(GEO_META_FIELDS, names(pheno))) pheno[[f]] <- NA_character_
  pheno[, GEO_META_FIELDS, drop = FALSE]
}

load_geo_series_cohort <- function(gse) {
  s <- read_series_matrix(file.path(RAW_DIR, "geo_matrix", paste0(gse, "_series_matrix.txt.gz")))
  list(expr_gene = collapse_probes_to_genes(s$expr, s$platform),
       metadata = geo_sample_metadata(s$pheno, gse), platform = s$platform)
}

# ---- recount3 ------------------------------------------------------------------------------
infer_gene_symbols_from_rse <- function(rse) {
  rd <- as.data.frame(SummarizedExperiment::rowData(rse))
  nms <- names(rd)
  symbol_cols <- nms[grepl("symbol|gene_name|gene.name|external_gene_name", nms, ignore.case = TRUE)]
  if (length(symbol_cols)) {
    sym <- standardize_gene_symbols(rd[[symbol_cols[1]]])
    if (sum(!is.na(sym)) >= max(10, nrow(rd) * 0.25)) return(sym)
  }
  id_cols <- nms[grepl("gene_id|geneid|ensembl", nms, ignore.case = TRUE)]
  if (length(id_cols)) {
    sym <- ensembl_to_symbol(rd[[id_cols[1]]])
    if (sum(!is.na(sym)) >= max(10, nrow(rd) * 0.25)) return(sym)
  }
  standardize_gene_symbols(rownames(rse))
}

load_recount3_cohort <- function(project) {
  rse <- readRDS(file.path(RAW_DIR, "recount3", paste0("recount3_", project, "_rse.rds")))
  assay_names <- SummarizedExperiment::assayNames(rse)
  preferred <- c("counts", "raw_counts", "count", "rail_counts", "gene_counts")
  assay_name <- if (any(preferred %in% assay_names)) preferred[preferred %in% assay_names][1] else assay_names[1]
  counts <- as.matrix(SummarizedExperiment::assay(rse, assay_name))
  storage.mode(counts) <- "numeric"
  sym <- infer_gene_symbols_from_rse(rse)
  keep <- !is.na(sym) & nzchar(sym)
  counts_gene <- rowsum(counts[keep, , drop = FALSE], group = sym[keep], reorder = FALSE)
  counts_gene <- filter_low_counts(counts_gene)
  edgeR::cpm(counts_gene, log = TRUE, prior.count = 1)
}

# ---- study-specific GEO supplementary matrices ("ARCHS4" route) ---------------------------
load_counts_by_ensembl <- function(path) {
  dt <- data.table::fread(path)
  gene_id <- dt[[1]]
  counts <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(counts) <- "numeric"
  sym <- ensembl_to_symbol(gene_id)
  keep <- !is.na(sym) & nzchar(sym)
  counts_gene <- rowsum(counts[keep, , drop = FALSE], group = sym[keep], reorder = FALSE)
  edgeR::cpm(filter_low_counts(counts_gene), log = TRUE, prior.count = 1)
}

load_gse208581 <- function() {
  load_counts_by_ensembl(file.path(RAW_DIR, "geo_supplementary", "GSE208581_raw_counts.tsv.gz"))
}

load_gse185263 <- function() {
  load_counts_by_ensembl(file.path(RAW_DIR, "geo_supplementary", "GSE185263_raw_counts.csv.gz"))
}

load_gse205672 <- function() {
  dt <- data.table::fread(file.path(RAW_DIR, "geo_supplementary",
                                    "GSE205672_AllSamples.GeneExpression.FPKM.txt.gz"))
  sample_cols <- setdiff(names(dt), c("gene_id", "SymbolID", "transcript_id(s)"))
  fpkm <- as.matrix(dt[, sample_cols, with = FALSE])
  suppressWarnings(storage.mode(fpkm) <- "numeric")
  fpkm[is.na(fpkm)] <- 0
  sym <- standardize_gene_symbols(dt$SymbolID)
  keep <- !is.na(sym) & nzchar(sym)
  fpkm_gene <- limma::avereps(fpkm[keep, , drop = FALSE], ID = sym[keep])
  fpkm_gene <- fpkm_gene[rowSums(fpkm_gene > 0, na.rm = TRUE) >= max(2, ceiling(ncol(fpkm_gene) * 0.10)), , drop = FALSE]
  expr_gene <- log2(fpkm_gene + 1)
  # keep the samples deposited in the GEO series (460 of the 466 matrix columns)
  sm <- readLines(gzfile(file.path(RAW_DIR, "geo_supplementary", "GSE205672_series_matrix.txt.gz")), warn = FALSE)
  deposited <- gsub('^"|"$', "", strsplit(sm[grepl("^!Sample_description", sm)][1], "\t", fixed = TRUE)[[1]][-1])
  expr_gene[, colnames(expr_gene) %in% deposited, drop = FALSE]
}

# ---- TCVGH / GSE216902 -----------------------------------------------------------------------
load_tcvgh <- function() {
  raw <- utils::read.delim(gzfile(file.path(RAW_DIR, "geo_supplementary",
                                            "GSE216902_All_read_count_matrix.txt.gz")),
                           check.names = FALSE, stringsAsFactors = FALSE)
  counts <- as.matrix(raw[, setdiff(names(raw), "Geneid"), drop = FALSE])
  storage.mode(counts) <- "numeric"
  sym <- ensembl_to_symbol(raw$Geneid)
  keep <- !is.na(sym) & nzchar(sym)
  counts_gene <- rowsum(counts[keep, , drop = FALSE], group = sym[keep], reorder = FALSE)
  edgeR::cpm(filter_low_counts(counts_gene), log = TRUE, prior.count = 1)
}

# patient number, sampling day and deposited condition label for each GSE216902
# sample, from the GEO SOFT file
tcvgh_soft_metadata <- function() {
  x <- readLines(gzfile(file.path(RAW_DIR, "geo_supplementary", "GSE216902_family.soft.gz")), warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", x)
  ends <- c(starts[-1] - 1L, length(x))
  do.call(rbind, lapply(seq_along(starts), function(i) {
    block <- x[starts[[i]]:ends[[i]]]
    title <- sub("^!Sample_title = ", "", grep("^!Sample_title = ", block, value = TRUE)[1])
    chars <- sub("^!Sample_characteristics_ch1 = ", "",
                 grep("^!Sample_characteristics_ch1 = ", block, value = TRUE))
    kv <- setNames(trimws(sub("^[^:]+:", "", chars)), tolower(trimws(sub(":.*$", "", chars))))
    data.frame(sample_id = sub(" \\(PBMC\\)$", "", title),
               patient_number = unname(kv["patient"]),
               visit_code = unname(kv["time"]),
               geo_condition = unname(kv["condition"]),
               stringsAsFactors = FALSE)
  }))
}
