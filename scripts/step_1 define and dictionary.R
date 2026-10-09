library(survival)
library(ggplot2)
library(patchwork)
library(digest)
library(data.table)
library(limma)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(recount3)
library(GSVA)
library(BiocParallel)
library(randomForest)
library(sva)
#BiocManager::install(c("hgu219.db", "hgu133plus2.db"))
library(hgu133plus2.db)
library(hgu219.db)
library(edgeR)
library(readr)


CODE_ROOT <- "D:/Gene_analysis/test"

env_or <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

DATA_DIR       <- file.path(CODE_ROOT, "data")
OUT_DIR        <- env_or("SEPSIS_OUT_DIR", file.path(CODE_ROOT, "outputs"))
WORK_DIR       <- env_or("SEPSIS_WORK_DIR", file.path(CODE_ROOT, "work"))
RAW_DIR        <- env_or("SEPSIS_RAW_DIR", file.path(WORK_DIR, "raw"))

# the scores recomputed in step 03 if they are there, otherwise the copy shipped in data/
RECOMPUTED_CHECKPOINTS <- file.path(WORK_DIR, "checkpoints")
CHECKPOINT_DIR <- env_or("SEPSIS_CHECKPOINT_DIR",
                         if (file.exists(file.path(RECOMPUTED_CHECKPOINTS, "axis_scores.csv.gz")))
                           RECOMPUTED_CHECKPOINTS else file.path(DATA_DIR, "checkpoints"))
SCORES_SOURCE  <- if (startsWith(CHECKPOINT_DIR, RECOMPUTED_CHECKPOINTS)) {
  "recomputed in step 03"
} else {
  "shipped in data/checkpoints"
}
EXPR_DIR       <- file.path(WORK_DIR, "expression")

TABLE_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")
VALUE_DIR <- file.path(OUT_DIR, "values")
LOG_DIR   <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, TABLE_DIR, FIG_DIR, VALUE_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- bundled inputs -------------------------------------------------------------
FRAMEWORK_CSV       <- file.path(DATA_DIR, "framework", "framework_membership.csv")
PROBE_MAP_DIR       <- file.path(DATA_DIR, "annotation")
CTS_LABEL_DIR       <- file.path(DATA_DIR, "labels")
TCVGH_OUTCOME_R     <- file.path(CODE_ROOT, "lib", "tcvgh_in_hospital_mortality.R")
PUBLIC_MANIFEST_CSV <- file.path(DATA_DIR, "public_data_manifest.csv")

##----Caculate rank score 計算分數-----------------------------------------------------------------
RE_AXIS       <- file.path(RECOMPUTED_CHECKPOINTS, "axis_scores.csv.gz")
RE_MODULE     <- file.path(RECOMPUTED_CHECKPOINTS, "module_scores.csv.gz")
RE_COVERAGE   <- file.path(RECOMPUTED_CHECKPOINTS, "module_gene_coverage.csv.gz")
RE_GEO_META   <- file.path(RECOMPUTED_CHECKPOINTS, "geo_sample_metadata.csv.gz")
RE_TCVGH_META <- file.path(RECOMPUTED_CHECKPOINTS, "tcvgh_sample_metadata.csv")

# ---- analysis constants -----------------------------------------------------------
# Stored axis codes are the column names in every score file; display codes are
# what the tables and figures print (stored M = M-D, N_MP = M-P, G = M [metabolic],
# T = D-T). The vector order is the order used in the tables and figures.
AXIS_CODES <- c("M", "N_MP", "L", "I", "A", "E", "H", "G", "T")
AXIS_STORAGE_ORDER <- c("A", "E", "G", "H", "I", "L", "M", "N_MP", "T")
AXIS_DISPLAY <- c(M = "M-D", N_MP = "M-P", L = "L", I = "I", A = "A",
                  E = "E", H = "H", G = "M", T = "D-T")
AXIS_NAME <- c(M = "myeloid-detrimental", N_MP = "myeloid-protective",
               L = "lymphoid", I = "interferon", A = "cytokine/chemokine",
               E = "endothelial", H = "heme/coagulation", G = "metabolic",
               T = "damage-tolerance")
AXIS_LABEL <- setNames(paste0(AXIS_DISPLAY[AXIS_CODES], ", ", AXIS_NAME[AXIS_CODES]),
                       AXIS_CODES)

# The 12 cohorts, keyed by the dataset identifier used in the score files.
COHORTS <- data.frame(
  dataset = c("GSE65682", "ARCHS4_GSE205672", "ARCHS4_GSE185263", "GSE134347",
              "ARCHS4_GSE208581", "recount3_SRP212956", "GSE26440",
              "recount3_SRP050000", "GSE95233", "GSE57065",
              "recount3_SRP049820", "TCVGH_GSE216902"),
  accession = c("GSE65682", "GSE205672", "GSE185263", "GSE134347", "GSE208581",
                "SRP212956", "GSE26440", "SRP050000", "GSE95233", "GSE57065",
                "GSE60424 / SRP049820", "GSE216902 (TCVGH)"),
  route = c("GEO", "ARCHS4", "ARCHS4", "GEO", "ARCHS4", "recount3", "GEO",
            "recount3", "GEO", "GEO", "recount3", "GEO"),
  platform = c("Microarray (GPL13667)", "RNA-seq", "RNA-seq", "Microarray (GPL17586)",
               "RNA-seq", "RNA-seq", "Microarray (GPL570)", "RNA-seq",
               "Microarray (GPL570)", "Microarray (GPL570)", "RNA-seq", "RNA-seq"),
  stringsAsFactors = FALSE
)
GEO_SERIES_COHORTS   <- c("GSE134347", "GSE57065", "GSE95233", "GSE26440", "GSE65682")
MORTALITY_COHORTS    <- c("GSE65682", "GSE26440", "GSE95233", "TCVGH")
LONGITUDINAL_COHORTS <- c("GSE57065", "GSE95233", "TCVGH_GSE216902")
CTS_COHORTS          <- c("GSE65682", "GSE95233", "GSE134347", "GSE26440", "GSE57065")
SOM_COHORTS          <- c("GSE65682", "GSE95233")
SOM_MODULES <- c("SoM_M1_detrimental_acute_myeloid_metabolic",
                 "SoM_M2_detrimental_granulopoiesis_neutrophil",
                 "SoM_M3_protective_monocyte_IFN_autophagy",
                 "SoM_M4_protective_antigen_lymphoid")

# ---- small helpers ------------------------------------------------------------------
read_csv <- function(path, ...) {
  if (!file.exists(path)) stop("Required input not found: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}
write_csv <- function(x, path, ...) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", ...)
}
clean_genes <- function(x) {
  x <- toupper(trimws(x))
  unique(x[!is.na(x) & nzchar(x)])
}
z_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  y <- (x - mean(x, na.rm = TRUE)) / s
  y[!is.finite(y)] <- NA_real_
  y
}
assert <- function(cond, ...) {
  if (!isTRUE(cond)) stop("ASSERTION FAILED: ", paste0(...), call. = FALSE)
  invisible(TRUE)
}
assert_unique_key <- function(x, cols, label) {
  key <- do.call(paste, c(x[cols], sep = "\r"))
  assert(!anyDuplicated(key), label, " has a duplicated key (", paste(cols, collapse = "+"), ")")
}
sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

# ---- per-script logging and value recording ----------------------------------------
# Each step calls start_script() first and finish_script() last. The key numbers a step computes
# are recorded with record_value() and written to outputs/values/<step>.csv, so that any printed
# result can be traced back to the step that produced it.
.values <- new.env(parent = emptyenv())
.values$rows <- list()
.script_name <- NULL
.log_con <- NULL

start_script <- function(name) {
  .script_name <<- name
  .log_con <<- file(file.path(LOG_DIR, paste0(name, ".log")), open = "wt", encoding = "UTF-8")
  sink(.log_con, split = TRUE)
  cat(sprintf("== %s | %s\n", name, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat("   scores: ", SCORES_SOURCE, " (", CHECKPOINT_DIR, ")\n", sep = "")
  invisible(TRUE)
}

record_value <- function(id, value, note = "") {
  assert(length(value) == 1L, "record_value(", id, ") needs a single value")
  .values$rows[[id]] <- data.frame(id = id, value = as.numeric(value),
                                   script = .script_name, note = note,
                                   stringsAsFactors = FALSE)
  invisible(value)
}

finish_script <- function() {
  if (length(.values$rows)) {
    v <- do.call(rbind, .values$rows)
    write_csv(v, file.path(VALUE_DIR, paste0(.script_name, ".csv")))
    cat(sprintf("   recorded %d key numbers\n", nrow(v)))
  }
  cat(sprintf("%s PASS\n", .script_name))
  sink()
  close(.log_con)
  invisible(TRUE)
}

# ---- figure devices that work on Windows, macOS and headless Linux -------------------
save_figure <- function(plot, stem, width, height, dpi = 300) {
  png_path <- file.path(FIG_DIR, paste0(stem, ".png"))
  pdf_path <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi,
                    bg = "white", device = ragg::agg_png, limitsize = FALSE)
  } else {
    ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi,
                    bg = "white", limitsize = FALSE)
  }
  pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, bg = "white",
                  device = pdf_dev, limitsize = FALSE)
  cat(sprintf("   figure written: %s (.png, .pdf)\n", stem))
  invisible(png_path)
}

# ---- shared loaders -------------------------------------------------------------------
load_framework <- function() {
  fw <- read_csv(FRAMEWORK_CSV)
  fw$gene_symbol <- toupper(trimws(fw$gene_symbol))
  assert(nrow(fw) == 881L, "framework must contain 881 gene-module memberships")
  assert(length(unique(fw$gene_symbol)) == 725L, "framework must contain 725 unique genes")
  assert(length(unique(fw$module)) == 51L, "framework must contain 51 modules")
  fw
}

load_axis_scores <- function() {
  ax <- read_csv(RE_AXIS)
  assert(nrow(ax) == 3083L, "axis score checkpoint must have 3,083 profiles")
  assert_unique_key(ax, c("dataset", "sample_id"), "axis scores")
  ax
}

load_module_scores <- function() {
  mo <- read_csv(RE_MODULE)
  assert(nrow(mo) == 3083L, "module score checkpoint must have 3,083 profiles")
  mo
}

load_cts_labels <- function(cohorts = CTS_COHORTS) {
  do.call(rbind, lapply(cohorts, function(co) {
    x <- read_csv(file.path(CTS_LABEL_DIR, paste0("cts_labels_", co, ".csv")))
    x$cohort <- co
    x[, c("Sample", "CTS", "cohort")]
  }))
}

invisible(TRUE)