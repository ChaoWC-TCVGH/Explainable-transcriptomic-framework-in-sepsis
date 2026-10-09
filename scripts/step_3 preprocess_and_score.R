library("limma")
library("edgeR")
library("data.table")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("SummarizedExperiment")

.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))
source(file.path(.script_dir, "lib", "scoring.R"))
source(file.path(.script_dir, "lib", "expression_io.R"))


dir.create(EXPR_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RECOMPUTED_CHECKPOINTS, recursive = TRUE, showWarnings = FALSE)

fw <- load_framework()
sets <- framework_sets(fw)

expr <- list()
geo_meta <- list()
for (gse in GEO_SERIES_COHORTS) {
  cat(sprintf("   reading %s series matrix\n", gse))
  x <- load_geo_series_cohort(gse)
  expr[[gse]] <- x$expr_gene
  geo_meta[[gse]] <- x$metadata
  cat(sprintf("      platform %s -> %d genes x %d samples\n", x$platform, nrow(x$expr_gene), ncol(x$expr_gene)))
}
for (p in c("SRP049820", "SRP050000", "SRP212956")) {
  expr[[paste0("recount3_", p)]] <- load_recount3_cohort(p)
}
expr[["ARCHS4_GSE185263"]] <- load_gse185263()
expr[["ARCHS4_GSE205672"]] <- load_gse205672()
expr[["ARCHS4_GSE208581"]] <- load_gse208581()
expr[["TCVGH_GSE216902"]] <- load_tcvgh()

assert(setequal(names(expr), COHORTS$dataset), "the 12 cohorts were not all loaded")
n_profiles <- sum(vapply(expr, ncol, integer(1)))
cat(sprintf("   %d cohorts, %d profiles\n", length(expr), n_profiles))
assert(n_profiles == 3083L, "expected 3,083 deposited profiles, found ", n_profiles)
for (ds in names(expr)) saveRDS(expr[[ds]], file.path(EXPR_DIR, paste0(ds, ".rds")))

scores <- score_all(expr, sets)

gz_write <- function(x, path) {
  con <- gzfile(path, "wt")
  utils::write.csv(x, con, row.names = FALSE)
  close(con)
}



gz_write(scores$axis[, c("dataset", "sample_id", AXIS_STORAGE_ORDER)], RE_AXIS)
gz_write(scores$module, RE_MODULE)
gz_write(scores$coverage, RE_COVERAGE)
gz_write(do.call(rbind, geo_meta), RE_GEO_META)

tm <- tcvgh_soft_metadata()
tm <- tm[match(colnames(expr[["TCVGH_GSE216902"]]), tm$sample_id), ]

source(TCVGH_OUTCOME_R)
label <- unique(tm[, c("patient_number", "geo_condition")])

write_csv(tm[, c("sample_id", "patient_number", "visit_code")], RE_TCVGH_META)

cat("   checkpoints written to ", RECOMPUTED_CHECKPOINTS, "\n", sep = "")
