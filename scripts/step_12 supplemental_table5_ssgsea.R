.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))

source(file.path(.script_dir, "lib", "scoring.R"))
start_script("step_12 supplemental_table5_ssgsea")

assert(requireNamespace("GSVA", quietly = TRUE) && "ssgseaParam" %in% getNamespaceExports("GSVA"),
       "GSVA >= 1.99 (with ssgseaParam) is required")
cat(sprintf("   GSVA %s\n", as.character(utils::packageVersion("GSVA"))))

fw <- load_framework()
expr <- readRDS(file.path(EXPR_DIR, "GSE65682.rds"))
gene_sets <- lapply(split(fw$gene_symbol, fw$module), clean_genes)
measured <- lapply(gene_sets, intersect, y = rownames(expr))
n_measured <- vapply(measured, length, integer(1))
gsm <- measured[n_measured >= 3]
srs <- unique(fw$module[fw$module_name == "Scavenger_Resolution_Support"])
record_value("ssgsea.genes", nrow(expr))
record_value("ssgsea.profiles", ncol(expr))
record_value("ssgsea.modules_scorable", length(gsm))
record_value("ssgsea.scavenger_resolution_support.measured_genes", n_measured[[srs]])

args <- list(exprData = expr, geneSets = gsm, minSize = 3, alpha = 0.25, normalize = TRUE, verbose = FALSE)
args <- args[names(args) %in% names(formals(GSVA::ssgseaParam))]
param <- do.call(GSVA::ssgseaParam, args)
ss <- as.matrix(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
module_z <- t(scale(t(ss)))
module_z[!is.finite(module_z)] <- 0
m2a <- tapply(fw$axis, fw$module, function(x) x[[1]])
axis_modules <- split(names(m2a), unname(m2a))
axis_raw <- do.call(rbind, lapply(AXIS_CODES, function(a) {
  colMeans(module_z[intersect(axis_modules[[a]], rownames(module_z)), , drop = FALSE], na.rm = TRUE)
}))
rownames(axis_raw) <- AXIS_CODES
axis_ss <- t(scale(t(axis_raw)))

mo <- load_module_scores(); mo <- mo[mo$dataset == "GSE65682", ]
ax <- load_axis_scores(); ax <- ax[ax$dataset == "GSE65682", ]
samples <- colnames(module_z)
mo <- mo[match(samples, mo$sample_id), ]; ax <- ax[match(samples, ax$sample_id), ]
module_rho <- vapply(rownames(module_z), function(m) cor(module_z[m, ], mo[[m]], method = "spearman"), numeric(1))
axis_rho <- vapply(AXIS_CODES, function(a) cor(axis_ss[a, ], ax[[a]], method = "spearman"), numeric(1))

record_value("ssgsea.module_rho.min", min(module_rho))
record_value("ssgsea.module_rho.median", median(module_rho))
record_value("ssgsea.module_rho.max", max(module_rho))
record_value("ssgsea.axis_rho.lowest_is_DT", as.numeric(names(which.min(axis_rho)) == "T"))
record_value("ssgsea.axis_rho.highest_is_L", as.numeric(names(which.max(axis_rho)) == "L"))

# ---- mortality refit in GSE65682 -------------------------------------------------------------------
d <- read_csv(file.path(TABLE_DIR, "mortality_analysis_set.csv"))
d <- d[d$cohort == "GSE65682", ]
assert(nrow(d) == 479L && sum(d$death) == 114L, "GSE65682 outcome set must be 479 / 114")
idx <- match(d$sample_id, colnames(axis_ss))
assert(!anyNA(idx), "ssGSEA scores missing for outcome patients")
fit <- function(x) {
  sm <- summary(glm(d$death ~ x, family = binomial()))$coefficients
  b <- sm["x", "Estimate"]; se <- sm["x", "Std. Error"]
  c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = sm["x", "Pr(>|z|)"])
}
st5 <- do.call(rbind, lapply(AXIS_CODES, function(a) {
  r <- fit(d[[a]]); s <- fit(axis_ss[a, idx])
  data.frame(axis = a, display = AXIS_DISPLAY[[a]], concordance_rho = axis_rho[[a]],
             rank_OR = r[["OR"]], rank_lo = r[["lo"]], rank_hi = r[["hi"]], rank_p = r[["p"]],
             ss_OR = s[["OR"]], ss_lo = s[["lo"]], ss_hi = s[["hi"]], ss_p = s[["p"]])
}))
st5$rank_q <- p.adjust(st5$rank_p, "BH")
st5$ss_q <- p.adjust(st5$ss_p, "BH")
write_csv(st5, file.path(TABLE_DIR, "SupplementalTable5_ssgsea_vs_rank_mean.csv"))
print(st5, row.names = FALSE, digits = 3)

t2 <- read_csv(file.path(TABLE_DIR, "Table2_mortality_associations.csv"))
assert(max(abs(st5$rank_OR - t2$GSE65682_OR[match(st5$axis, t2$axis)])) < 1e-9,
       "rank-mean ORs must reproduce the GSE65682 column of Table 2")
for (i in seq_len(nrow(st5))) {
  a <- st5$axis[i]
  record_value(paste0("st5.rho.", a), st5$concordance_rho[i])
  for (s in c("OR", "lo", "hi", "q")) {
    record_value(sprintf("st5.rank.%s.%s", a, s), st5[[paste0("rank_", s)]][i])
    record_value(sprintf("st5.ssgsea.%s.%s", a, s), st5[[paste0("ss_", s)]][i])
  }
}
record_value("ssgsea.directions_preserved", sum((st5$rank_OR - 1) * (st5$ss_OR - 1) > 0))
sig_r <- st5$axis[st5$rank_q <= 0.05]; sig_s <- st5$axis[st5$ss_q <= 0.05]
record_value("ssgsea.n_significant_rank", length(sig_r))
record_value("ssgsea.same_significant_axes", as.numeric(setequal(sig_r, sig_s)))
cat("   q <= 0.05 rank-mean: ", paste(AXIS_DISPLAY[sig_r], collapse = ", "),
    " | ssGSEA: ", paste(AXIS_DISPLAY[sig_s], collapse = ", "), "\n", sep = "")
finish_script()
