.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))

start_script("step_13 supplemental_table4_figure3_tNLR")

fw <- load_framework()
candidates <- list(
  neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "CEACAM8", "MMP8", "LCN2", "DEFA4", "ELANE", "CAMP", "LTF"),
  lymphocyte = c("CD3D", "CD3E", "CD2", "TRAC", "IL7R", "LCK", "CD27", "CCR7", "SELL"))
markers <- lapply(candidates, setdiff, y = unique(fw$gene_symbol))
assert(setequal(markers$neutrophil, c("FCGR3B", "CSF3R", "CXCR2", "MMP8")) &&
       setequal(markers$lymphocyte, c("CD2", "TRAC", "CD27")), "off-panel marker set drifted")
record_value("tnlr.markers.neutrophil", length(markers$neutrophil))
record_value("tnlr.markers.lymphocyte", length(markers$lymphocyte))

row_z <- function(x) { z <- t(scale(t(x))); z[!is.finite(z)] <- 0; z }
marker_counts <- list()
build_tnlr <- function(expr, cohort) {
  z <- row_z(expr)
  avail <- lapply(markers, intersect, y = rownames(z))
  # as in the original analysis: each lineage needs at least two measured markers
  # (GPL13667, the GSE65682 platform, has no TRAC probe, so GSE65682 uses 4 + 2)
  assert(all(lengths(avail) >= 2L), cohort, ": fewer than two measured markers for a lineage")
  marker_counts[[cohort]] <<- data.frame(cohort = cohort, neutrophil = length(avail$neutrophil),
                                         lymphocyte = length(avail$lymphocyte),
                                         missing = paste(setdiff(unlist(markers), unlist(avail)), collapse = ","))
  cat(sprintf("   %-10s measured markers neutrophil/lymphocyte = %d/%d%s\n", cohort, length(avail$neutrophil),
              length(avail$lymphocyte), if (nzchar(marker_counts[[cohort]]$missing)) paste0(" (missing ", marker_counts[[cohort]]$missing, ")") else ""))
  neut <- z_safe(colMeans(z[avail$neutrophil, , drop = FALSE]))
  lymph <- z_safe(colMeans(z[avail$lymphocyte, , drop = FALSE]))
  data.frame(sample_id = colnames(expr), tNLR = z_safe(neut - lymph), stringsAsFactors = FALSE)
}
src <- c(GSE65682 = "GSE65682", GSE26440 = "GSE26440", GSE95233 = "GSE95233", TCVGH = "TCVGH_GSE216902")
tnlr <- do.call(rbind, lapply(names(src), function(co) build_tnlr(readRDS(file.path(EXPR_DIR, paste0(src[[co]], ".rds"))), co)))

mc <- do.call(rbind, marker_counts)

record_value("tnlr.cohorts_with_all_seven_markers", sum(mc$neutrophil == 4L & mc$lymphocyte == 3L),
             "GSE65682 (GPL13667) has no TRAC probe")
missing_markers <- mc[nzchar(mc$missing), ]
record_value("tnlr.markers_absent_from_framework", as.numeric(!any(unlist(markers) %in% fw$gene_symbol)))
record_value("tnlr.only_missing_marker_is_TRAC_in_GSE65682",
             as.numeric(nrow(missing_markers) == 1L && missing_markers$cohort == "GSE65682" &&
                        missing_markers$missing == "TRAC"))
record_value("tnlr.GSE65682.lymphocyte_markers", mc$lymphocyte[mc$cohort == "GSE65682"], "CD2 and CD27")

d <- read_csv(file.path(TABLE_DIR, "mortality_analysis_set.csv"))
d <- merge(d, tnlr, by = "sample_id", all.x = TRUE, sort = FALSE)
assert(nrow(d) == 665L && !anyNA(d$tNLR), "tNLR must be available for all 665 patients")

rho <- do.call(rbind, lapply(MORTALITY_COHORTS, function(co) {
  z <- d[d$cohort == co, ]
  data.frame(cohort = co, axis = AXIS_CODES,
             rho = vapply(AXIS_CODES, function(a) cor(z[[a]], z$tNLR, method = "spearman"), numeric(1)))
}))
med_rho <- tapply(rho$rho, rho$axis, median)

or_row <- function(sm, term = "x") {
  b <- sm[term, 1]; se <- sm[term, if ("se(coef)" %in% colnames(sm)) "se(coef)" else "Std. Error"]
  c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = sm[term, "Pr(>|z|)"])
}
out <- list(); cells <- list()
for (a in AXIS_CODES) {
  z <- transform(d, x = d[[a]])
  b0 <- or_row(summary(clogit(death ~ x + strata(cohort), data = z, method = "exact"))$coefficients)
  b1 <- or_row(summary(clogit(death ~ x + tNLR + strata(cohort), data = z, method = "exact"))$coefficients)
  out[[a]] <- data.frame(axis = a, display = AXIS_DISPLAY[[a]], median_rho = med_rho[[a]],
                         base_OR = b0[["OR"]], base_lo = b0[["lo"]], base_hi = b0[["hi"]], base_p = b0[["p"]],
                         tnlr_OR = b1[["OR"]], tnlr_lo = b1[["lo"]], tnlr_hi = b1[["hi"]], tnlr_p = b1[["p"]])
  for (co in MORTALITY_COHORTS) {
    zc <- z[z$cohort == co, ]
    c0 <- or_row(summary(glm(death ~ x, data = zc, family = binomial()))$coefficients)
    c1 <- or_row(summary(glm(death ~ x + tNLR, data = zc, family = binomial()))$coefficients)
    cells[[length(cells) + 1L]] <- data.frame(axis = a, cohort = co, base_OR = c0[["OR"]], base_lo = c0[["lo"]],
      base_hi = c0[["hi"]], tnlr_OR = c1[["OR"]], tnlr_lo = c1[["lo"]], tnlr_hi = c1[["hi"]])
  }
}
st4 <- do.call(rbind, out)
st4$base_q <- p.adjust(st4$base_p, "BH")
st4$tnlr_q <- p.adjust(st4$tnlr_p, "BH")
cells <- do.call(rbind, cells)
cells$direction_preserved <- (cells$base_OR - 1) * (cells$tnlr_OR - 1) > 0
cells$ci_spans_1_before_and_after <- cells$base_lo < 1 & cells$base_hi > 1 & cells$tnlr_lo < 1 & cells$tnlr_hi > 1
write_csv(st4, file.path(TABLE_DIR, "SupplementalTable4_tnlr_adjustment.csv"))

print(st4, row.names = FALSE, digits = 3)

t2 <- read_csv(file.path(TABLE_DIR, "Table2_mortality_associations.csv"))
assert(max(abs(st4$base_OR - t2$pooled_OR[match(st4$axis, t2$axis)])) < 1e-9,
       "unadjusted pooled ORs must reproduce Table 2")
for (i in seq_len(nrow(st4))) {
  a <- st4$axis[i]
  record_value(paste0("st4.median_rho.", a), st4$median_rho[i])
  for (s in c("OR", "lo", "hi", "q")) {
    record_value(sprintf("st4.base.%s.%s", a, s), st4[[paste0("base_", s)]][i])
    record_value(sprintf("st4.tnlr.%s.%s", a, s), st4[[paste0("tnlr_", s)]][i])
  }
}
record_value("tnlr.pooled_directions_preserved", sum((st4$base_OR - 1) * (st4$tnlr_OR - 1) > 0))
record_value("tnlr.same_significant_axes", as.numeric(setequal(st4$axis[st4$base_q <= 0.05], st4$axis[st4$tnlr_q <= 0.05])))
record_value("tnlr.n_significant_after", sum(st4$tnlr_q <= 0.05))
record_value("tnlr.cohort_cells", nrow(cells))
record_value("tnlr.cohort_directions_preserved", sum(cells$direction_preserved))
flips <- cells[!cells$direction_preserved, ]
record_value("tnlr.flips_are_L_GSE26440_L_TCVGH_M_GSE95233",
             as.numeric(setequal(paste(flips$axis, flips$cohort), c("L GSE26440", "L TCVGH", "G GSE95233"))))
record_value("tnlr.flips_ci_span_1_before_and_after", as.numeric(all(flips$ci_spans_1_before_and_after)))

# ---- Supplemental Figure 3 --------------------------------------------------------------------------
long <- rbind(
  data.frame(axis = st4$display, model = "Unadjusted (Table 2)", OR = st4$base_OR, lo = st4$base_lo, hi = st4$base_hi, q = st4$base_q),
  data.frame(axis = st4$display, model = "Adjusted for tNLR", OR = st4$tnlr_OR, lo = st4$tnlr_lo, hi = st4$tnlr_hi, q = st4$tnlr_q))
long$axis <- factor(long$axis, levels = rev(AXIS_DISPLAY[AXIS_CODES]))
long$model <- factor(long$model, levels = c("Unadjusted (Table 2)", "Adjusted for tNLR"))
long$sig <- ifelse(long$q <= 0.05, "q <= 0.05", "q > 0.05")
p <- ggplot(long, aes(x = OR, y = axis, colour = model, shape = sig)) +
  geom_vline(xintercept = 1, linewidth = 0.5, colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = position_dodge(width = 0.62), linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.62), size = 2.7, fill = "white") +
  scale_x_continuous(trans = "log", breaks = c(0.4, 0.6, 0.8, 1, 1.5, 2, 3, 4)) +
  scale_colour_manual(values = c("Unadjusted (Table 2)" = "#4D4D4D", "Adjusted for tNLR" = "#E02930"), name = NULL) +
  scale_shape_manual(values = c("q <= 0.05" = 16, "q > 0.05" = 21), name = NULL) +
  labs(x = "Cohort-stratified pooled odds ratio per SD (95% CI)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 10), legend.position = "bottom")

ggsave(
  filename = file.path(FIG_DIR, "SupplementalFigure3_NLR_adjustment.tiff"),  # 路徑跟副檔名依需求調整
  plot = p,
  width = 6.6,
  height = 4.6,
  units = "in",        # 需確認：原本save_figure用的單位是in還是cm
  dpi = 300,
  compression = "lzw"  # 如果輸出格式是tiff才需要這個參數
)

finish_script()
