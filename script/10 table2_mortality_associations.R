
.script_dir <- getwd()
source(file.path(.script_dir, "scripts","01 define and dictionary.R"))

start_script("10 table2_mortality_associations")

d <- read_csv(file.path(TABLE_DIR, "mortality_analysis_set.csv"))
assert(nrow(d) == 665L && sum(d$death) == 156L, "mortality set must be 665 patients / 156 deaths")

or_ci <- function(b, se) c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))
fit_cohort <- function(z, a) {
  sm <- summary(glm(death ~ x, data = transform(z, x = z[[a]]), family = binomial()))$coefficients
  c(or_ci(sm["x", "Estimate"], sm["x", "Std. Error"]), p = sm["x", "Pr(>|z|)"])
}
fit_pooled <- function(z, a) {
  sm <- summary(clogit(death ~ x + strata(cohort), data = transform(z, x = z[[a]]), method = "exact"))$coefficients
  c(or_ci(sm["x", "coef"], sm["x", "se(coef)"]), p = sm["x", "Pr(>|z|)"])
}

rows <- list()
for (a in AXIS_CODES) {
  r <- data.frame(axis = a, display = AXIS_DISPLAY[[a]], name = AXIS_NAME[[a]], stringsAsFactors = FALSE)
  for (co in MORTALITY_COHORTS) {
    f <- fit_cohort(d[d$cohort == co, ], a)
    r[[paste0(co, "_OR")]] <- f[["OR"]]; r[[paste0(co, "_lo")]] <- f[["lo"]]; r[[paste0(co, "_hi")]] <- f[["hi"]]
  }
  f <- fit_pooled(d, a)
  r$pooled_OR <- f[["OR"]]; r$pooled_lo <- f[["lo"]]; r$pooled_hi <- f[["hi"]]; r$pooled_p <- f[["p"]]
  rows[[a]] <- r
}
t2 <- do.call(rbind, rows)
t2$pooled_q <- p.adjust(t2$pooled_p, method = "BH")
write_csv(t2, file.path(TABLE_DIR, "Table2_mortality_associations.csv"))

fmt <- function(o, l, h) sprintf("%.2f (%.2f-%.2f)", o, l, h)
show <- data.frame(Axis = paste0(t2$display, ", ", t2$name))
for (co in MORTALITY_COHORTS) show[[co]] <- fmt(t2[[paste0(co, "_OR")]], t2[[paste0(co, "_lo")]], t2[[paste0(co, "_hi")]])
show$Pooled <- fmt(t2$pooled_OR, t2$pooled_lo, t2$pooled_hi)
show$P <- signif(t2$pooled_p, 2); show$q <- signif(t2$pooled_q, 2)
print(show, row.names = FALSE)

for (i in seq_len(nrow(t2))) {
  a <- t2$axis[i]
  for (co in MORTALITY_COHORTS) for (s in c("OR", "lo", "hi"))
    record_value(sprintf("table2.%s.%s.%s", co, a, s), t2[[paste0(co, "_", s)]][i])
  for (s in c("OR", "lo", "hi", "p", "q")) record_value(sprintf("table2.pooled.%s.%s", a, s), t2[[paste0("pooled_", s)]][i])
}
sig <- t2$axis[t2$pooled_q <= 0.05]
record_value("table2.n_significant_axes", length(sig))
record_value("table2.significant_set_is_MD_MP_L_E_M_DT", as.numeric(setequal(sig, c("M", "N_MP", "L", "E", "G", "T"))))
record_value("table2.max_q_significant", max(t2$pooled_q[t2$pooled_q <= 0.05]))
cat("   pooled q <= 0.05: ", paste(AXIS_DISPLAY[sig], collapse = ", "), "\n", sep = "")
finish_script()
