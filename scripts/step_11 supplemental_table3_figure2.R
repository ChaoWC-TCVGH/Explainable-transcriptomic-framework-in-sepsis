.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))

start_script("step_11 supplemental_table3_figure2")

ax <- load_axis_scores()
manifest <- read_csv(file.path(TABLE_DIR, "sample_manifest.csv"))
d <- merge(ax[, c("dataset", "sample_id", AXIS_CODES)], manifest, by = c("dataset", "sample_id"))
d <- d[d$dataset %in% LONGITUDINAL_COHORTS & d$case_status %in% c("sepsis", "septic_shock") &
         is.finite(d$time_from_baseline), ]

paired <- list(); deltas <- list()
for (co in LONGITUDINAL_COHORTS) {
  z <- d[d$dataset == co, ]
  pids <- unique(z$patient_id)
  for (pid in pids) {
    p <- z[z$patient_id == pid, ]
    p <- p[order(p$time_from_baseline, p$sample_id), ]
    for (j in seq_len(nrow(p))) for (a in AXIS_CODES)
      deltas[[length(deltas) + 1L]] <- data.frame(dataset = co, patient_id = pid, time = p$time_from_baseline[j],
                                                   axis = a, delta = p[[a]][j] - p[[a]][1])
  }
  for (a in AXIS_CODES) {
    dd <- unlist(lapply(pids, function(pid) {
      p <- z[z$patient_id == pid, ]
      p <- p[order(p$time_from_baseline, p$sample_id), ]
      if (nrow(p) < 2L) NULL else p[[a]][nrow(p)] - p[[a]][1]
    }))
    exact_ok <- length(dd) > 0L && !any(dd == 0) && !anyDuplicated(abs(dd))
    wt <- suppressWarnings(wilcox.test(dd, mu = 0, exact = exact_ok))
    paired[[length(paired) + 1L]] <- data.frame(dataset = co, axis = a, n_patients = length(dd),
                                                median_delta = median(dd), p = wt$p.value)
  }
}
paired <- do.call(rbind, paired)
paired$q_BH <- ave(paired$p, paired$dataset, FUN = function(x) p.adjust(x, "BH"))
write_csv(paired, file.path(TABLE_DIR, "SupplementalTable3_longitudinal_paired_tests.csv"))
print(paired, row.names = FALSE, digits = 3)

tag <- c(GSE57065 = "GSE57065", GSE95233 = "GSE95233", TCVGH_GSE216902 = "TCVGH")
for (i in seq_len(nrow(paired))) {
  record_value(sprintf("st3.%s.%s.median", tag[[paired$dataset[i]]], paired$axis[i]), paired$median_delta[i])
  record_value(sprintf("st3.%s.%s.q", tag[[paired$dataset[i]]], paired$axis[i]), paired$q_BH[i])
}
for (co in LONGITUDINAL_COHORTS) record_value(paste0("longitudinal.n.", tag[[co]]), paired$n_patients[paired$dataset == co][1])
record_value("longitudinal.n.total", sum(tapply(paired$n_patients, paired$dataset, `[`, 1)))
sig <- function(a) sum(paired$q_BH[paired$axis == a] <= 0.05)
fall <- function(a) sum(paired$median_delta[paired$axis == a] < 0)
rise <- function(a) sum(paired$median_delta[paired$axis == a] > 0)
record_value("longitudinal.metabolic.falls_in_cohorts", fall("G"))
record_value("longitudinal.metabolic.significant_in_cohorts", sig("G"))
record_value("longitudinal.lymphoid.rises_in_cohorts", rise("L"))
record_value("longitudinal.lymphoid.significant_in_cohorts", sig("L"))

# ---- Supplemental Figure 2 (descriptive pooled medians) --------------------------------------
deltas <- do.call(rbind, deltas)
pooled <- aggregate(delta ~ time + axis, deltas, FUN = median)

day_pos <- c(`0` = 0, `1` = 1, `2` = 2, `7` = 4)
day_lab <- c(`0` = "0", `1` = "1", `2` = "2", `7` = "7")
assert(setequal(as.character(unique(pooled$time)), names(day_pos)), "days after first sample must be 0, 1, 2 and 7")
pooled$plot_x <- unname(day_pos[as.character(pooled$time)])
pooled$axis_label <- factor(pooled$axis, levels = AXIS_CODES, labels = AXIS_LABEL[AXIS_CODES])
cols <- setNames(c("#9E6A9E", "#8B6D5C", "#1B7F5A", "#63A859", "#4C78A8", "#F28E2B", "#76B7B2", "#E15759", "#EDC948"),
                 levels(pooled$axis_label))
p <- ggplot(pooled, aes(plot_x, delta, color = axis_label, group = axis_label)) +
  geom_hline(yintercept = 0, color = "#626262", linewidth = 0.55) +
  geom_line(linewidth = 1.05) + geom_point(size = 3.0) +
  scale_color_manual(values = cols) +
  scale_x_continuous(breaks = unname(day_pos), labels = unname(day_lab), limits = c(-0.12, 4.25)) +
  coord_cartesian(ylim = c(-1.0, 0.92), clip = "off") +
  labs(x = "Days after first sample", y = "Median within-patient change\n(dataset z-score)", color = NULL) +
  theme_minimal(base_size = 15) +
  theme(axis.text = element_text(size = 14.2), legend.text = element_text(size = 14.2), legend.position = "bottom",
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  guides(color = guide_legend(ncol = 2, byrow = TRUE))

ggsave(
  filename = file.path(FIG_DIR, "SupplementalFigure2.tiff"),  # 路徑跟副檔名依需求調整
  plot = p,
  width = 8.8,
  height = 6.5875,
  units = "in",        # 需確認：原本save_figure用的單位是in還是cm
  dpi = 300,
  compression = "lzw"  # 如果輸出格式是tiff才需要這個參數
)
finish_script()
