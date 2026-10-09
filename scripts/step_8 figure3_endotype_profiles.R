
.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))
start_script("step_8 figure3_endotype_profiles")

ax <- load_axis_scores()
mo <- load_module_scores()
manifest <- read_csv(file.path(TABLE_DIR, "sample_manifest.csv"))
meta <- read_csv(RE_GEO_META)
cts <- load_cts_labels()

scores <- merge(ax[, c("dataset", "sample_id", AXIS_CODES)], manifest, by = c("dataset", "sample_id"))
sepsis <- scores[!(scores$case_status %in% c("control", "nonseptic")), ]

cohort_contrasts <- function(dat, classifier, group_col, reference) {
  rows <- list()
  groups <- unique(dat[[group_col]])
  groups <- groups[!is.na(groups) & nzchar(groups)]
  for (co in unique(dat$dataset)) {
    d <- dat[dat$dataset == co, ]
    if (!reference %in% d[[group_col]]) next
    for (g in groups) {
      if (!g %in% d[[group_col]]) next
      for (a in AXIS_CODES) {
        rows[[length(rows) + 1L]] <- data.frame(
          classifier = classifier, dataset = co, group = g, reference = reference, axis = a,
          n_group = sum(d[[group_col]] == g), n_reference = sum(d[[group_col]] == reference),
          contrast = mean(d[[a]][d[[group_col]] == g], na.rm = TRUE) -
            mean(d[[a]][d[[group_col]] == reference], na.rm = TRUE), stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

mars <- merge(sepsis[sepsis$dataset == "GSE65682", ],
              meta[meta$dataset == "GSE65682", c("geo_accession", "endotype_class")],
              by.x = "sample_id", by.y = "geo_accession")
mars <- mars[mars$endotype_class %in% paste0("Mars", 1:4), ]
cts_dat <- merge(sepsis, cts, by.x = c("dataset", "sample_id"), by.y = c("cohort", "Sample"))
cts_dat <- cts_dat[cts_dat$CTS %in% paste0("CTS", 1:3), ]
som <- mo[mo$dataset %in% SOM_COHORTS, c("dataset", "sample_id", SOM_MODULES)]
som$som_severity <- rowMeans(som[, SOM_MODULES[1:2]]) - rowMeans(som[, SOM_MODULES[3:4]])
som_dat <- merge(sepsis[sepsis$dataset %in% SOM_COHORTS, ], som[, c("dataset", "sample_id", "som_severity")],
                 by = c("dataset", "sample_id"))
som_dat$SoM <- NA_character_
for (co in unique(som_dat$dataset)) {
  ii <- som_dat$dataset == co
  br <- unique(quantile(som_dat$som_severity[ii], c(0, 1/3, 2/3, 1), na.rm = TRUE))
  assert(length(br) == 4L, "non-unique SoM tertile breaks in ", co)
  som_dat$SoM[ii] <- as.character(cut(som_dat$som_severity[ii], breaks = br, include.lowest = TRUE,
                                      labels = c("SoM-Low", "SoM-Mid", "SoM-High")))
}

all_c <- rbind(cohort_contrasts(mars, "MARS", "endotype_class", "Mars2"),
               cohort_contrasts(cts_dat, "CTS", "CTS", "CTS1"),
               cohort_contrasts(som_dat, "SoM", "SoM", "SoM-Low"))
agg <- aggregate(contrast ~ classifier + group + reference + axis, all_c, FUN = mean)
names(agg)[names(agg) == "contrast"] <- "cohort_balanced_contrast"


record_value("fig3.n.CTS", nrow(cts_dat))
record_value("fig3.cohorts.CTS", length(unique(cts_dat$dataset)))
record_value("fig3.n.MARS", nrow(mars))
record_value("fig3.n.SoM", nrow(som_dat))
record_value("fig3.cohorts.SoM", length(unique(som_dat$dataset)))
cat(sprintf("   CTS n=%d (%d cohorts); MARS n=%d; SoM n=%d (%d cohorts)\n", nrow(cts_dat),
            length(unique(cts_dat$dataset)), nrow(mars), nrow(som_dat), length(unique(som_dat$dataset))))
v <- function(cl, g, a) agg$cohort_balanced_contrast[agg$classifier == cl & agg$group == g & agg$axis == a]
for (cl_g in list(c("MARS", "Mars1"), c("MARS", "Mars3"), c("MARS", "Mars4"), c("CTS", "CTS2"),
                  c("CTS", "CTS3"), c("SoM", "SoM-Mid"), c("SoM", "SoM-High"))) {
  for (a in AXIS_CODES) record_value(sprintf("fig3.%s.%s.%s", cl_g[1], cl_g[2], a), v(cl_g[1], cl_g[2], a))
}
# "Increasing SoM tertiles showed higher M-D, metabolic, endothelial and lower L, M-P, I scores"
for (a in c("M", "G", "E")) record_value(paste0("fig3.SoM.increasing.higher.", a),
                                         as.numeric(0 < v("SoM", "SoM-Mid", a) && v("SoM", "SoM-Mid", a) < v("SoM", "SoM-High", a)))
for (a in c("L", "N_MP", "I")) record_value(paste0("fig3.SoM.increasing.lower.", a),
                                            as.numeric(0 > v("SoM", "SoM-Mid", a) && v("SoM", "SoM-Mid", a) > v("SoM", "SoM-High", a)))

# ---- Figure 3 ----------------------------------------------------------------------------------
profile_plot <- function(classifier, groups, title, group_labels = groups) {
  z <- agg[agg$classifier == classifier & agg$group %in% groups, ]
  assert(setequal(unique(z$group), groups), classifier, ": every group, including the reference, must be present")
  z$group <- factor(z$group, levels = groups, labels = group_labels)
  z$axis_lab <- factor(z$axis, levels = rev(AXIS_CODES), labels = AXIS_LABEL[rev(AXIS_CODES)])
  z$direction <- ifelse(abs(z$cohort_balanced_contrast) < 1e-12, "Reference",
                        ifelse(z$cohort_balanced_contrast > 0, "Higher", "Lower"))
  ggplot(z, aes(x = cohort_balanced_contrast, y = axis_lab)) +
    geom_vline(xintercept = 0, color = "#6E6E6E", linewidth = 0.5) +
    geom_segment(aes(x = 0, xend = cohort_balanced_contrast, yend = axis_lab), color = "#BDBDBD", linewidth = 1.0) +
    geom_point(aes(color = direction), size = 4.2) +
    facet_wrap(~ group, nrow = 1) +
    scale_color_manual(values = c(Higher = "#C9253A", Lower = "#2878B5", Reference = "#8B8B8B"), guide = "none") +
    labs(title = title, x = "Axis-z difference", y = NULL) +
    theme_minimal(base_size = 17) +
    theme(plot.title = element_text(face = "bold", size = 20, hjust = 0), plot.title.position = "plot",
          strip.text = element_text(face = "bold", size = 18), axis.text = element_text(size = 16.5),
          axis.title.x = element_text(size = 17), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "#EAEAEA"), plot.margin = ggplot2::margin(6, 10, 10, 6))
}
fig <- profile_plot("CTS", c("CTS1", "CTS2", "CTS3"), "A. CTS") /
  profile_plot("MARS", c("Mars2", "Mars1", "Mars3", "Mars4"), "B. MARS",
               group_labels = c("MARS2", "MARS1", "MARS3", "MARS4")) /
  profile_plot("SoM", c("SoM-Low", "SoM-Mid", "SoM-High"), "C. SoM") +
  plot_layout(heights = c(1, 1.05, 1))

ggsave(
  filename = file.path(FIG_DIR, "Figure3_endotype_profiles.tiff"),  # 路徑跟副檔名依需求調整
  plot = fig,
  width = 15.5,
  height = 17.0,
  units = "in",        # 需確認：原本save_figure用的單位是in還是cm
  dpi = 300,
  compression = "lzw"  # 如果輸出格式是tiff才需要這個參數
)

finish_script()
