
.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))
start_script("step_7 figure2_coverage_and_correlation")

fw <- load_framework()
cv <- read_csv(RE_COVERAGE)
assert(nrow(cv) == 612L && length(unique(cv$dataset)) == 12L, "coverage checkpoint must be 12 cohorts x 51 modules")
cv$present_gene_symbols[is.na(cv$present_gene_symbols)] <- ""

# ---- (A) coverage --------------------------------------------------------------------------
axis_cov <- do.call(rbind, lapply(unique(cv$dataset), function(ds) do.call(rbind, lapply(AXIS_CODES, function(a) {
  z <- cv[cv$dataset == ds & cv$axis == a, ]
  axis_genes <- clean_genes(fw$gene_symbol[fw$axis == a])
  present <- clean_genes(unlist(strsplit(paste(z$present_gene_symbols[nzchar(z$present_gene_symbols)],
                                               collapse = ";"), ";", fixed = TRUE)))
  data.frame(dataset = ds, axis = a, axis_genes = length(axis_genes),
             measured_genes = sum(axis_genes %in% present),
             coverage = sum(axis_genes %in% present) / length(axis_genes), stringsAsFactors = FALSE)
}))))
cov_sum <- do.call(rbind, lapply(AXIS_CODES, function(a) {
  z <- axis_cov$coverage[axis_cov$axis == a]
  assert(length(z) == 12L, a, " must have 12 cohorts")
  data.frame(axis = a, label = AXIS_LABEL[[a]], mean_pct = 100 * mean(z), sd_pct = 100 * stats::sd(z))
}))

print(cov_sum, row.names = FALSE, digits = 3)
for (i in seq_len(nrow(cov_sum))) {
  record_value(paste0("coverage.mean_pct.", cov_sum$axis[i]), cov_sum$mean_pct[i])
  record_value(paste0("coverage.sd_pct.", cov_sum$axis[i]), cov_sum$sd_pct[i])
}
ord_mean <- cov_sum$axis[order(cov_sum$mean_pct)]
ord_sd <- cov_sum$axis[order(-cov_sum$sd_pct)]
record_value("coverage.lowest_mean_is_A", as.numeric(ord_mean[1] == "A"))
record_value("coverage.highest_mean_is_L", as.numeric(ord_mean[9] == "L"))
record_value("coverage.two_lowest_are_A_E", as.numeric(setequal(ord_mean[1:2], c("A", "E"))))
record_value("coverage.two_most_variable_are_A_E", as.numeric(setequal(ord_sd[1:2], c("A", "E"))))

# ---- (B) between-axis correlation ------------------------------------------------------------
ax <- load_axis_scores()
per <- lapply(split(ax, ax$dataset), function(d)
  cor(as.matrix(d[, AXIS_CODES]), method = "spearman", use = "pairwise.complete.obs"))
FIG2B_COHORTS <- c("GSE134347", "GSE26440", "GSE57065", "GSE95233")
published <- Reduce(`+`, per[FIG2B_COHORTS]) / length(FIG2B_COHORTS)
median12 <- apply(simplify2array(per), c(1, 2), median)

record_value("fig2b.published_matrix.MD_metabolic", published["M", "G"],
             "mean of within-cohort Spearman over GSE134347, GSE26440, GSE57065, GSE95233 (published Figure 2B)")
record_value("fig2b.published_matrix.n_cohorts", length(FIG2B_COHORTS))
record_value("fig2b.published_matrix.cohorts_are_GSE134347_GSE26440_GSE57065_GSE95233",
             as.numeric(setequal(FIG2B_COHORTS, c("GSE134347", "GSE26440", "GSE57065", "GSE95233"))))
record_value("fig2b.median12.MD_metabolic", median12["M", "G"],
             "median of within-cohort Spearman over all 12 cohorts; for comparison only; not shown in the figure")
cat(sprintf("   M-D vs metabolic: published matrix %.3f | median over 12 cohorts %.3f\n",
            published["M", "G"], median12["M", "G"]))

# ---- Figure 2 ---------------------------------------------------------------------------------
plot_a <- cov_sum
plot_a$axis_label <- factor(plot_a$axis, levels = rev(AXIS_CODES), labels = rev(AXIS_LABEL[AXIS_CODES]))
plot_a$lower <- pmax(0, (plot_a$mean_pct - plot_a$sd_pct) / 100)
plot_a$upper <- pmin(1, (plot_a$mean_pct + plot_a$sd_pct) / 100)
p_cov <- ggplot(plot_a, aes(mean_pct / 100, axis_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#A0A0A0", linewidth = 0.55) +
  geom_col(width = 0.66, fill = "#4C78A8") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.22, linewidth = 0.85, color = "#263746") +
  scale_x_continuous(limits = c(0, 1.06), breaks = seq(0, 1, by = 0.2),
                     labels = function(x) sprintf("%.0f%%", 100 * x), expand = expansion(mult = c(0, 0))) +
  labs(title = "A", x = "Unique-gene coverage (%)", y = NULL) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 26, hjust = 0), plot.title.position = "plot",
        axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.margin = ggplot2::margin(8, 10, 10, 8))
corr_long <- expand.grid(row_axis = AXIS_CODES, col_axis = AXIS_CODES, stringsAsFactors = FALSE)
corr_long$rho <- mapply(function(r, c) published[r, c], corr_long$row_axis, corr_long$col_axis)
corr_long$row_label <- factor(corr_long$row_axis, levels = rev(AXIS_CODES), labels = AXIS_LABEL[rev(AXIS_CODES)])
corr_long$col_label <- factor(corr_long$col_axis, levels = AXIS_CODES, labels = AXIS_DISPLAY[AXIS_CODES])
corr_long$text_color <- ifelse(abs(corr_long$rho) >= 0.60, "white", "#17202A")
p_corr <- ggplot(corr_long, aes(col_label, row_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.65) +
  geom_text(aes(label = sprintf("%.2f", rho), color = text_color), size = 5.80, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho") +
  labs(title = "B", x = NULL, y = NULL) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 26, hjust = 0), plot.title.position = "plot",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18), axis.text.y = element_text(size = 18),
        legend.position = "bottom", panel.grid = element_blank(), plot.margin = ggplot2::margin(8, 8, 10, 10)) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = grid::unit(3.0, "in"),
                               barheight = grid::unit(0.17, "in")))

ggsave(
  filename = file.path(FIG_DIR, "Figure2_axis_coverage_correlation.tiff"),  # 路徑跟副檔名依需求調整
  plot = p_cov / p_corr + plot_layout(heights = c(1, 1.06)),
  width = 13.2,
  height = 10.4,
  units = "in",        # 需確認：原本save_figure用的單位是in還是cm
  dpi = 300,
  compression = "lzw"  # 如果輸出格式是tiff才需要這個參數
)
finish_script()

