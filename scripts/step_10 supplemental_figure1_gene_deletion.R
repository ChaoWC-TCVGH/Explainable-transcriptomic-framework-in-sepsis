.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))
source(file.path(.script_dir, "lib", "scoring.R"))

start_script("step_10 supplemental_figure1_gene_deletion")

fw <- load_framework()
tagged <- grepl("CTS|SoM|HI-?DEF|HiDEF", fw$membership_source, ignore.case = TRUE)
union_genes <- sort(unique(fw$gene_symbol[tagged]))
assert(sum(tagged) == 151L && length(union_genes) == 121L, "classifier union must be 151 memberships / 121 genes")

build_panel <- function(mem, delete = character()) {
  md <- mem[!mem$gene_symbol %in% delete, ]
  gs <- lapply(split(md$gene_symbol, md$module), clean_genes)
  gs <- gs[vapply(gs, length, integer(1)) >= 3]
  ma <- tapply(md$axis, md$module, function(x) x[[1]])
  list(gene_sets = gs, axis_modules = split(names(ma), unname(ma)),
       n_genes = length(unique(md$gene_symbol)), n_modules = length(gs))
}
score_panel <- function(ranked, panel) {
  gs <- lapply(panel$gene_sets, function(g) intersect(clean_genes(g), rownames(ranked)))
  gs <- gs[vapply(gs, length, integer(1)) >= 3]
  mz <- do.call(rbind, lapply(gs, function(g) colMeans(ranked[g, , drop = FALSE], na.rm = TRUE)))
  rownames(mz) <- names(gs)
  mz <- zscore_rows(mz)
  am <- do.call(rbind, lapply(AXIS_CODES, function(a) {
    mm <- intersect(panel$axis_modules[[a]], rownames(mz))
    assert(length(mm) > 0L, "no scorable module left for axis ", a)
    colMeans(mz[mm, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(am) <- AXIS_CODES
  list(axis = t(zscore_rows(am)), scored = rownames(mz))
}

full <- build_panel(fw)
del <- build_panel(fw, union_genes)
record_value("deletion.genes_retained_panel", del$n_genes)
record_value("deletion.modules_structurally_eligible_panel", del$n_modules)
cat(sprintf("   full panel %d genes / %d modules; deleted panel %d genes / %d modules\n",
            full$n_genes, full$n_modules, del$n_genes, del$n_modules))

cohorts <- c("GSE134347", "GSE26440", "GSE57065", "GSE65682", "GSE95233")
ck <- load_axis_scores()
mp_modules <- unique(fw$module[fw$axis == "N_MP"])
corr <- list(); sanity <- list(); mp_scored <- integer()
for (co in cohorts) {
  expr <- readRDS(file.path(EXPR_DIR, paste0(co, ".rds")))
  ranked <- scale(apply(expr, 2, rank, ties.method = "average", na.last = "keep"))
  f <- score_panel(ranked, full)
  g <- score_panel(ranked, del)
  mp_scored[co] <- length(intersect(g$scored, mp_modules))
  ref <- ck[ck$dataset == co, ]
  ref <- ref[match(rownames(f$axis), ref$sample_id), ]
  for (a in AXIS_CODES) {
    corr[[length(corr) + 1L]] <- data.frame(dataset = co, axis = a,
      spearman_rho = cor(f$axis[, a], g$axis[, a], method = "spearman"))
    sanity[[length(sanity) + 1L]] <- data.frame(dataset = co, axis = a,
      max_abs_diff_vs_checkpoint = max(abs(f$axis[, a] - ref[[a]])))
  }
  cat(sprintf("   %-10s %4d samples; M-P modules scorable after deletion: %d\n", co, ncol(expr), mp_scored[co]))
}
corr <- do.call(rbind, corr)
sanity <- do.call(rbind, sanity)
write_csv(corr, file.path(TABLE_DIR, "deletion_full_vs_deleted_spearman.csv"))
write_csv(sanity, file.path(TABLE_DIR, "deletion_full_rescore_vs_checkpoint.csv"))
assert(max(sanity$max_abs_diff_vs_checkpoint) < 1e-6, "full-panel rescore does not reproduce the checkpoint axis scores")

agg <- do.call(rbind, lapply(AXIS_CODES, function(a) {
  v <- corr$spearman_rho[corr$axis == a]
  data.frame(axis = a, mean = mean(v), min = min(v), max = max(v))
}))
write_csv(agg, file.path(TABLE_DIR, "SupplementalFigure1_source_data.csv"))
print(agg, row.names = FALSE, digits = 3)
eight <- agg[agg$axis != "N_MP", ]
mp <- agg[agg$axis == "N_MP", ]
record_value("deletion.cohorts", length(cohorts))
for (i in seq_len(nrow(agg))) {
  record_value(paste0("deletion.mean.", agg$axis[i]), agg$mean[i])
  record_value(paste0("deletion.min.", agg$axis[i]), agg$min[i])
}
record_value("deletion.eight_axes.lowest_mean", min(eight$mean))
record_value("deletion.eight_axes.highest_mean", max(eight$mean))
record_value("deletion.eight_axes.lowest_cohort_value", min(eight$min))
record_value("deletion.MP.mean", mp$mean)
record_value("deletion.MP.min", mp$min)
record_value("deletion.MP.max", mp$max)
record_value("deletion.MP.scorable_modules_min", min(mp_scored))
record_value("deletion.MP.scorable_modules_max", max(mp_scored))

# ---- Supplemental Figure 1 -------------------------------------------------------------------
axis_plot_name <- AXIS_NAME
axis_plot_name[["H"]] <- "heme-coagulation"
agg$label <- factor(paste0(AXIS_DISPLAY[agg$axis], "  ", axis_plot_name[agg$axis]),
                    levels = rev(paste0(AXIS_DISPLAY[AXIS_CODES], "  ", axis_plot_name[AXIS_CODES])))
p <- ggplot(agg, aes(x = mean, y = label)) +
  geom_col(width = 0.65, fill = "#4E79A7") +
  geom_errorbarh(aes(xmin = min, xmax = max), height = 0, colour = "#2B2B2B", linewidth = 0.6) +
  geom_point(aes(x = min), colour = "#2B2B2B", size = 1.9) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 1.0)) +
  labs(x = "Full versus union-deleted Spearman correlation", y = NULL,
       subtitle = "Bar, five-cohort mean; line, minimum to maximum; point, minimum") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9.5), plot.subtitle = element_text(size = 8.5, colour = "#555555"),
        plot.margin = ggplot2::margin(8, 14, 6, 6))
ggsave(
  filename = file.path(FIG_DIR, "SupplementalFigure1.tiff"),  # 路徑跟副檔名依需求調整
  plot = p,
  width = 6.6,
  height = 4.05,
  units = "in",        # 需確認：原本save_figure用的單位是in還是cm
  dpi = 300,
  compression = "lzw"  # 如果輸出格式是tiff才需要這個參數
)
finish_script()
