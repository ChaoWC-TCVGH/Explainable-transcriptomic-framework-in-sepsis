.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))

start_script("step_9 figure4_mars_cts_correspondence")

ax <- load_axis_scores()
meta <- read_csv(RE_GEO_META)
cts <- load_cts_labels("GSE65682")

mars <- meta[meta$dataset == "GSE65682" & meta$endotype_class %in% paste0("Mars", 1:4),
             c("geo_accession", "endotype_class")]
names(mars) <- c("sample_id", "MARS")
d <- merge(mars, setNames(cts[cts$CTS %in% paste0("CTS", 1:3), c("Sample", "CTS")], c("sample_id", "CTS")),
           by = "sample_id")
assert(nrow(d) == 479L, "MARS-CTS set must contain 479 GSE65682 patients")

tab <- table(MARS = d$MARS, CTS = d$CTS)
write_csv(as.data.frame.matrix(tab) |> cbind(MARS = rownames(tab), x = _), file.path(TABLE_DIR, "Figure4A_MARS_CTS_crosstab.csv"))
print(tab)
record_value("fig4.n_total", sum(tab))
record_value("fig4.n_displayed_MARS1to3", sum(tab[c("Mars1", "Mars2", "Mars3"), ]))
for (m in rownames(tab)) {
  record_value(paste0("fig4.row_total.", m), sum(tab[m, ]))
  for (k in colnames(tab)) record_value(sprintf("fig4.count.%s.%s", m, k), tab[m, k])
}
record_value("fig4.pct.Mars1.CTS2", 100 * tab["Mars1", "CTS2"] / sum(tab["Mars1", ]))
record_value("fig4.pct.Mars3.CTS3", 100 * tab["Mars3", "CTS3"] / sum(tab["Mars3", ]))
record_value("fig4.pct.Mars1.CTS1", 100 * tab["Mars1", "CTS1"] / sum(tab["Mars1", ]))

chi <- suppressWarnings(chisq.test(tab))
n <- sum(tab)
cramers_v <- sqrt(as.numeric(chi$statistic) / (n * (min(dim(tab)) - 1)))
set.seed(1)
B <- 5000
perm_chi <- replicate(B, suppressWarnings(chisq.test(table(d$MARS, sample(d$CTS)))$statistic))
perm_p <- (1 + sum(perm_chi >= as.numeric(chi$statistic))) / (B + 1)
record_value("fig4.cramers_v", cramers_v)
record_value("fig4.permutation_p", perm_p)
cat(sprintf("   Cramer's V = %.3f; permutation P = %.2g (B = %d)\n", cramers_v, perm_p, B))

# ---- discordant subgroup contrasts ------------------------------------------------------------
dd <- merge(d, ax[ax$dataset == "GSE65682", c("sample_id", AXIS_CODES)], by = "sample_id")
contrasts <- list(c("Mars1", "CTS1", "CTS2"), c("Mars3", "CTS1", "CTS3"))
drift <- do.call(rbind, lapply(contrasts, function(cc) {
  z <- dd[dd$MARS == cc[1], ]
  do.call(rbind, lapply(AXIS_CODES, function(a) data.frame(
    MARS = cc[1], discordant = cc[2], concordant = cc[3], axis = a,
    n_discordant = sum(z$CTS == cc[2]), n_concordant = sum(z$CTS == cc[3]),
    diff = mean(z[[a]][z$CTS == cc[2]]) - mean(z[[a]][z$CTS == cc[3]]), stringsAsFactors = FALSE)))
}))

for (i in seq_len(nrow(drift))) record_value(sprintf("fig4b.%s.%s", drift$MARS[i], drift$axis[i]), drift$diff[i])
for (m in c("Mars1", "Mars3")) {
  z <- drift[drift$MARS == m, ][1, ]
  record_value(paste0("fig4b.n_discordant.", m), z$n_discordant)
  record_value(paste0("fig4b.n_concordant.", m), z$n_concordant)
}

# ---- Figure 4 ------------------------------------------------------------------------------------
t3 <- as.data.frame(tab, stringsAsFactors = FALSE)
t3 <- t3[t3$MARS != "Mars4", ]
t3$row_prop <- t3$Freq / ave(t3$Freq, t3$MARS, FUN = sum)
t3$dominant <- ave(t3$Freq, t3$MARS, FUN = function(z) z == max(z)) == 1
t3$label <- ifelse(t3$Freq == 0, "", sprintf("%d\n(%.0f%%)", t3$Freq, 100 * t3$row_prop))
t3$MARS <- factor(t3$MARS, levels = rev(paste0("Mars", 1:3)), labels = rev(paste0("MARS", 1:3)))
t3$CTS <- factor(t3$CTS, levels = paste0("CTS", 1:3))
p_grid <- ggplot(t3, aes(CTS, MARS)) +
  geom_tile(fill = "white", color = "#757575", linewidth = 0.9) +
  geom_tile(data = t3[t3$dominant & t3$Freq > 0, ], fill = "#E5EFF9", color = "#4F4F4F", linewidth = 1.0) +
  geom_text(aes(label = label, fontface = ifelse(dominant & Freq > 0, "bold", "plain")),
            size = 6.2, lineheight = 0.92, color = "#17202A") +
  labs(x = "Package-derived CTS group", y = "Deposited MARS group", title = "A") +
  coord_equal() + theme_minimal(base_size = 16) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 22, hjust = 0),
        plot.title.position = "plot", axis.text = element_text(size = 16), axis.title = element_text(size = 16))
drift$axis_label <- factor(drift$axis, levels = rev(AXIS_CODES),
                           labels = paste0(AXIS_DISPLAY[rev(AXIS_CODES)], "  ", AXIS_NAME[rev(AXIS_CODES)]))
drift$facet <- sprintf("%s: CTS1 (n=%d) vs %s (n=%d)", sub("Mars", "MARS", drift$MARS),
                       drift$n_discordant, drift$concordant, drift$n_concordant)
drift$direction <- ifelse(drift$diff >= 0, "Higher", "Lower")
p_drift <- ggplot(drift, aes(diff, axis_label, color = direction)) +
  geom_segment(aes(x = 0, xend = diff, yend = axis_label), linewidth = 1.0, color = "#B8B8B8") +
  geom_point(size = 3.7) + geom_vline(xintercept = 0, color = "#555555", linewidth = 0.6) +
  facet_wrap(~facet, ncol = 1) +
  scale_color_manual(values = c(Higher = "#B2182B", Lower = "#2166AC"), name = NULL) +
  labs(x = "Mean difference (z-score)", y = NULL, title = "B") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14), strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 22, hjust = 0), plot.title.position = "plot")

ggsave(
  filename = file.path(FIG_DIR, "Figure4_MARS_CTS.tiff"),  # 路徑跟副檔名依需求調整
  plot = p_grid + p_drift + plot_layout(widths = c(1.3, 1.65)),
  width = 12.6,
  height = 5.6,
  units = "in",        # 需確認：原本save_figure用的單位是in還是cm
  dpi = 300,
  compression = "lzw"  # 如果輸出格式是tiff才需要這個參數
)
finish_script()
