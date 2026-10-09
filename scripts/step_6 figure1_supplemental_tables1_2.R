
.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))
start_script("step_6 figure1_supplemental_tables1_2.R")

fw <- load_framework()

# ---- Supplemental Table 1 ------------------------------------------------------------
st1 <- do.call(rbind, lapply(AXIS_CODES, function(a) {
  z <- fw[fw$axis == a, ]
  data.frame(code = AXIS_DISPLAY[[a]], axis = AXIS_NAME[[a]], modules = length(unique(z$module)),
             memberships = nrow(z), unique_genes = length(unique(z$gene_symbol)), stringsAsFactors = FALSE)
}))
st1 <- rbind(st1, data.frame(code = "", axis = "Total (unique)", modules = length(unique(fw$module)),
                             memberships = nrow(fw), unique_genes = length(unique(fw$gene_symbol))))
write_csv(st1, file.path(TABLE_DIR, "SupplementalTable1_axis_structure.csv"))
print(st1, row.names = FALSE)
for (i in seq_along(AXIS_CODES)) {
  a <- AXIS_CODES[i]
  record_value(paste0("st1.modules.", a), st1$modules[i])
  record_value(paste0("st1.memberships.", a), st1$memberships[i])
  record_value(paste0("st1.unique_genes.", a), st1$unique_genes[i])
}
record_value("framework.modules", length(unique(fw$module)))
record_value("framework.memberships", nrow(fw))
record_value("framework.unique_genes", length(unique(fw$gene_symbol)))
record_value("framework.axes", length(unique(fw$axis)))

# ---- Supplemental Table 2 (module display names as printed) -------------------------------
display_module <- function(x) {
  x <- gsub("Monocyte_Protective_Reparative", "Monocyte_Reparative", x, fixed = TRUE)
  gsub("Neutrophil_Protective_Regulatory", "Neutrophil_Regulatory", x, fixed = TRUE)
}
st2 <- do.call(rbind, lapply(split(fw, fw$module), function(z) {
  data.frame(axis = z$axis[1], biological_axis = AXIS_NAME[[z$axis[1]]],
             module_code = z$module[1], module = display_module(z$module_name[1]),
             n_genes = length(unique(z$gene_symbol)),
             genes = paste(sort(unique(z$gene_symbol)), collapse = ", "), stringsAsFactors = FALSE)
}))
st2 <- st2[order(match(st2$axis, AXIS_CODES), -st2$n_genes), ]
write_csv(st2, file.path(TABLE_DIR, "SupplementalTable2_module_gene_composition.csv"))

# ---- provenance of classifier-derived genes ----------------------------------------------
##敏感性測試 當刪除掉原始論文內會有的axis 
tag_cts <- grepl("CTS", fw$membership_source, ignore.case = TRUE)
tag_som <- grepl("SoM", fw$membership_source, ignore.case = TRUE)
tag_hid <- grepl("HI-?DEF|HiDEF", fw$membership_source, ignore.case = TRUE)
tagged <- tag_cts | tag_som | tag_hid
union_genes <- sort(unique(fw$gene_symbol[tagged]))
record_value("provenance.tagged_memberships", sum(tagged))
record_value("provenance.union_genes", length(union_genes))
record_value("provenance.other_genes", length(unique(fw$gene_symbol)) - length(union_genes))

# the three published source lists as locked into the framework
src_cts <- unique(fw$gene_symbol[fw$membership_source == "CTS_reference_docx_classifier_gene"])
src_som <- unique(fw$gene_symbol[fw$membership_source == "SoM_reference_docx_42_gene_signature"])
src_hid <- unique(fw$gene_symbol[fw$module %in% c("HD_MD", "Np1", "Np2", "HD_LP")])
record_value("provenance.source_list.CTS", length(src_cts))
record_value("provenance.source_list.SoM", length(src_som))
record_value("provenance.source_list.HIDEF", length(src_hid))
record_value("provenance.source_lists_union", length(unique(c(src_cts, src_som, src_hid))),
             "all genes of the four HI-DEF modules plus the CTS and SoM genes, including genes without a classifier tag; the figure legend reports the 121 tagged genes")
named_modules <- unique(fw$module_name[grepl("\\((HI-DEF-like|SoM-M[0-9]|CTS2-aligned|SRS-aligned)\\)", fw$module_name)])
record_value("provenance.modules_with_classifier_label", length(named_modules))
cat(sprintf("   classifier union: %d genes from %d tagged memberships; %d modules carry a provenance label\n",
            length(union_genes), sum(tagged), length(named_modules)))

# ---- 比對關鍵字 並刪除------------------------------
kept <- fw[!fw$gene_symbol %in% union_genes, ]
kept_sizes <- tapply(kept$gene_symbol, kept$module, function(g) length(unique(g)))
record_value("deletion.genes_retained", length(unique(kept$gene_symbol)))
record_value("deletion.modules_structurally_eligible", sum(kept_sizes >= 3))
attr_rows <- do.call(rbind, lapply(AXIS_CODES, function(a) {
  genes <- unique(fw$gene_symbol[fw$axis == a])
  lost <- intersect(genes, union_genes)
  on <- function(tag) unique(fw$gene_symbol[tagged & tag])
  hid_only <- setdiff(intersect(lost, on(tag_hid)), union(on(tag_cts), on(tag_som)))
  data.frame(axis = a, display = AXIS_DISPLAY[[a]], axis_genes = length(genes), genes_lost = length(lost),
             pct_lost = 100 * length(lost) / length(genes),
             lost_hidef_tagged = length(intersect(lost, on(tag_hid))),
             lost_cts_tagged = length(intersect(lost, on(tag_cts))),
             lost_som_tagged = length(intersect(lost, on(tag_som))),
             lost_hidef_only = length(hid_only),
             modules = length(unique(fw$module[fw$axis == a])),
             modules_structurally_eligible = sum(kept_sizes[intersect(names(kept_sizes), unique(fw$module[fw$axis == a]))] >= 3),
             stringsAsFactors = FALSE)
}))
write_csv(attr_rows, file.path(TABLE_DIR, "deletion_gene_attribution_by_axis.csv"))
print(attr_rows, row.names = FALSE, digits = 3)
mp <- attr_rows[attr_rows$axis == "N_MP", ]
record_value("deletion.MP.genes_lost", mp$genes_lost)
record_value("deletion.MP.axis_genes", mp$axis_genes)
record_value("deletion.MP.lost_hidef_tagged", mp$lost_hidef_tagged)
record_value("deletion.MP.lost_hidef_only", mp$lost_hidef_only)
record_value("deletion.MP.lost_cts_tagged", mp$lost_cts_tagged)
record_value("deletion.MP.hidef_share_of_lost_pct", 100 * mp$lost_hidef_tagged / mp$genes_lost)
record_value("deletion.MP.hidef_lost_pct_of_axis_genes", 100 * mp$lost_hidef_tagged / mp$axis_genes)
record_value("deletion.MP.lost_pct_of_axis_genes", mp$pct_lost)
record_value("deletion.MP.modules", mp$modules)

# ---- mean module-to-own-axis correlation (pooled 3,083 profiles) --------------------------
mo <- load_module_scores()
ax <- load_axis_scores()
key <- paste(mo$dataset, mo$sample_id)
ax <- ax[match(key, paste(ax$dataset, ax$sample_id)), ]
module_axis <- tapply(fw$axis, fw$module, function(x) x[[1]])
coh <- do.call(rbind, lapply(AXIS_CODES, function(a) {
  ms <- intersect(names(module_axis)[module_axis == a], names(mo))
  r <- vapply(ms, function(m) suppressWarnings(cor(mo[[m]], ax[[a]], use = "pairwise")), numeric(1))
  cm <- suppressWarnings(cor(as.matrix(mo[, ms, drop = FALSE]), use = "pairwise"))
  data.frame(axis = a, display = AXIS_DISPLAY[[a]], mean_module_to_axis_r = mean(r, na.rm = TRUE),
             mean_inter_module_r = mean(cm[upper.tri(cm)], na.rm = TRUE), stringsAsFactors = FALSE)
}))
write_csv(coh, file.path(TABLE_DIR, "axis_coherence.csv"))
print(coh, row.names = FALSE, digits = 3)
for (i in seq_len(nrow(coh))) record_value(paste0("coherence.module_to_axis.", coh$axis[i]), coh$mean_module_to_axis_r[i])
record_value("coherence.module_to_axis.min_axis_is_DT", as.numeric(coh$axis[which.min(coh$mean_module_to_axis_r)] == "T"))
record_value("coherence.module_to_axis.max_axis_is_L", as.numeric(coh$axis[which.max(coh$mean_module_to_axis_r)] == "L"))

# ---- Figure 1 ----------------------------------------------------------------------------------
mod <- aggregate(gene_symbol ~ axis + module + module_name, data = fw, FUN = function(g) length(unique(g)))
names(mod)[names(mod) == "gene_symbol"] <- "n_genes"
assert(nrow(mod) == 51L, "Figure 1 needs 51 modules")
mod$axis <- factor(mod$axis, levels = AXIS_CODES)
mod <- mod[order(mod$axis, -mod$n_genes), ]
mod$idx <- ave(seq_len(nrow(mod)), mod$axis, FUN = seq_along) - 1
mod$subcol <- mod$idx %% 2
mod$subrow <- mod$idx %/% 2
rows_a <- pmax(tapply(mod$subrow, mod$axis, max) + 1, 2)
gap <- 0.55; top_of <- setNames(numeric(), character()); cum <- 0
for (a in AXIS_CODES) { top_of[a] <- cum; cum <- cum + rows_a[[a]] + gap }
block_top <- top_of + 1; block_bot <- top_of + rows_a
grid_top <- max(block_bot)
mod$y <- grid_top - (top_of[as.character(mod$axis)] + mod$subrow + 1) + 1
ytop <- grid_top
mod$label <- sprintf("%s (%d)", display_module(mod$module_name), mod$n_genes)
must_have <- c("Scavenger_Resolution_Support (3)", "Monocyte_Reparative (HI-DEF-like) (34)",
               "Lymphoid_Protective (HI-DEF-like) (29)", "BNK_Protect (19)", "TypeII_IFNg (15)",
               "TNF_IL1_Storm (19)", "Glycocalyx_EndoActiv (24)", "AntagonisticPleiotropy (18)",
               "DDR_Tolerance (16)", "SASP_Decompensation (16)", "AntiInflam_Resol (17)")
assert(all(must_have %in% mod$label), "Figure 1 labels drifted: ", paste(setdiff(must_have, mod$label), collapse = " | "))

axcol <- setNames(c("#E02930", "#0070FD", "#1D8937", "#714BB2", "#EE655A", "#12A8F4",
                    "#BB7F15", "#088377", "#898889"), AXIS_CODES)
mod$fill <- axcol[as.character(mod$axis)]
axis_long <- c(M = "M-D  Myeloid-detrimental", N_MP = "M-P  Myeloid-protective", L = "L  Lymphoid",
               I = "I  Interferon", A = "A  Cytokine/chemokine", E = "E  Endothelial",
               H = "H  Heme/coagulation", G = "M  Metabolic", T = "D-T  Damage-tolerance")
axinfo <- do.call(rbind, lapply(AXIS_CODES, function(a) {
  ymin_a <- grid_top - block_bot[[a]] + 1 - 0.5
  ymax_a <- grid_top - block_top[[a]] + 1 + 0.5
  data.frame(axis = a, ymin = ymin_a, ymax = ymax_a, yc = (ymin_a + ymax_a) / 2,
             fill = axcol[[a]], label = axis_long[[a]], stringsAsFactors = FALSE)
}))
inp_x <- 1.90; inp_w <- 3.60; pw_x <- 5.70; pw_w <- 2.70
ax_x <- 8.75; ax_hw <- 1.50; tile_hw <- 0.14; modL_x <- 10.80; tile_vh <- 0.50
char_w <- 0.0890
modR_x <- modL_x + tile_hw + 0.10 + max(nchar(mod$label[mod$subcol == 0])) * char_w + 0.35 + tile_hw
right_edge <- modR_x + tile_hw + 0.10 + max(nchar(mod$label[mod$subcol == 1])) * char_w + 0.25
mod$tx <- ifelse(mod$subcol == 0, modL_x, modR_x)
ymid <- mean(range(mod$y))
inputs <- data.frame(
  label = c("CTS source list\n(18 genes)", "HI-DEF source list\n(104 genes)",
            "SUBSPACE-SoM source list\n(42 genes)",
            "Literature-curated\nsepsis mechanisms\ne.g. disease tolerance, heme\ntoxicity, metabolic effect",
            "Reactome / MSigDB\npathway-library genes\ne.g. IL6/JAK/STAT3,\ncomplement, glycolysis"),
  fill = c("#FD7D00", "#FD7D00", "#FD7D00", "#43CD4C", "#0070FC"),
  text_col = c("grey8", "grey8", "grey8", "grey8", "white"), stringsAsFactors = FALSE)
grid_y_top <- ytop + 0.5
grid_span <- grid_y_top - min(axinfo$ymin)
slot_h <- 0.80 * grid_span / nrow(inputs)
inputs$y <- grid_y_top - slot_h * (seq_len(nrow(inputs)) - 1) - slot_h / 2
inputs$h <- slot_h * 0.86
legend_df <- data.frame(label = c("Classifier-derived (hypothesis-anchored)", "Literature-curated mechanism",
                                  "Reactome / MSigDB pathway library"),
                        fill = c("#FD7D00", "#43CD4C", "#0070FC"), frac = c(0.855, 0.925, 0.995))
legend_df$y <- grid_y_top - legend_df$frac * grid_span
leg_x0 <- inp_x - inp_w / 2; leg_x1 <- leg_x0 + 0.44; leg_vh <- 0.55

p1 <- ggplot() +
  geom_rect(data = inputs, aes(xmin = inp_x - inp_w / 2, xmax = inp_x + inp_w / 2, ymin = y - h / 2, ymax = y + h / 2),
            fill = inputs$fill, color = "grey35", linewidth = 0.35) +
  geom_text(data = inputs, aes(x = inp_x, y = y, label = label), size = 5.8, lineheight = 0.98,
            color = inputs$text_col, fontface = "bold") +
  geom_segment(data = inputs, aes(x = inp_x + inp_w / 2, y = y, xend = pw_x - pw_w / 2, yend = ymid),
               arrow = arrow(length = grid::unit(0.11, "cm"), type = "closed"), color = "grey50", linewidth = 0.28) +
  geom_rect(data = legend_df, aes(xmin = leg_x0, xmax = leg_x1, ymin = y - leg_vh, ymax = y + leg_vh),
            fill = legend_df$fill, color = "grey35", linewidth = 0.3) +
  geom_text(data = legend_df, aes(x = leg_x1 + 0.20, y = y, label = label), hjust = 0, size = 5.8,
            color = "grey12", fontface = "bold") +
  annotate("rect", xmin = pw_x - pw_w / 2, xmax = pw_x + pw_w / 2, ymin = ymid - 2.70, ymax = ymid + 2.70,
           fill = "#D1E5F0", color = "grey40", linewidth = 0.4) +
  annotate("text", x = pw_x, y = ymid, size = 4.80, fontface = "bold", lineheight = 1.02,
           label = "Knowledge-based\nmodule construction\n(by biological mechanism\nor pathway, NOT by\nclassifier label)") +
  geom_segment(data = axinfo, aes(x = pw_x + pw_w / 2, y = ymid, xend = ax_x - ax_hw, yend = yc),
               arrow = arrow(length = grid::unit(0.09, "cm"), type = "closed"), color = "grey60", linewidth = 0.24) +
  geom_rect(data = axinfo, aes(xmin = ax_x - ax_hw, xmax = ax_x + ax_hw, ymin = ymin, ymax = ymax),
            fill = axinfo$fill, color = "grey20", linewidth = 0.35) +
  geom_text(data = axinfo, aes(x = ax_x - ax_hw + 0.10, y = yc, label = label), color = "white",
            fontface = "bold", size = 5.26, hjust = 0) +
  geom_rect(data = mod, aes(xmin = tx - tile_hw, xmax = tx + tile_hw, ymin = y - tile_vh, ymax = y + tile_vh),
            fill = mod$fill, color = "white", linewidth = 0.2) +
  geom_text(data = mod, aes(x = tx + tile_hw + 0.10, y = y, label = label), hjust = 0, size = 4.55,
            color = "grey12", fontface = "bold") +
  annotate("text", x = inp_x, y = ytop + 2.2, label = "Inputs\n(provenance)", fontface = "bold", size = 6.4, lineheight = 0.95) +
  annotate("text", x = pw_x, y = ytop + 2.2, label = "Biological\ngrouping", fontface = "bold", size = 6.4, lineheight = 0.95) +
  annotate("text", x = ax_x, y = ytop + 2.2, label = "9 biological axes", fontface = "bold", size = 6.4) +
  annotate("text", x = (modL_x + right_edge) / 2, y = ytop + 2.2, label = "51 biological modules (gene count)",
           fontface = "bold", size = 6.4) +
  scale_x_continuous(limits = c(0, right_edge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(legend_df$y) - 1.2, ytop + 3.9)) +
  theme_void(base_size = 13) +
  theme(legend.position = "none", plot.background = element_rect(fill = "white", color = NA),
        plot.margin = ggplot2::margin(6, 8, 6, 10))
fig_w <- right_edge * 1.02
ggsave(
  filename = file.path(FIG_DIR,"Figure1_framework.tiff"),  # FIGURE_DIR換成save_figure內部用的實際輸出資料夾
  plot = p1,
  width = fig_w,
  height = fig_w / 1.777,
  units = "in",        # 依你原本fig_w定義的單位調整（in/cm/px擇一）
  dpi = 300,
  compression = "lzw"
)
write_csv(mod[, c("axis", "module", "module_name", "n_genes", "label")], file.path(TABLE_DIR,"Figure1_source_data.csv"))
finish_script()
