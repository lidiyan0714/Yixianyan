# =============================================================================
# Figure 5C / Figure S12A regeneration + Table S10 (detailed) regeneration
# =============================================================================
# Purpose
# -------
# Reproduces, in R, the figure and table that were regenerated for the
# pancreatitis multi-omics manuscript's scRNA-seq DEG panels, using ONE
# consistent, transparent standard applied to BOTH the table and the figure:
#
#     Seurat FindMarkers output, per cluster (C0-C10), ABP vs Control and
#     HTGP vs Control separately.
#     Significance filter:  p_val_adj < 0.01  &  |avg_log2FC| > 0.25
#
# Because Table S10 and Figure 5C/S12A are built from the exact same
# filtered gene sets, their per-cluster totals match by construction.
#
# Data source
# -----------
# scRNAseq/11_diff_analysis_results/{cluster}_{group}_{group}_vs_con_diff_gene.csv
#   group = "ASP" (-> ABP) or "WL" (-> HTGP)
#   columns: gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj
#
# Output
# ------
#   Figure5C_v4.pdf / Figure5C_v4.png   (or FigureS12A_* if TITLE changed below)
#   Table S10 DEGs in each cell type_v2.xlsx  (Summary / ABP / HTGP sheets)
#
# Output text/fonts
# -----------------
# All text is set in Arial, and the PDF is saved via the cairo_pdf device
# so every letter is embedded as a real scalable/selectable outline (never
# a rasterized bitmap) -- gene names can be selected and copied
# individually in a PDF viewer and stay crisp at any zoom. Windows ships
# Arial natively, so running this script as-is on Windows (this script's
# intended platform) renders true Arial with no extra setup. If Arial is
# not found (e.g. on Linux/macOS without it installed), the code below
# falls back to "sans" (the system default sans-serif) so the script still
# runs; install the "arial" font family on that machine to get true Arial.
#
# The figure is sized for A4 landscape (297 x 210 mm = 11.69 x 8.27 in) so
# it drops straight into a manuscript/print layout.
#
# Dependencies: tidyverse, ggplot2, gridExtra, openxlsx, systemfonts (for
# reliable font-availability checking)
#   install.packages(c("tidyverse","gridExtra","openxlsx","systemfonts"))
# =============================================================================

library(tidyverse)
library(gridExtra)
library(grid)
library(openxlsx)

# pick Arial if installed on this machine, else fall back to plain "sans"
FONT_FAMILY <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  if ("Arial" %in% systemfonts::system_fonts()$family) FONT_FAMILY <- "Arial"
} else {
  FONT_FAMILY <- "Arial"  # optimistic default on Windows; falls back silently if absent
}

# ---------------------------------------------------------------------------
# 0. Configuration -- EDIT THESE TWO PATHS FOR YOUR MACHINE
# ---------------------------------------------------------------------------
DATA_DIR <- "./scRNAseq/11_diff_analysis_results"   # folder with the 22 CSVs
OUT_DIR  <- "./Figures and tables"                   # where outputs are written

# Figure title / output basename -- change these two lines to regenerate
# Figure S12A instead of Figure 5C (they use the exact same data/standard).
TITLE        <- "Figure 5C (regenerated)"
OUT_BASENAME <- "Figure5C_v4_R"
# TITLE        <- "Figure S12A (regenerated)"
# OUT_BASENAME <- "FigureS12A_v4_R"

PADJ_TH <- 0.01
LFC_TH  <- 0.25
N_LABELS <- 5    # top |avg_log2FC| genes labeled per cluster per group

CLUSTERS <- paste0("C", 0:10)

CELL_TYPE <- c(
  C0  = "NK cell",
  C1  = "Naive CD4+ T cell",
  C2  = "Inflammatory monocytes (CD14+)",
  C3  = "B cell",
  C4  = "Activated/memory CD4+ T cell",
  C5  = "Monocytes/macrophages",
  C6  = "Plasma cell",
  C7  = "T cell (CD8+ naive/early activated)",
  C8  = "Non-classical monocytes (CD16+)",
  C9  = "Plasmacytoid dendritic cell (pDC)",
  C10 = "Hematopoietic stem/progenitor cell (HSPC)"
)

CLUSTER_COLORS <- c(
  C0 = "#B0B0B0", C1 = "#5FBFC0", C2 = "#2E9E6E", C3 = "#31447A", C4 = "#E8776B",
  C5 = "#7C8FC4", C6 = "#7FCBB0", C7 = "#B7E39A", C8 = "#7A5230", C9 = "#C9B79C", C10 = "#3D6FB4"
)

RED  <- "#F08072"   # ABP, significant
TEAL <- "#5FBFC0"   # HTGP, significant
GRAY <- "#BFBFBF"   # not significant

GROUP_CODE <- c(ABP = "ASP", HTGP = "WL")

# ---------------------------------------------------------------------------
# 1. Load all 22 CSVs, mark significance
# ---------------------------------------------------------------------------
load_one <- function(cluster, group) {
  code <- GROUP_CODE[[group]]
  path <- file.path(DATA_DIR, sprintf("%s_%s_%s_vs_con_diff_gene.csv", cluster, code, code))
  df <- suppressMessages(read_csv(path, show_col_types = FALSE))
  names(df)[1] <- "Gene_symbol"
  df %>%
    mutate(
      Cluster = cluster,
      Group = group,
      sig = (p_val_adj < PADJ_TH) & (abs(avg_log2FC) > LFC_TH),
      Direction = ifelse(avg_log2FC > 0, "Up", "Down")
    )
}

all_data <- map_dfr(CLUSTERS, function(cl) {
  bind_rows(load_one(cl, "ABP"), load_one(cl, "HTGP"))
})

# ---------------------------------------------------------------------------
# 2. Per-cluster stats: totals, up/down, same-direction ABP-HTGP overlap
# ---------------------------------------------------------------------------
stats <- map_dfr(CLUSTERS, function(cl) {
  abp  <- all_data %>% filter(Cluster == cl, Group == "ABP", sig)
  htgp <- all_data %>% filter(Cluster == cl, Group == "HTGP", sig)
  abp_up    <- abp   %>% filter(Direction == "Up")   %>% pull(Gene_symbol)
  abp_down  <- abp   %>% filter(Direction == "Down") %>% pull(Gene_symbol)
  htgp_up   <- htgp  %>% filter(Direction == "Up")   %>% pull(Gene_symbol)
  htgp_down <- htgp  %>% filter(Direction == "Down") %>% pull(Gene_symbol)
  tibble(
    Cluster = cl,
    Cell_type = CELL_TYPE[[cl]],
    ABP_Up = length(abp_up), ABP_Down = length(abp_down), ABP_Total = nrow(abp),
    ABP_Tested = all_data %>% filter(Cluster == cl, Group == "ABP") %>% nrow(),
    HTGP_Up = length(htgp_up), HTGP_Down = length(htgp_down), HTGP_Total = nrow(htgp),
    HTGP_Tested = all_data %>% filter(Cluster == cl, Group == "HTGP") %>% nrow(),
    Overlap_Up = length(intersect(abp_up, htgp_up)),
    Overlap_Down = length(intersect(abp_down, htgp_down))
  )
})
print(stats)

# ---------------------------------------------------------------------------
# 3. Table S10 v2 (detailed) -- Summary / ABP / HTGP sheets
# ---------------------------------------------------------------------------
detail <- all_data %>%
  filter(sig) %>%
  mutate(Comparison = paste0(Group, " vs Control"),
         Cell_type = CELL_TYPE[Cluster]) %>%
  select(Cluster, Cell_type, Comparison, Gene_symbol, avg_log2FC, Direction,
         p_val, p_val_adj, pct.1, pct.2) %>%
  arrange(Cluster, p_val_adj)

wb <- createWorkbook()
addWorksheet(wb, "Summary"); writeData(wb, "Summary", stats)
addWorksheet(wb, "ABP");     writeData(wb, "ABP", detail %>% filter(Comparison == "ABP vs Control"))
addWorksheet(wb, "HTGP");    writeData(wb, "HTGP", detail %>% filter(Comparison == "HTGP vs Control"))
saveWorkbook(wb, file.path(OUT_DIR, "Table S10 DEGs in each cell type_v2_R.xlsx"), overwrite = TRUE)

# ---------------------------------------------------------------------------
# 4. Build the scatter panel
# ---------------------------------------------------------------------------
set.seed(42)
col_width <- 1.0
gap <- 0.12
x_centers <- setNames((seq_along(CLUSTERS) - 1) * (col_width + gap), CLUSTERS)

all_data <- all_data %>%
  mutate(
    xc = x_centers[Cluster],
    jitter = (runif(n()) - 0.5) * (col_width * 0.86),
    dot_x = xc + jitter   # the ACTUAL x-position each gene is plotted at
  )

plot_df <- all_data %>%
  mutate(
    x = dot_x,
    color_group = case_when(
      sig & Group == "ABP"  ~ "ABP (significant)",
      sig & Group == "HTGP" ~ "HTGP (significant)",
      TRUE ~ "Not significant"
    )
  )
plot_df$color_group <- factor(plot_df$color_group,
                               levels = c("Not significant", "ABP (significant)", "HTGP (significant)"))
plot_df <- plot_df %>% arrange(color_group)  # draw gray first, colored on top

# top |logFC| labeled genes PER CLUSTER, pooled across BOTH groups (not 5
# per group) -- keep dot_x (their true plotted position) so a SHORT leader
# line can be drawn from the label to the actual dot.
#
# NOTE: we deliberately do NOT use ggrepel::geom_text_repel() here. ggrepel's
# collision-avoidance runs across the WHOLE plot canvas in one pass, which
# pushed labels far outside their own cluster's column into neighboring
# clusters (long diagonal leader lines crossing the entire figure). Instead
# we mirror the already-verified Python/matplotlib design: each label is
# nudged just a SHORT distance away from its own dot (not all the way to a
# fixed far column position), then clipped so it still can't cross into a
# neighboring cluster's column. This keeps leader lines short and keeps
# labels/lines confined to their own column, including for the first (C0)
# and last (C10) columns which have no neighboring column to absorb
# overflow.
label_offset <- col_width * 0.14
label_margin <- col_width * 0.06
label_df <- all_data %>%
  filter(sig) %>%
  group_by(Cluster) %>%
  slice_max(order_by = abs(avg_log2FC), n = N_LABELS) %>%
  ungroup() %>%
  mutate(
    raw_lx = dot_x + ifelse(Group == "ABP", label_offset, -label_offset),
    lx = pmin(pmax(raw_lx, xc - col_width/2 + label_margin), xc + col_width/2 - label_margin),
    ly = avg_log2FC
  )

# background column shading (alternate gray/white)
shade_df <- tibble(Cluster = CLUSTERS, xc = x_centers[CLUSTERS]) %>%
  mutate(fill = ifelse(row_number() %% 2 == 1, "#F7F7F7", "#FFFFFF"))

p <- ggplot() +
  geom_rect(data = shade_df, aes(xmin = xc - col_width/2, xmax = xc + col_width/2,
                                  ymin = -13, ymax = 13), fill = shade_df$fill, color = NA) +
  geom_point(data = plot_df, aes(x = x, y = avg_log2FC, color = color_group),
             size = 0.35, alpha = 0.5, show.legend = TRUE) +
  scale_color_manual(values = c("Not significant" = GRAY,
                                 "ABP (significant)" = RED,
                                 "HTGP (significant)" = TEAL), name = NULL) +
  geom_segment(data = label_df, aes(x = dot_x, y = avg_log2FC, xend = lx, yend = ly),
               color = "#999999", linewidth = 0.25, alpha = 0.8) +
  geom_text(data = label_df, aes(x = lx, y = ly, label = Gene_symbol),
            size = 1.9, fontface = "italic", family = FONT_FAMILY, hjust = 0.5, vjust = 0.5) +
  # cluster ID strip
  geom_rect(data = shade_df, aes(xmin = xc - col_width/2 + 0.02, xmax = xc + col_width/2 - 0.02,
                                  ymin = -0.55, ymax = 0.55), fill = CLUSTER_COLORS[shade_df$Cluster]) +
  geom_text(data = shade_df, aes(x = xc, y = 0, label = Cluster),
            color = "white", fontface = "bold", family = FONT_FAMILY, size = 2.6) +
  scale_x_continuous(breaks = NULL) +
  coord_cartesian(ylim = c(-13, 13), xlim = c(-0.7, max(x_centers) + 0.7), clip = "off") +
  labs(
    title = paste0(TITLE, ". Differentially expressed genes across cell clusters in ABP and HTGP relative to controls."),
    subtitle = paste0(
      "Standard (shared with Table S10 v2): Seurat FindMarkers, p_val_adj < 0.01 and |avg_log2FC| > 0.25.  ",
      "Gray = not significant; colored = significant.\n",
      "Top 5 genes by |avg_log2FC| are labeled per cluster (pooled across ABP and HTGP); full DEG counts ",
      "and up/down/overlap statistics are in Table S10 v2 (Summary sheet)."
    ),
    x = NULL,
    y = "avg_log2FC  (ABP vs Control, above 0   |   HTGP vs Control, below 0)"
  ) +
  theme_minimal(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = FONT_FAMILY),
    plot.title = element_text(face = "plain", size = 8.5, hjust = 0, family = FONT_FAMILY),
    plot.subtitle = element_text(size = 6.8, hjust = 0, family = FONT_FAMILY),
    axis.text = element_text(size = 6, family = FONT_FAMILY),
    axis.title = element_text(size = 7, family = FONT_FAMILY),
    legend.text = element_text(size = 6.5, family = FONT_FAMILY),
    legend.position = "right"
  )

# ---------------------------------------------------------------------------
# 5. Save the plot (no separate reference table -- that information is
#    already in Table S10 v2's Summary sheet, saved above)
#
# A4 landscape (297 x 210 mm = 11.69 x 8.27 in). cairo_pdf embeds real
# scalable/selectable font outlines (not a rasterized bitmap font), so
# every gene name can be selected and copied individually in a PDF viewer.
# ---------------------------------------------------------------------------
ggsave(file.path(OUT_DIR, paste0(OUT_BASENAME, ".pdf")), p,
       width = 11.69, height = 8.27, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, paste0(OUT_BASENAME, ".png")), p,
       width = 11.69, height = 8.27, units = "in", dpi = 300)

cat("Saved:", file.path(OUT_DIR, paste0(OUT_BASENAME, ".pdf")), "| font used:", FONT_FAMILY, "\n")
cat("Saved:", file.path(OUT_DIR, "Table S10 DEGs in each cell type_v2_R.xlsx"), "\n")
