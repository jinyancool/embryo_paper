pkgs <- c(
  "fs",
  "futile.logger",
  "configr",
  "stringr",
  "ggpubr",
  "ggthemes",
  "jhtools",
  "glue",
  "ggsci",
  "patchwork",
  "tidyverse",
  "dplyr",
  "Seurat",
  "paletteer",
  "cowplot",
  "ComplexHeatmap",
  "circlize"
)
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/fig2")
workdir |> fs::dir_create() |> setwd()

yaml_fn <- "~/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn = yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn = yaml_fn, "stage")[c("E9.5", "E11.5", "E13.5")]

my_theme1 <- theme_classic(base_size = 8) +
  theme(
    legend.key.size = unit(3, "mm"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = .5)
  )

## fig2a-c: specific m/z expression in the selected tissue of mouse in each stage -----
mtb_info <- read_rds("~/ref/kegg/mouse/mmu_mtb_pth_cpd.rds")

xlsx_fn1 <- glue(
  "/cluster/home/jhuang/projects/collabrators/data/wangwenjie/mouse/metabolism/mouse_adjusted/DZLM2023110146_DZLM2024030584-b1-张进-王文杰-空间代谢组-项目报告/2.定性结果/Qualitative.xlsx"
)
neg_anot1 <- readxl::read_excel(xlsx_fn1, sheet = "neg-all") |>
  mutate(mz_id = paste0("neg-", mz)) |> dplyr::filter(KEGG %in% mtb_info$cpd_id, !is.na(Metabolites))
pos_anot1 <- readxl::read_excel(xlsx_fn1, sheet = "pos-all") |>
  mutate(mz_id = paste0("pos-", mz)) |> dplyr::filter(KEGG %in% mtb_info$cpd_id, !is.na(Metabolites))
anot1 <- rbind(neg_anot1, pos_anot1)
anot2 <- anot1 |> dplyr::filter(KEGG %in% mtb_info$cpd_id)

samples <- c("E95", "E115", "E135")

ord_e95 <- c(
  "Forebrain",
  "Midbrain",
  "Hindbrain",
  "Spinal cord",
  "Caudal neuropore",
  "Branchial arch",
  "Somite",
  "Heart",
  "AGM",
  "Mesenchyme",
  "Liver",
  "Lung",
  "Gut",
  "Cavity"
)
ord_e115 <- c(
  "Forebrain",
  "Midbrain",
  "Hindbrain",
  "Spinal cord",
  "Epidermis",
  "Branchial arch",
  "Jaw and tooth",
  "Forelimb",
  "Hindlimb",
  "Somite",
  "Heart",
  "AGM",
  "Liver",
  "Gut",
  "Cavity",
  "Embryo membrane"
)
ord_e135 <- c(
  "Forebrain",
  "Midbrain",
  "Hindbrain",
  "Diencephalon",
  "Spinal cord",
  "Epidermis",
  "Ear",
  "Jaw and tooth",
  "Hindlimb",
  "Cartilage",
  "Muscle",
  "Heart",
  "Blood vessel",
  "Gonad",
  "Kidney",
  "Liver",
  "Lung",
  "Gut"
)
ord_lst <- list(E95 = ord_e95, E115 = ord_e115, E135 = ord_e135)

rds_fn1 <-
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig2ac_mouse_mz_htp_dat_norm_lst.rds"
dat_norm_lst <- read_rds(rds_fn1)
for (samp in samples) {
  dat_norm = dat_norm_lst[[samp]][, ord_lst[[samp]]]
  lvls <- colnames(dat_norm)
  ## right anot ----
  nm <- rownames(dat_norm)
  tst1 = rbind(neg_anot1, pos_anot1) |> dplyr::filter(mz_id %in% nm) |>
    dplyr::select(all_of(c("mz_id", "Metabolites"))) |> dplyr::distinct() |>
    dplyr::group_by(mz_id) |>
    summarise(meta = paste0(Metabolites, collapse = "; "))
  mz_idx <- which(nm %in% tst1$mz_id)
  lab_mks <- tst1 |> as.data.frame() |> column_to_rownames("mz_id") |>
    .[nm[mz_idx], ] |> str_wrap(., width = 30)
  right_anot1 <- rowAnnotation(
    link = anno_mark(
      at = mz_idx,
      labels = lab_mks,
      labels_gp = gpar(fontsize = 4),
      link_width = unit(5, "mm"),
      link_height = unit(.2, "mm")
    )
  )
  top_anot11 <-
    ComplexHeatmap::HeatmapAnnotation(
      tissue = fct(lvls),
      show_legend = F,
      show_annotation_name = F,
      height = unit(1, 'mm'),
      width = unit(60, "mm"),
      col = list(tissue = cols_tissue)
    )

  col_fun2 <- colorRamp2(c(-3, 0, 3), c("lightblue", "gray100", "#aa3333"))

  htp11 <- ComplexHeatmap::Heatmap(
    as.matrix(dat_norm),
    top_annotation = top_anot11,
    show_column_names = T,
    show_row_names = F,
    row_names_gp = gpar(fontsize = 2),
    name = "z-score",
    show_row_dend = F,
    column_names_gp = gpar(fontsize = 6),
    show_column_dend = F,
    column_names_side = "top",
    column_names_rot = 30,
    right_annotation = right_anot1,
    cluster_columns = F,
    cluster_rows = F,
    col = col_fun2,
    height = unit(10, "cm"),
    width = unit(6, "cm")
  )
  pdf(
    glue("fig2ac_mz_top15_heatmap_{samp}_v2.pdf"),
    width = 4.5,
    height = 6
  )
  print(htp11)
  dev.off()
}

## fig2d: mouse selected tissue specific spatial distribution -----
rds_fn4 <-
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mmu_visium_obj_lst.rds"
mmu_visium_lst = read_rds(rds_fn4)

for (samp in c("E9.5", "E11.5", "E13.5")) {
  obj1 <- mmu_visium_lst[[samp]]
  sel_tissues <- c("CNS", "Heart", "Liver", "AGM")
  if (samp == "E9.5") {
    scale_fct <- 1.6
    plot_width <- 2
    plot_height <- 2
  } else if (samp == "E11.5") {
    scale_fct <- 1.1
    plot_width <- 4
    plot_height <- 4
  } else {
    scale_fct <- 2.0
    plot_width <- 6
    plot_height <- 6
    sel_tissues[4] <- "G and K"
  }
  plst1 <- lapply(sel_tissues, \(sel) {
    obj1$labels_plot <- case_when(obj1$labels_fig2d %in% sel ~ sel, TRUE ~ "others")
    Seurat::SpatialDimPlot(
      obj1,
      group.by = "labels_plot",
      label = F,
      pt.size = scale_fct,
      image.alpha = 0,
      stroke = NA,
      cols = c(cols_tissue[sel], "gray99")
    ) &
      my_theme1 & labs(title = sel) & Seurat::NoLegend() &
      theme(plot.title = element_text(size = 8),
            plot.margin = margin(c(0, 0, 0, 0), unit = "cm")) &
      Seurat::NoAxes() & coord_fixed()
  })
  fig2d <- plst1 |> patchwork::wrap_plots(nrow = 1)
  ggsave(
    glue("fig2d_mouse_spatial_tissue_{samp}_v250322.pdf"),
    fig2d,
    width = plot_width * 4,
    height = plot_height,
    units = "cm"
  )
}

## fig2g: metabolites changes among stages in liver -----
csv_fn1 <-
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig2g_mouse_mz_liver_long_fc.csv"
long_fc <- read_csv(csv_fn1) |> mutate(name = fct(as.character(name), levels = paste0("E", c(9.5, 11.5, 13.5))))
line1 = ggplot2::ggplot(long_fc, aes(
  x = name,
  y = value,
  group = feature,
  color = group
)) +
  ggplot2::geom_line() + scale_color_manual(values = unname(cols_tissue[1:6])) +
  ggplot2::facet_wrap( ~ group, ncol = 3) + my_theme1 +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = .5
    )
  ) +
  labs(x = "", y = "Relative abundance")
ggsave("fig2g_liver_mtb_line1.pdf",
       line1,
       width = 6,
       height = 3)

## fig2h: moran's I changes -----
csv_fn2 <-
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig2h_mouse_mz_moran_trends.csv"
outcsv <- read_csv(csv_fn2)
draw_line <- function(outcsv) {
  sigtissue <- outcsv  |>
    dplyr::select(feature, anot1) |>
    dplyr::distinct() |>
    group_by(anot1) |>
    summarise(n = n()) |>
    group_by(anot1) |>
    dplyr::filter(anot1 != "others") |>
    arrange(desc(n)) |>
    pull(anot1)
  outcsv |>
    dplyr::filter(anot1 != "others") |>
    mutate(anot1 = factor(anot1, levels = sigtissue)) |>
    ggplot(aes(
      x = name,
      y = value,
      group = feature,
      color = anot1
    )) +
    geom_line() + scale_color_manual(values = cols_tissue) +
    facet_wrap( ~ anot1) +
    ylab(label = "Moran's Index") +
    xlab(label = "Stage") +
    my_theme1 +
    theme(
      legend.position = "none",
      axis.line = element_blank(),
      panel.border = element_rect(
        linewidth = .5,
        fill = NA,
        color = "black"
      ),
      axis.text.x = element_text(
        angle = 90,
        hjust = 0.5,
        vjust = 0.5
      )
    )
}
pdf(glue("./fig2h_mz_moran_trend_tissue.pdf"))
draw_line(outcsv)
dev.off()

## fig2i: bar plot of each tissue ----
draw_bar <- function(outcsv) {
  outcsv  |>
    dplyr::select(feature, anot1) |>
    dplyr::distinct() |>
    group_by(anot1) |>
    summarise(n = n()) |>
    arrange(desc(n)) |>
    mutate(label = factor(anot1, levels = anot1)) |>
    ggplot(aes(x = label, y = n, fill = label)) +
    geom_bar(stat = "identity") +
    my_theme1 + scale_fill_manual(values = cols_tissue) +
    theme(legend.position = "none",
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0
          )) +
    labs(x = "", y = "Metabolite Number")
}
pdf(glue("./fig2i_mouse_mz_num_tissue.pdf"))
  draw_bar(outcsv)
dev.off()
