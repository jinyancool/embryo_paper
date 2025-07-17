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
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/fig3")
workdir |> fs::dir_create() |> setwd()

yaml_fn <- "~/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn = yaml_fn, "tissue")

my_theme1 <- theme_classic(base_size = 8) +
  theme(
    legend.key.size = unit(3, "mm"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = .5)
  )

## fig3a-c: human m/z tissue specific markers -----
mtb_info <- read_rds("~/ref/kegg/human/hsa_mtb_pth_cpd.rds")
xlsx_fn1 <-
  glue(
    "/cluster/home/jhuang/projects/collabrators/data/wangwenjie/human/metabolism/human_adjusted/DZLM2024030584-b2_DZLM2024030572/Qualitative.xlsx"
  )
neg_anot1 <- readxl::read_excel(xlsx_fn1, sheet = "neg-all") |>
  mutate(mz_id = paste0("neg-", mz)) |> dplyr::filter(KEGG %in% mtb_info$cpd_id, !is.na(Metabolites))
pos_anot1 <- readxl::read_excel(xlsx_fn1, sheet = "pos-all") |>
  mutate(mz_id = paste0("pos-", mz)) |> dplyr::filter(KEGG %in% mtb_info$cpd_id, !is.na(Metabolites))

samples <- c("yao1", "yao2", "yao5")

ord_yao1 <- c(
  "Forebrain",
  "Midbrain",
  "Hindbrain",
  "Spinal cord",
  "Optic vesicle",
  "Jaw and tooth",
  "Branchial arch",
  "Somite",
  "Heart",
  "Liver",
  "Gut",
  "Embryo membrane"
)
ord_yao2 <- c(
  "Forebrain",
  "Midbrain",
  "Hindbrain",
  "Spinal cord",
  "Optic vesicle",
  "Jaw and tooth",
  "Branchial arch",
  "Somite",
  "Heart",
  "AGM",
  "Liver",
  "Lung",
  "Gut",
  "Umbilical cord"
)
ord_yao5 <- c(
  "Forebrain",
  "Midbrain",
  "Hindbrain",
  "Diencephalon",
  "Spinal cord",
  "Ear",
  "Jaw and tooth",
  "Forelimb",
  "Hindlimb",
  "Cartilage",
  "Muscle",
  "Heart",
  "Blood vessel",
  "Liver",
  "Gut"
)
ord_lst <- list(yao1 = ord_yao1, yao2 = ord_yao2, yao5 = ord_yao5)

rds_fn3 <-
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig3ac_human_mz_dat_norm_lst.rds"
dat_norm_lst <- read_rds(rds_fn3)
for (samp in samples) {
  dat_norm = dat_norm_lst[[samp]]
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
      labels_gp = gpar(fontsize = 6),
      link_width = unit(3, "mm"),
      link_height = unit(.05, "mm")
    )
  )
  top_anot11 <-
    ComplexHeatmap::HeatmapAnnotation(
      tissue = fct(lvls),
      show_legend = F,
      show_annotation_name = F,
      col = list(tissue = cols_tissue)
    )

  col_fun2 <- colorRamp2(c(-3, 0, 3), c("lightblue", "gray100", "#aa3333"))

  htp11 <- ComplexHeatmap::Heatmap(
    as.matrix(dat_norm),
    top_annotation = top_anot11,
    show_column_names = T,
    show_row_names = F,
    row_names_gp = gpar(fontsize = 3),
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
    glue("fig3ac_mz_top15_heatmap_{samp}_v250322.pdf"),
    width = 5,
    height = 5
  )
  print(htp11)
  dev.off()
}

## fig3d: human selected tissue specific spatial distribution -----
rds_fn4 <-
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/human_visium_all_gene_lst.rds"
hsa_visium_lst = read_rds("../rds/human_visium_all_gene_lst.rds")
names(hsa_visium_lst) <- c("CS12", "CS14", "CS18")

for (samp in c("CS12", "CS14", "CS18")) {
  obj1 <- hsa_visium_lst[[samp]]
  if (samp == "CS12") {
    scale_fct <- 1.8
    plot_width <- 3
    plot_height <- 3
  } else if (samp == "CS14") {
    scale_fct <- 1.3
    plot_width <- 4
    plot_height <- 4
  } else {
    scale_fct <- 1.3
    plot_width <- 8
    plot_height <- 8
  }
  plst1 <- lapply(c("CNS", "Heart", "Liver", "Somite"), \(sel) {
    obj1$labels_plot <- case_when(obj1$labels_fig3d %in% sel ~ sel, TRUE ~ "others")
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
  plst1 <- plst1 |> patchwork::wrap_plots(nrow = 1)
  ggsave(
    glue("fig3d_human_spatial_tissue_{samp}_v250322.pdf"),
    plst1,
    width = plot_width * 4,
    height = plot_height,
    units = "cm"
  )
}

## fig3g: changes among stages of human metabolites in liver -----
csv_fn1 <- "/cluster/home/ztao_jh/projects/embryo/analysis/zhangjing/human/rnaseq/new_label/ren_time_point_2/neg/Liver.csv"
log_fc <- read_csv(csv_fn1)
log_fc <- log_fc |>
  mutate(stage = case_when(
    grepl("yao1", name) ~ "CS12",
    grepl("yao2", name) ~ "CS14",
    TRUE ~ "CS18"
  )) |>
  mutate(stage = fct(as.character(stage), levels = paste0("CS1", c(2, 4, 8)))) |>
  dplyr::filter(!grepl("Trend7", group))

line1 = ggplot2::ggplot(log_fc, aes(
  x = stage,
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
ggsave("fig3g_liver_mtb_line1.pdf",
       line1,
       width = 6,
       height = 3)

## fig2h: moran's I changes -----
csv_fn2 <-
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig3h_human_mz_moran_trends.csv"
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
      x = stage,
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
pdf(glue("./fig3h_mz_moran_trend_tissue.pdf"))
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
            vjust = .5
          )) +
    labs(x = "", y = "Metabolite Number")
}
pdf(glue("./fig3i_human_mz_num_tissue.pdf"))
draw_bar(outcsv)
dev.off()
