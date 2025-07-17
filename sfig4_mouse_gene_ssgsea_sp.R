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
  "viridis"
)
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/sfig4")
workdir |> fs::dir_create() |> setwd()

my_theme1 <- theme_classic(base_size = 8) +
  theme(
    legend.key.size = unit(3, "mm"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

## figS4a: visium data specific gene expression -----
rds_fn1 <- "../rds/mmu_visium_lst.rds"
seu_lst = read_rds(rds_fn1)

e95_genes1 <- c("Six3", "Trim54", "Bhmt")
e115_genes1 <- c("Dlx1", "Adprhl1", "Cyp2c70")
e135_genes1 <- c("Neurog2", "Myl7", "Kng2")
genes_lst1 <- list(E9.5 = e95_genes1,
                   E11.5 = e115_genes1,
                   E13.5 = e135_genes1)
plst1 <- lapply(paste0("E", c(9.5, 11.5, 13.5)), \(stg) {
  if (stg == "E9.5") {
    scale_fct <- 4.5
    plot_width <- 3
    plot_height <- 3
    sel_img <- "slice1"
  } else if (stg == "E11.5") {
    scale_fct <- 2.4
    plot_width <- 3
    plot_height <- 3
    sel_img <- "slice1.3"
  } else {
    scale_fct <- 4.5
    plot_width <- 3
    plot_height <- 3
    sel_img <- "slice1.2"
  }

  feats <- genes_lst1[[stg]]
  p1 <- Seurat::SpatialFeaturePlot(
    seu_lst[[stg]],
    features = feats,
    ncol = 3,
    pt.size.factor = scale_fct,
    images = sel_img
  ) & my_theme1 &
    Seurat::NoAxes() & coord_fixed()
  ggsave(
    glue("sfig4ac_mouse_gene_expr_{stg}.pdf"),
    p1,
    width = plot_width * 3,
    height = plot_height,
    unit = "in"
  )
})

## figS4d: visium data specific metabolic genes ----
e95_genes2 <- c("Gad2", "Kcnj5", "Bhmt")
e115_genes2 <- c("Kcnk10", "Kcnj5", "Apof")
e135_genes2 <- c("Gad1", "Cox7a1", "Miox")
genes_lst2 <- list(E9.5 = e95_genes2,
                   E11.5 = e115_genes2,
                   E13.5 = e135_genes2)
plst2 <- lapply(paste0("E", c(9.5, 11.5, 13.5)), \(stg) {
  if (stg == "E9.5") {
    scale_fct <- 3
    plot_width <- 3
    plot_height <- 3
    sel_img <- "slice1"
  } else if (stg == "E11.5") {
    scale_fct <- 1.6
    plot_width <- 3
    plot_height <- 3
    sel_img <- "slice1.3"
  } else {
    scale_fct <- 2
    plot_width <- 3
    plot_height <- 3
    sel_img <- "slice1.2"
  }

  feats <- genes_lst2[[stg]]
  Seurat::DefaultAssay(seu_lst[[stg]]) <- "SCT"
  p1 <- Seurat::SpatialFeaturePlot(
    seu_lst[[stg]],
    features = feats,
    ncol = 3,
    slot = "data",
    pt.size.factor = scale_fct,
    images = sel_img
  ) & my_theme1 &
    Seurat::NoAxes() & coord_fixed()
  ggsave(
    glue("sfig4df_mouse_mtb_gene_expr_{stg}.pdf"),
    p1,
    width = plot_width * 3,
    height = plot_height,
    unit = "in"
  )
})

## figS4g-i: ssgsea score of KEGG pathways -----
kegg_info <- read_csv("~/ref/kegg/mouse/kegg_mmu_all_pth_genes.csv")
kegg_pth_id <- kegg_info[, c(1, 3)] |> dplyr::distinct()
rds_fn2 <- "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/sfig4gi_ssgsea_score_df_lst1.rds"
df_lst1 <- read_rds(rds_fn2)
sel_pth1 <-
  dplyr::filter(
    kegg_pth_id,
    grepl(
      "(GABA|Glutamater|Taurine|Cardiac|Insulin sig|cAMP|Pentose pho|Purine|Pyrimi)",
      pth_name
    )
  )
p_lst1 <- lapply(1:nrow(sel_pth1), \(idx) {
  p_df1 <- df_lst1 |> dplyr::select(all_of(c(
    "imagerow", "imagecol", "stage", sel_pth1[[1]][idx]
  ))) |>
    dplyr::rename("pth" = sel_pth1[[1]][idx])
  p_df1 <- mutate(
    p_df1,
    imagecol = case_when(
      stage == "E11.5" ~ imagecol,
      stage == "E9.5" ~ imagecol,
      stage == "E13.5" ~ -1 * imagecol
    ),
    imagerow = case_when(
      stage == "E11.5" ~ -1 * imagerow,
      stage == "E9.5" ~ -1 * imagerow,
      stage == "E13.5" ~ imagerow
    )
  ) |>
    mutate(stage = fct(as.character(stage), levels = c("E9.5", "E11.5", "E13.5")))
  ggplot2::ggplot(p_df1, aes(x = imagecol, y = imagerow, color = pth)) +
    geom_point(size = .5) + viridis::scale_color_viridis() + facet_wrap( ~ stage, scale = "free") +
    my_theme1 + Seurat::NoAxes() + labs(title = sel_pth1[[2]][idx], color = "")
})
pdf("sfig4gi_ssgsea_score.pdf",
    width = 6,
    height = 3)
print(p_lst1)
dev.off()

## figS9d: signaling score ----
sel_pth2 <-
  dplyr::filter(kegg_pth_id,
                grepl("(Wnt sig|Hedgehog|TGF|Notch sig|HIF|mTOR)", pth_name))
p_lst2 <- lapply(1:nrow(sel_pth2), \(idx) {
  p_df1 <- df_lst1 |> dplyr::select(all_of(c(
    "imagerow", "imagecol", "stage", sel_pth2[[1]][idx]
  ))) |>
    dplyr::rename("pth" = sel_pth2[[1]][idx])
  p_df1 <- mutate(
    p_df1,
    imagecol = case_when(
      stage == "E11.5" ~ imagecol,
      stage == "E9.5" ~ imagecol,
      stage == "E13.5" ~ -1 * imagecol
    ),
    imagerow = case_when(
      stage == "E11.5" ~ -1 * imagerow,
      stage == "E9.5" ~ -1 * imagerow,
      stage == "E13.5" ~ imagerow
    )
  ) |>
    mutate(stage = fct(as.character(stage), levels = c("E9.5", "E11.5", "E13.5")))
  ggplot2::ggplot(p_df1, aes(x = imagecol, y = imagerow, color = pth)) +
    geom_point(size = .5) + viridis::scale_color_viridis() + facet_wrap( ~ stage, scale = "free") +
    my_theme1 + Seurat::NoAxes() + labs(title = sel_pth2[[2]][idx], color = "")
})
pdf("sfig9df_ssgsea_score.pdf",
    width = 6,
    height = 3)
print(p_lst2)
dev.off()
