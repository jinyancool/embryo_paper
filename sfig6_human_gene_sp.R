pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "viridis")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/sfig6")
workdir %>% fs::dir_create() %>% setwd()

my_theme1 <- theme_classic(base_size = 8) + 
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"), 
        axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"))

## figS5a: visium data specific gene expression -----
rds_fn1 <- "../rds/hsa_visium_obj_lst.rds"
seu_lst = read_rds(rds_fn1)

yao1_genes1 <- c("PAX5", "LRRC10", "SULT1E1")
yao2_genes1 <- c("FOXL2", "LRRC10", "CYP3A7")
yao5_genes1 <- c("FOXL2", "NKX2-5", "SLC17A1")
genes_lst1 <- list(yao1 = yao1_genes1, yao2 = yao2_genes1, yao5 = yao5_genes1)
plst1 <- lapply(paste0("yao", c(1, 2, 5)), \(stg) {
  if(stg == "yao1") {
    scale_fct <- 3
    plot_width <- 3
    plot_height <- 3
    # sel_img <- "slice1"
  } else if (stg == "yao2") {
    scale_fct <- 1.8
    plot_width <- 3
    plot_height <- 3
    # sel_img <- "slice1.3"
  } else {
    scale_fct <- 1.2
    plot_width <- 3
    plot_height <- 3
    # sel_img <- "slice1.2"
  }
  
  feats <- genes_lst1[[stg]]
  p1 <- Seurat::SpatialFeaturePlot(seu_lst[[stg]], features = feats, ncol = 3, 
                                   pt.size.factor = scale_fct, stroke = NA, #images = sel_img
                                   ) & my_theme1 & 
    Seurat::NoAxes() & coord_fixed()
  ggsave(glue("sfig6ac_human_gene_expr_{stg}.pdf"), p1, 
         width = plot_width * 3, height = plot_height, unit = "in")
})

## figS5d: visium data specific metabolic genes ----
yao1_genes2 <- c("NOS2", "COX6A2", "SLC13A5")
yao2_genes2 <- c("SLC6A1", "COX6A2", "KLB")
yao5_genes2 <- c("SLC4A4", "SLC5A1", "SLC23A1")
genes_lst2 <- list(yao1 = yao1_genes2, yao2 = yao2_genes2, yao5 = yao5_genes2)
plst2 <- lapply(paste0("yao", c(1, 2, 5)), \(stg) {
  if(stg == "yao1") {
    scale_fct <- 3
    plot_width <- 3
    plot_height <- 3
    # sel_img <- "slice1"
  } else if (stg == "yao2") {
    scale_fct <- 1.8
    plot_width <- 3
    plot_height <- 3
    # sel_img <- "slice1.3"
  } else {
    scale_fct <- 1.2
    plot_width <- 3
    plot_height <- 3
    # sel_img <- "slice1.2"
  }
  
  feats <- genes_lst2[[stg]]
  Seurat::DefaultAssay(seu_lst[[stg]]) <- "SCT"
  p1 <- Seurat::SpatialFeaturePlot(seu_lst[[stg]], features = feats, ncol = 3, slot = "data", 
                                   pt.size.factor = scale_fct, stroke = NA, #images = sel_img
                                   ) & my_theme1 & 
    Seurat::NoAxes() & coord_fixed()
  ggsave(glue("sfig6df_human_mtb_gene_expr_{stg}.pdf"), p1, 
         width = plot_width * 3, height = plot_height, unit = "in")
})

## ssgsea score re-scale -----

