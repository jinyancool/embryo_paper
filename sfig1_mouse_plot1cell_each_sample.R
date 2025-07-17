pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "paletteer", "cowplot", "ComplexHeatmap", "circlize", "plot1cell")  
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/figS1")
workdir |> fs::dir_create() |> setwd()

yaml_fn <- "~/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn= yaml_fn, "stage") |> .[c("E9.5", "E11.5", "E13.5")] 

## figS1c: plot1cell, umap + outer circle -----
mmu_visium_lst <- list()
for(idx in c("E9.5", "E11.5", "E13.5")){
  # intersect_mtb_genes <- intersect(rownames(seu_lst3[[idx]]), mmu_mtb_genes)
  obj1 <- seu_lst3[[idx]]
  obj1 <- obj1 |> Seurat::SCTransform(assay = "Spatial") |>
    Seurat::RunPCA(npcs = 100) |> Seurat::RunUMAP(dims = 1:30) |> Seurat::RunTSNE(dims = 1:10) |>
    Seurat::FindNeighbors(dims = 1:30)
  Idents(obj1) <- obj1$tissuetype
  mmu_visium_lst[[idx]] <- obj1
}
## figS1c: all genes of each sample ----
rds_fn2 <- "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mmu_visium_lst.rds"
mmu_visium_lst <- read_rds(rds_fn2)
for(idx in c("E9.5", "E11.5", "E13.5")) {
  obj1 <- mmu_visium_lst[[idx]]
  circ_dat <- prepare_circlize_data(obj1, scale = .65)
  clust_cols <- cols_tissue[order(unique(circ_dat$tissuetype))] # caution: the color must be orderd by the names, as a-z
  lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid",
                                      legend_gp = gpar(fill = stg_cols),
                                      title_position = "topleft", title = "stage")
  pdf(glue("figS1c_plot1cell_all_genes_{idx}.pdf"), width = 8, height = 8)
  plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
                col.use = clust_cols ,bg.color = 'white',
                kde2d.n = 200, repel = T, label.cex = .9)
  
  add_track(circ_dat, group = "stage", colors = stg_cols[idx], track_num = 2)
  draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
  dev.off()
}
## figS1c: metabolic genes of each sample ----
rds_fn3 <- 
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mmu_visium_mtb_gene_obj_lst.rds"
seu_new3_lst <- read_rds(rds_fn3)
for(idx in 1:length(seu_new3_lst)) {
  circ_dat <- prepare_circlize_data(seu_new3_lst[[idx]], scale = .65)
  clust_cols = cols_tissue[sort(unique(circ_dat$tissuetype))]
  lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                      legend_gp = gpar(fill = stg_cols), 
                                      title_position = "topleft", title = "stage")
  pdf(glue("figS1c_plot1cell_mtb_genes_{names(seu_lst3)[idx]}.pdf"), width = 8, height = 8)
  plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
                col.use = clust_cols ,bg.color = 'white',
                kde2d.n = 200, repel = T, label.cex = .6)
  
  add_track(circ_dat, group = "stage", colors = stg_cols[names(seu_lst3)[idx]], track_num = 2) 
  
  draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
  dev.off()
}

