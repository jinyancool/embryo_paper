pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "paletteer", "cowplot", "ComplexHeatmap", "circlize", "plot1cell")  
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/figS2")
workdir %>% fs::dir_create() %>% setwd()

yaml_fn <- "/cluster/home/danyang_jh/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn = yaml_fn, "stage")[c("CS12", "CS14", "CS18")]
## figS2c: plot1cell of all genes ----
rds_fn6 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/human_visium_all_gene_lst.rds"
seu_new2_lst = read_rds(rds_fn6)
for(idx in 1:length(seu_new2_lst)){
  circ_dat <- prepare_circlize_data(seu_new2_lst[[idx]], scale = .65)
  clust_cols = cols_tissue[sort(unique(circ_dat$tissue))]
  lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                      legend_gp = gpar(fill = stg_cols), 
                                      title_position = "topleft", title = "stage")
  pdf(glue("figS2c_plot1cell_all_genes_{names(seu_new2_lst)[idx]}.pdf"), width = 8, height = 8)
  plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
                col.use = clust_cols ,bg.color = 'white',
                kde2d.n = 200, repel = T, label.cex = .9)
  add_track(circ_dat, group = "stage", colors = stg_cols[unique(seu_new2_lst[[idx]]$stage)], track_num = 2) 
  draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
  dev.off()
}

## figS2c: plot1cell of mtb genes ----
rds_fn7 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/human_visium_all_gene_lst.rds"
seu_new3_lst = read_rds(rds_fn7)
for(idx in 1:length(seu_new3_lst)){
  circ_dat <- prepare_circlize_data(seu_new3_lst[[idx]], scale = .65)
  clust_cols = cols_tissue[sort(unique(circ_dat$tissue))]
  lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                      legend_gp = gpar(fill = stg_cols), 
                                      title_position = "topleft", title = "stage")
  pdf(glue("figS3_human_plot1cell_mtb_genes_{unique(seu_new3_lst[[idx]]$stage)}.pdf"), width = 8, height = 8)
  plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
                col.use = clust_cols ,bg.color = 'white',
                kde2d.n = 200, repel = T, label.cex = .6)
  
  add_track(circ_dat, group = "stage", colors = stg_cols[unique(seu_new3_lst[[idx]]$stage)], track_num = 2) 
  
  draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
  dev.off()
}

## metabolic genes of each sample ---- 
rds_fn7 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/hsa_visium_mtb_gene_lst.rds"
seu_new3_lst = read_rds(rds_fn7)
for(idx in 1:length(seu_lst)){
  circ_dat <- prepare_circlize_data(seu_new3_lst[[idx]], scale = .65)
  clust_cols = cols_tissue[sort(unique(circ_dat$tissue))]
  lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid",
                                      legend_gp = gpar(fill = stg_cols),
                                      title_position = "topleft", title = "stage")
  pdf(glue("figS2c_human_plot1cell_mtb_genes_{names(seu_lst)[idx]}_v2.pdf"), width = 8, height = 8)
  plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
                col.use = clust_cols ,bg.color = 'white',
                kde2d.n = 200, repel = T, label.cex = .6)

  add_track(circ_dat, group = "stage", colors = stg_cols[unique(seu_lst[[idx]]$stage)], track_num = 2)

  draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
  dev.off()
}


