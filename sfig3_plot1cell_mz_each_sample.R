pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "paletteer", "cowplot", "ComplexHeatmap", "circlize", "plot1cell")  
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/figS3")
workdir |> fs::dir_create() |> setwd()

yaml_fn <- "~/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
cols_stg <- jhtools::show_me_the_colors(config_fn= yaml_fn, "stage")

## figS3a: plot1cell of mouse m/z merged data -----
rds_fn1 <- 
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mouse_mz_mrg_lst.rds"
mz_mrg_lst <- read_rds(rds_fn1)

stg_cols <- cols_stg[names(mz_mrg_lst)]
for(idx in 1:length(mz_mrg_lst)){
  circ_dat <- prepare_circlize_data(mz_mrg_lst[[idx]], scale = .7)
  circ_dat <- circ_dat[!is.na(rownames(circ_dat)), ]
  clust_cols = cols_tissue[sort(unique(circ_dat$tissuetype))]
  lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                      legend_gp = gpar(fill = stg_cols), 
                                      title_position = "topleft", title = "stage")
  pdf(glue("figS3_plot1cell_all_mz_{names(mz_mrg_lst)[idx]}.pdf"), width = 5.5, height = 5.5)
  plot_circlize(circ_dat, do.label = T, pt.size = 0.4, 
                col.use = clust_cols ,bg.color = 'white', contour.nlevels = 100, 
                kde2d.n = 2000, repel = T, label.cex = .7)
  
  add_track(circ_dat, group = "stage", colors = stg_cols[names(mz_mrg_lst)[idx]], track_num = 2) 
  
  draw(lgd_squre, x = unit(25, "mm"), y = unit(8, "mm"), just = c("right", "bottom"))
  dev.off()
}

rds_fn2 <- "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mouse_mz_obj_merged.rds"
mz_obj_mrg = read_rds(rds_fn2)
circ_dat <- prepare_circlize_data(mz_obj_mrg, scale = .7)
clust_cols = cols_tissue[sort(unique(circ_dat$tissuetype))]
lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                    legend_gp = gpar(fill = stg_cols), 
                                    title_position = "topleft", title = "stage")
pdf(glue("fig1c_plot1cell_all_mz_mrg.pdf"), width = 5.5, height = 5.5)
plot_circlize(circ_dat, do.label = T, pt.size = 0.4, 
              col.use = clust_cols ,bg.color = 'white', contour.nlevels = 100, 
              kde2d.n = 2000, repel = T, label.cex = .7)

add_track(circ_dat, group = "stage", colors = stg_cols[sort(unique(circ_dat$stage))], track_num = 2) 

draw(lgd_squre, x = unit(25, "mm"), y = unit(8, "mm"), just = c("right", "bottom"))
dev.off()

## figS3c: plot1cell of human m/z merged data -----
rds_fn3 <- "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/human_mz_mrg_lst.rds"
mz_mrg_lst = read_rds(rds_fn3)
for(idx in c("yao1", "yao2", "yao5")){
  obj1 <- mz_mrg_lst[[idx]][, !is.na(mz_mrg_lst[[idx]]$tissue)]
  Idents(obj1) <- obj1$tissue
  circ_dat <- prepare_circlize_data(obj1, scale = .7)
  circ_dat <- circ_dat[!is.na(rownames(circ_dat)), ]
  clust_cols = cols_tissue[sort(unique(circ_dat$tissue))]
  lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid",
                                      legend_gp = gpar(fill = stg_cols),
                                      title_position = "topleft", title = "stage")
  pdf(glue("figS3c_plot1cell_all_mz_{unique(obj1$stage)}.pdf"), width = 5.5, height = 5.5)
  plot_circlize(circ_dat, do.label = T, pt.size = 0.4,
                col.use = clust_cols ,bg.color = 'white', contour.nlevels = 100,
                kde2d.n = 2000, repel = T, label.cex = .7)

  add_track(circ_dat, group = "stage", colors = stg_cols[unique(obj1$stage)], track_num = 2)

  draw(lgd_squre, x = unit(20, "mm"), y = unit(8, "mm"), just = c("right", "bottom"))
  dev.off()
}



