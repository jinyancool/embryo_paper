pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "viridis")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/sfig12")
workdir %>% fs::dir_create() %>% setwd()

yaml_fn <- "/cluster/home/danyang_jh/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn = yaml_fn, "stage")[c("CS12", "CS14", "CS18")]

my_theme1 <- theme_classic(base_size = 8) + 
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"), 
        axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"))

## figS12a: RCTD results1, main types distribution -----
rds_fn1 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/sfig12a_rctd_res_df_lst1.rds"
df_lst1 <- read_rds(rds_fn1)
lapply(c("E9.5", "E11.5", "E13.5"), \(stg) {
  if(stg == "E13.5") {
    scl = .5
    plot_width = 4
    plot_height = 4
  } else if (stg == "E11.5") {
    scl = .5
    plot_width = 4
    plot_height = 3
  } else {
    scl = .8
    plot_width = 3
    plot_height = 2.5
  }
  hcc2t_ptn_ratio <- df_lst1[[stg]]
  pie_p <- ggplot2::ggplot() +
    scatterpie::geom_scatterpie(aes(x = row, y = col), pie_scale = scl, col = NA, data = hcc2t_ptn_ratio, 
                                cols = colnames(hcc2t_ptn_ratio)[-c(1:3)], long_format = F) +
    ggplot2::theme_void(base_size = 8) + 
    ggsci::scale_fill_igv() + coord_fixed() + 
    theme(legend.key.size = unit(3, "mm"), legend.position = "right") 
  ggsave(glue("sfig12a_rctd_prop_pie_{stg}.pdf"), pie_p, width = plot_width, height = plot_height)
})

## figS12b: RCTD results1, detailed types distribution in liver region -----
rds_fn2 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/sfig12b_liver_rctd_res.rds"
rctd_dat <- read_rds(rds_fn2)
pie_p <- ggplot2::ggplot() +
  scatterpie::geom_scatterpie(aes(x = row, y = col), pie_scale = 1, col = NA, data = rctd_dat, 
                              cols = colnames(rctd_dat)[-c(1:3)], long_format = F) +
  ggplot2::theme_void(base_size = 6) + 
  ggsci::scale_fill_igv(alpha = .8) + 
  theme(legend.key.size = unit(3, "mm"), legend.position = "right") 
ggsave(glue("sfig12b_liver_rctd_prop_pie.pdf"), pie_p, width = 4, height = 3)

## figS12c: RCTD weights, liver of each stage -----
rds_fn3 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/sfig12c_visium_lst4.rds"
visium_lst4 <- read_rds(rds_fn3)
p_lst1 <- list()
for(feat in sel_types) {
  plst1 <- lapply(visium_lst4, \(seu) {
    if(grepl("E9.5", unique(seu$sample))) {
      pt_size = 4
    } else if (unique(seu$sample) == "E11.5_S3E2") {
      pt_size = 3
    } else if (unique(seu$sample) == "E11.5_S3E2_2"){
      pt_size = 6
    } else if (unique(seu$sample) == "E13.5_S1E1_1st") {
      pt_size = 3.5
    } else if (unique(seu$sample) == "E13.5_S1E1_2nd") {
      pt_size = 4.8
    }
    feat1 <- Seurat::SpatialFeaturePlot(seu[, !is.na(seu[[feat]])], image.alpha = 0.5, features = feat, 
                                        pt.size.factor = pt_size, stroke = NA, combine = F, crop = T)
    feat1 <- feat1 %>% lapply(., \(p) p + Seurat::DarkTheme() + NoGrid() + NoAxes()) %>% 
      patchwork::wrap_plots() &
      paletteer::scale_fill_paletteer_c("grDevices::Viridis") & 
      theme(legend.position = "right", legend.key.size = unit(3, "mm"), 
            text = element_text(size = 6)) & 
      labs(title = feat, fill = "proportion", x = "", y = "")
  })
  p_lst1[[feat]] <- patchwork::wrap_plots(plst1, nrow = 1)
}
pdf("sfig12c_rctd_proportion_spatial_feats.pdf", height = 2, width = 8)
print(p_lst1)
dev.off()

