pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "paletteer", "cowplot", "ComplexHeatmap", "circlize", "parallel", "monocle")  
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/fig6")
workdir %>% fs::dir_create() %>% setwd()

yaml_fn <- "/cluster/home/danyang_jh/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn = yaml_fn, "stage")[c("E9.5", "E11.5", "E13.5")]

my_theme1 <- theme_classic(base_size = 8) + 
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"), 
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = .5))
my_theme2 <- theme_classic(base_size = 8) + theme(legend.key.size = unit(3, 'mm')) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = .5), 
        axis.ticks = element_blank(), axis.title = element_blank(), 
        panel.grid = element_blank(), panel.border = element_rect(linewidth = .5, fill = NA)) 

## fig6b: monocle2 of gene and m/z data ----
rds_fn1 = 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/align_new/mtb_new/pseudotime_2502011/seu1_brain2.rds"
seu1_brain2 = read_rds(rds_fn1)
feat1 = Seurat::SpatialFeaturePlot(seu1_brain2, features = "Pseudotime", 
                                   pt.size.factor = 1.3, image.alpha = 1) & 
  my_theme1 & coord_fixed() & viridis::scale_fill_viridis()
ggsave("fig6b_gene_pseudotime_spatial.pdf", feat1, width = 4, height = 2)
rds_fn2 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/align_new/mtb_new/pseudotime_2502011/monocle2_cds_vld_2e-1.rds"
cds_vld_gene <- read_rds(rds_fn2)
cds_vld_gene$Pseudotime <- seu1_brain2$Pseudotime
p_traj1 <- plot_cell_trajectory(cds_vld_gene, color_by = "Pseudotime", cell_size = 5e-1) + 
  theme_classic(base_size = 8) + viridis::scale_color_viridis() + 
  theme(legend.position = "right", legend.key.size = unit(3, "mm")) + ggplot2::coord_fixed()
ggsave(glue::glue("fig6b_pseudotime_monocle2_gene.pdf"), p_traj1, 
       width = 4, height = 2)

### cds object, m/z based -----
rds_fn3 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/align_new/mtb_new/pseudotime_2502011/cds_vld_mz.rds"
cds_vld_mz <- read_rds(rds_fn3)
p_traj2 <- plot_cell_trajectory(cds_vld_mz, color_by = "Pseudotime", cell_size = 5e-1) + 
  theme_classic(base_size = 8) + viridis::scale_color_viridis() + 
  theme(legend.position = "right", legend.key.size = unit(3, "mm")) + ggplot2::coord_fixed()
ggsave(glue::glue("fig6b_pseudotime_monocle2_mz.pdf"), p_traj2, 
       width = 4, height = 2)

rds_fn4 = 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig6b_mz_seu_obj.rds"
mz_brain = read_rds(rds_fn4)
feat2 = Seurat::SpatialFeaturePlot(mz_brain, features = "Pseudotime", 
                                   pt.size.factor = 1.3) & 
  my_theme1 & coord_fixed() & viridis::scale_fill_viridis()
ggsave("fig6b_mz_pseudotime_spatial.pdf", feat2, width = 4, height = 2)

## fig6g: signaling and mtb pathway co-occurrence -----
comp_df <- tibble(signaling = "Notch signaling pathway", 
                  metabolic = c("Glycolysis / Gluconeogenesis", "Alanine, aspartate and glutamate metabolism", 
                                "Sphingolipid metabolism"))
rds_fn3 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig4gi_obj_lst_v2.rds"
obj_lst <- read_rds(rds_fn3)
plst3 <- list()
for(samp in c("E115", "E135")) {
  obj1 <- obj_lst[[samp]]
  if(samp == "E115") {
    cord_df = Seurat::GetTissueCoordinates(obj1) %>% .[colnames(obj1), ] %>% tibble() %>% 
      mutate(x = imagecol, y = -1 * imagerow)
  } else {
    cord_df = Seurat::GetTissueCoordinates(obj1) %>% .[colnames(obj1), ] %>% tibble() %>% 
      mutate(x = -1 * imagecol, y = imagerow)
    
  }
  
  cord1 <- cbind(cord_df, obj1@meta.data)
  
  plst3[[samp]] <- lapply(1:nrow(comp_df), \(idx) {
    feat1 <- comp_df[["signaling"]][idx] %>% as.character() 
    feat2 <- comp_df[["metabolic"]][idx] %>% as.character() 
    df4p <- cord1 %>% dplyr::select(all_of(c("x", "y", feat1, feat2))) %>% 
      dplyr::rename("feat1" = feat1, "feat2" = feat2) %>% 
      mutate(top = case_when((feat1 > quantile(.$feat1, .85)) & (feat2 > quantile(.$feat2, .85)) ~ "top 15%", 
                             TRUE ~ "no"))
    p <- ggplot2::ggplot() + 
      ggplot2::geom_point(
        data = df4p, aes(x = x, y = y, color = feat1), show.legend = F, 
        alpha = .5, size = .2
      ) + labs(color = feat1) + 
      scale_color_continuous(low = "gray90", high = "#ff808f") +
      my_theme2 + coord_fixed() + labs(color = "signaling") + 
      ggnewscale::new_scale_color() + 
      ggplot2::geom_point(
        data = df4p, aes(color = feat2, x = x, y = y), show.legend = F, 
        alpha = .5, size = .2) + 
      paletteer::scale_color_paletteer_c("pals::kovesi.linear_bgy_10_95_c74") + 
      my_theme2 + coord_fixed() + labs(color = "metabolic") + 
      ggnewscale::new_scale_color() + 
      geom_point(data = dplyr::filter(df4p, top == "top 15%"), 
                 mapping = aes(x = x, y = y, color = top), 
                 size = .3, color = "#f20020") + coord_fixed()
    p1 <- p + 
      ggnewscale::new_scale_fill() +
      ggplot2::stat_density_2d_filled(
        data = df4p %>% dplyr::filter(top == "top 15%"),
        mapping = aes(fill = ..ndensity.., #alpha = ..ndensity.., 
                      x = x, y = y), alpha = .5, 
        geom = "raster", contour = F, show.legend = T
      ) + 
      ggplot2::scale_fill_gradientn(colors = c("white", "#F19E62FF", "yellow")) + 
      coord_fixed() + theme(legend.position = "right")
    p2 = p1 + 
      ggnewscale::new_scale_fill() +
      geom_density_2d(
        aes(x = x, y = y), 
        data = df4p %>% dplyr::filter(top == "top 15%"),
        contour_var	= "ndensity", alpha = .8, 
        show.legend = T, linewidth = .2
      ) + 
      coord_fixed() + labs(title = glue("{feat1}\n&\n{feat2}")) + 
      theme(plot.title = element_text(hjust = .5))
    
  })
}
fig6g <- (plst3[[1]][[1]] | plst3[[2]][[1]])/(plst3[[1]][[2]] | plst3[[2]][[2]]) / (plst3[[1]][[3]] | plst3[[2]][[3]])
pdf(glue("fig6g_e115_e135_signal_mz_co_expr.pdf"), width = 7, height = 11)
print(fig6g)
dev.off()



