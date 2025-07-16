pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "paletteer", "cowplot", "ComplexHeatmap", "circlize", "parallel")  
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/fig4")
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

# fig4a: mouse regulon activity of scenic -----
write_csv(mouseobj_full_frame, "../rds/fig4a_tf_mouseobj_full_frame.csv")
csv_fn1 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig4a_tf_mouseobj_full_frame.csv"
mouseobj_full_frame <- read_csv(csv_fn1)
setis <- c("HNF4A(+)", "NEUROD2(+)", "GATA4(+)") %>% stringr::str_to_sentence()
parallel::mclapply(setis, function(seti){
  p <- mouseobj_full_frame %>%
    dplyr::filter(name == {{seti}}) %>%
    ggplot(aes(x = imagecol, y = imagerow, color = value)) +
    geom_point(size = 0.1) +
    scale_color_viridis_c() +
    my_theme1 + 
    facet_wrap(~ tissue, ncol = 3) +
    coord_fixed() +
    ggtitle(seti) +
    Seurat::DarkTheme() + NoGrid() + NoAxes()
  ggsave(glue("./fig4a_mouse_tf_{seti}_v250322.pdf"), width = 7, height = 3)
}, mc.cores = 3)

# fig4b: human regulon activity of scenic -----
setis <- c("HNF4A(+)", "NEUROD2(+)", "GATA4(+)")
csv_fn2 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig4b_tf_humanobj_full_frame.csv"
humanobj_full_frame = read_csv(csv_fn2)
parallel::mclapply(setis, function(seti){
  p <- humanobj_full_frame %>%
    dplyr::filter(name == {{seti}}) %>%
    ggplot(aes(x = imagecol, y = imagerow, color = value)) +
    geom_point(size = 0.01) +
    scale_color_viridis_c() +
    my_theme1 + 
    facet_wrap(~ tissue, nrow = 1) +
    coord_fixed() + 
    ggtitle(seti) + 
    Seurat::DarkTheme() + NoGrid() + NoAxes()
  ggsave(glue("./fig4b_human_tf_{seti}_v250322.pdf"), width = 7, height = 3, unit = "in")
}, mc.cores = 3)

# fig4g-i: co-occurrence of regulons and metabolic pathway scores, in mouse -----
tbl1 <- tibble(tf = "Hnf4a", metabolic = c("Purine metabolism", "Pyrimidine metabolism"))
tbl2 <- tibble(tf = "Sox2", metabolic = c("Alanine, aspartate and glutamate metabolism", "Sphingolipid metabolism"))
tbl3 <- tibble(tf = "Gata4", metabolic = c("Biosynthesis of unsaturated fatty acids", "Citrate cycle (TCA cycle)"))
comp_df <- bind_rows(tbl1, tbl2, tbl3)

rds_fn3 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig4gi_obj_lst_v2.rds"
obj_lst <- read_rds(rds_fn3)
obj1 <- obj_lst[["E115"]]
cord_df = Seurat::GetTissueCoordinates(obj1) %>% .[colnames(obj1), ] %>% tibble() %>% 
  mutate(x = imagecol, y = -1 * imagerow)
cord1 <- cbind(cord_df, obj1@meta.data)

plst3 <- lapply(1:nrow(comp_df), \(idx) {
  feat1 <- comp_df[["tf"]][idx] %>% as.character() %>% gsub(" / ", "_", x = .)
  feat2 <- comp_df[["metabolic"]][idx] %>% as.character() %>% gsub(" / ", "_", x = .)
  df4p <- cord1 %>% dplyr::select(all_of(c("x", "y", feat1, feat2))) %>% 
    dplyr::rename("feat1" = feat1, "feat2" = feat2) %>% 
    mutate(top = case_when((feat1 > quantile(.$feat1, .85)) & (feat2 > quantile(.$feat2, .85)) ~ "top 15%", 
                           TRUE ~ "no"))
  sp1 <- ggplot2::ggplot() + ggplot2::geom_point(
    data = df4p, aes(x = x, y = y, color = feat1), show.legend = F, 
    alpha = 1, size = .2
  ) + labs(color = feat1) + 
    scale_color_continuous(low = "gray90", high = "#ff808f") +
    my_theme2 + coord_fixed() + 
    labs(title = feat1, color = "")
  sp2 <- ggplot2::ggplot() + ggplot2::geom_point(
    data = df4p, aes(x = x, y = y, color = feat2), show.legend = F, 
    alpha = 1, size = .2) + labs(color = feat1, color = "") + 
    paletteer::scale_color_paletteer_c("pals::kovesi.linear_bgy_10_95_c74") + 
    my_theme2 + coord_fixed() + 
    labs(title = glue("{feat2}"))
  
  p <- ggplot2::ggplot() + 
    ggplot2::geom_point(
      data = df4p, aes(x = x, y = y, color = feat1), show.legend = F, 
      alpha = .5, size = .2
    ) + labs(color = feat1) + 
    scale_color_continuous(low = "gray90", high = "#ff808f") +
    my_theme2 + coord_fixed() + labs(color = "tf") + 
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
    coord_fixed() + labs(title = "intersection of top 15%") + 
    theme(plot.title = element_text(hjust = .5))
  sp1 + sp2 + p2 + patchwork::plot_layout(ncol = 3, width = c(1, 1, 1))

})
pdf(glue("fig4gi_e115_tf_mz_co_expr.pdf"), width = 8, height = 3)
print(plst3)
dev.off()



