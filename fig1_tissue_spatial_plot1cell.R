pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "paletteer", "cowplot", "ComplexHeatmap", "circlize", "plot1cell")  
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/fig1")
workdir %>% fs::dir_create() %>% setwd()

yaml_fn <- "/cluster/home/danyang_jh/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn= yaml_fn, "stage")

my_theme1 <- theme_classic(base_size = 8) + 
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"), 
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = .5))

## figure1b: tissue atlas, with VISIUM data -----
rds_fn1 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mmu_visium_obj_lst.rds"
mmu_visium <- read_rds(rds_fn1)
samples_mmu <- c("E95", "E115", "E135")
### tissuetype in mouse -----
plst1 <- lapply(samples_mmu, \(samp) {
  obj1 = mmu_visium[[samp]]
  if(samp == "E95") {
    scale_fct <- 1.6
    plot_width <- 2
    plot_height <- 2
  } else if (samp == "E115") {
    scale_fct <- 1.1
    plot_width <- 4
    plot_height <- 4
  } else {
    scale_fct <- 2.0
    plot_width <- 6
    plot_height <- 6
  }
  dim1 = Seurat::SpatialDimPlot(obj1, pt.size.factor = scale_fct, image.alpha = 0, 
                         group.by = "tissuetype") & 
    ggplot2::scale_fill_manual(values = cols_tissue) & my_theme1 & 
    Seurat::NoAxes() & Seurat::NoLegend()
  ggsave(glue("fig1b_mouse_tissuetype_spatial_{samp}.pdf"), dim1, width = plot_width, height = plot_height)
  return(dim1)
}) %>% setNames(nm = samples_mmu)
### legend only ----
tbl1 = lapply(samples_mmu, \(samp) {
  tissuetypes <- mmu_visium[[samp]] %>% .$tissuetype %>% unique()
  tibble(tissuetype = tissuetypes)
}) %>% bind_rows()
pt1 <- ggplot(tbl1, aes(x = 1, y = tissuetype, color = tissuetype)) + 
  geom_point(size = .4) + ggplot2::scale_color_manual(values = cols_tissue) + 
  my_theme1 + coord_fixed() + labs(color = "")
legend <- cowplot::get_legend(pt1)
plst1[["legend"]] <- cowplot::plot_grid(NULL, legend, NULL, ncol = 3, rel_widths = c(.1, 1, .1))
ggsave("fig1b_legend_mouse_tissuetype_spatial.pdf", plst1[["legend"]], width = 2, height = 2)
### combination -----
fig1b_mouse <- cowplot::plot_grid(plotlist = c(plst1), nrow = 1, rel_widths = c(.8, 1, 1.2, 1))
ggsave("fig1b_mouse_tissuetype.pdf", fig1b_mouse, width = unit(8, "cm"), height = unit(3, "cm"))

## fig1c: human tissue atlas, VISIUM data -----
samples_hsa <- c("yao1", "yao2", "yao5")
rds_fn2 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/hsa_visium_obj_lst.rds"
hsa_visium <- read_rds(rds_fn2)
plst1 <- lapply(samples_hsa, \(samp) {
  obj1 = hsa_visium[[samp]]
  if(samp == "yao1") {
    scale_fct <- 1.8
    plot_width <- 2
    plot_height <- 2
  } else if (samp == "yao2") {
    scale_fct <- 1.3
    plot_width <- 4
    plot_height <- 4
  } else {
    scale_fct <- 1.5
    plot_width <- 6
    plot_height <- 6
  }
  dim1 = Seurat::SpatialDimPlot(obj1, pt.size.factor = scale_fct, image.alpha = 0, 
                                group.by = "tissue") & 
    ggplot2::scale_fill_manual(values = cols_tissue) & my_theme1 & 
    Seurat::NoAxes() & Seurat::NoLegend()
  ggsave(glue("fig1c_human_tissue_spatial_{samp}.pdf"), dim1, width = plot_width, height = plot_height)
  return(dim1)
}) %>% setNames(nm = samples_hsa)
### legend only ----
tbl1 = lapply(samples_hsa, \(samp) {
  tissues <- hsa_visium[[samp]] %>% .$tissue %>% unique()
  tibble(tissue = tissues)
}) %>% bind_rows() %>% dplyr::distinct()
pt1 <- ggplot(tbl1, aes(x = 1, y = tissue, color = tissue)) + 
  geom_point(size = .4) + ggplot2::scale_color_manual(values = cols_tissue) + 
  my_theme1 + coord_fixed() + labs(color = "")
legend <- cowplot::get_legend(pt1)
plst1[["legend"]] <- cowplot::plot_grid(NULL, legend, NULL, ncol = 3, rel_widths = c(.1, 1, .1))
ggsave("fig1c_legend_human_tissue_spatial.pdf", plst1[["legend"]], width = 2, height = 2)
### combination -----
fig1c_human <- cowplot::plot_grid(plotlist = c(plst1), nrow = 1, rel_widths = c(.6, 1, 1.2, 1))
ggsave("fig1c_human_tissue.pdf", fig1c_human, width = unit(8, "cm"), height = unit(3, "cm"))


## fig1d: plot1cell -----
### mouse visium data, all genes of all samples -----
rds_fn4 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mmu_visium_merged_obj.rds"
visium_mmu_mrg <- read_rds(rds_fn4)
circ_dat <- prepare_circlize_data(visium_mmu_mrg, scale = .7)
clust_cols = cols_tissue[sort(unique(circ_dat$tissuetype))]
lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                    legend_gp = gpar(fill = stg_cols), 
                                    title_position = "topleft", title = "stage")
pdf(glue("fig1d_mouse_plot1cell_all_genes.pdf"), width = 8, height = 8)
plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
              col.use = clust_cols ,bg.color = 'white',
              kde2d.n = 200, repel = T, label.cex = .9)
add_track(circ_dat, group = "stage", colors = stg_cols[order(names(stg_cols))], 
          track_num = 2) 
draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
dev.off()
### mouse visium data, metabolic genes of all samples -----
rds_fn5 <-           
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mmu_visium_mtb_gene_merged_obj.rds"
seu_new3_mrg <- read_rds(rds_fn5)
circ_dat <- prepare_circlize_data(seu_new3_mrg, scale = .7)
clust_cols = cols_tissue[sort(unique(circ_dat$tissuetype))]
lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                    legend_gp = gpar(fill = stg_cols), 
                                    title_position = "topleft", title = "stage")
pdf(glue("fig1d_mouse_plot1cell_mtb_genes.pdf"), width = 8, height = 8)
plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
              col.use = clust_cols, bg.color = 'white',
              kde2d.n = 200, repel = T, label.cex = .9)
add_track(circ_dat, group = "stage", colors = stg_cols, track_num = 2) 
draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
dev.off()


### mouse m/z data, merged spots with visium coordinates -----
rds_fn2 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/mouse_mz_obj_merged.rds"
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

### human visium data, all genes of all samples -----
rds_fn6 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/human_visium_obj_merged.rds"
seu_new2_mrg <- read_rds(rds_fn6)
circ_dat <- prepare_circlize_data(seu_new2_mrg, scale = .7)
clust_cols = cols_tissue[sort(unique(circ_dat$tissue))]
lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                    legend_gp = gpar(fill = stg_cols), 
                                    title_position = "topleft", title = "stage")
pdf(glue("fig1d_human_plot1cell_all_genes.pdf"), width = 8, height = 8)
plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
              col.use = clust_cols ,bg.color = 'white',
              kde2d.n = 200, repel = T, label.cex = .9)

add_track(circ_dat, group = "stage", colors = stg_cols[sort(names(stg_cols))], track_num = 2) 

draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
dev.off()

### human visium data, metabolic genes of all samples -----
rds_fn8 <- "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/hsa_visium_mtb_gene_merged_obj.rds"
seu_new3_mrg <- read_rds(rds_fn8)
circ_dat <- prepare_circlize_data(seu_new3_mrg, scale = .7)
clust_cols = cols_tissue[sort(unique(circ_dat$tissue))]
lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                    legend_gp = gpar(fill = stg_cols), 
                                    title_position = "topleft", title = "stage")
pdf(glue("fig1d_human_plot1cell_mtb_genes_merged.pdf"), width = 8, height = 8)
plot_circlize(circ_dat, do.label = T, pt.size = 0.5,
              col.use = clust_cols ,bg.color = 'white',
              kde2d.n = 200, repel = T, label.cex = .9)
add_track(circ_dat, group = "stage", colors = stg_cols[sort(names(stg_cols))], track_num = 2) 
draw(lgd_squre, x = unit(40, "mm"), y = unit(12, "mm"), just = c("right", "bottom"))
dev.off()

### human m/z data, all merged samples -----
rds_fn4 = "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/human_mz_obj_merged.rds"
seu_mrg2 = read_rds(rds_fn4)
circ_dat <- prepare_circlize_data(seu_mrg2, scale = .7)

clust_cols = cols_tissue[sort(unique(circ_dat$tissue))]
lgd_squre <- ComplexHeatmap::Legend(at = names(stg_cols), type = "grid", 
                                    legend_gp = gpar(fill = stg_cols), 
                                    title_position = "topleft", title = "stage")
pdf(glue("fig1d_human_plot1cell_all_mz_mrg.pdf"), width = 5.5, height = 5.5)
plot_circlize(circ_dat, do.label = T, pt.size = 0.4, 
              col.use = clust_cols ,bg.color = 'white', contour.nlevels = 100, 
              kde2d.n = 2000, repel = T, label.cex = .7)

add_track(circ_dat, group = "stage", colors = stg_cols[sort(unique(circ_dat$stage))], track_num = 2) 

draw(lgd_squre, x = unit(25, "mm"), y = unit(8, "mm"), just = c("right", "bottom"))
dev.off()



