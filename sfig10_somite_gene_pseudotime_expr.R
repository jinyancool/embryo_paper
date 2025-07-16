pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "viridis", "RColorBrewer")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/sfig10")
workdir %>% fs::dir_create() %>% setwd()

yaml_fn <- "/cluster/home/danyang_jh/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn = yaml_fn, "stage")[c("CS12", "CS14", "CS18")]

my_theme1 <- theme_classic(base_size = 8) + 
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"), 
        axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"))

## figS10d: monocle2 results of E11.5 somite, visium data -----
rds_fn1 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/align_new/mtb_new/pseudotime/cds_vld.rds"
cds_vld = read_rds(rds_fn1)

p1 <- monocle::plot_cell_trajectory(cds_vld, cell_size= .2, color_by = "State") + 
  theme_classic(base_size= 8) + scale_color_manual(values = unname(stg_cols)) + 
  theme( legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black")) + 
  coord_fixed()
ggsave("sfig10d_ddrtree_state.pdf", p1, width = 3, height = 2)

rds_fn2 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig5d_mouse_e11.5_gene_seu1_somite2.rds"
seu1_somite2 <- read_rds(rds_fn2)
seu1_somite2$State <- cds_vld$State

sp1 <- Seurat::SpatialDimPlot(seu1_somite2, group.by = 'State', pt.size.factor = 2,
                       label = F, label.size = 2, label.color = "black", label.box = F,
                       repel = T) & 
  theme_classic(base_size= 8) & #labs(title = "The State of E11.5 somite") &
  scale_fill_manual(values = unname(stg_cols)) & 
  theme(legend.key.size = unit(3, "mm")) & coord_fixed()
ggsave("sfig10d_e11.5_somite_state_spatial.pdf", sp1, width = 3, height = 2)

## figS10e: metabolic genes heatmap -----
rds_fn3 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/sfig10e_e11.5_somite_htp_df.rds"
df <- read_rds(rds_fn3)

pdf("sfig10e_e11.5_somite_heatmap1.pdf", width = 6, height = 8)
visCluster(object = df, plot.type = "both")
dev.off()

## figS10f: gene enrichment of metabolic gene cluster 4 ----
xlsx_fn1 <- 
  "/cluster/home/danyang_jh/projects/collabrators/analysis/wangwenjie/align_new/mtb_new/pseudotime_250116/mtb_gene_4clusters/enrich_kegg_res_lst.xlsx"
enrich_res <- readxl::read_excel(xlsx_fn1, sheet = 4)
enrich_res <- enrich_res %>% 
  dplyr::mutate(bg_num = str_split(BgRatio, "/", simplify = T)[, 1] %>% as.numeric()) %>% 
  dplyr::mutate(set_num = str_split(GeneRatio, "/", simplify = T)[, 2] %>% as.numeric()) %>% 
  mutate(ratio = Count/bg_num) %>% mutate(generatio = Count/set_num) %>% 
  dplyr::filter(category == "Metabolism") %>% 
  dplyr::slice_head(n = 10) %>% dplyr::arrange(generatio) %>% 
  mutate(Description = fct(as.character(Description)))
pt1 <- ggplot2::ggplot(enrich_res, aes(x = generatio, y = Description, color = -log10(p.adjust), size = Count)) + 
  geom_point() + viridis::scale_color_viridis() + my_theme1 + 
  ggplot2::scale_y_discrete(label = function(x) str_wrap(x, width = 30)) + 
  theme(axis.line = element_blank(), panel.border = element_rect(linewidth = .5, fill = NA, color = "black")) + 
  labs(x = "Gene Ratio", y = "", title = "Metabolic pathways in cluster4")
ggsave("sfig10f_kegg_mtb_enrich_clust4.pdf", pt1, width = 5, height = 3)

## figS10g-j: local focus of specific genes in the selected pathways -----
GetAllCoordinates <- function(.data) {
  .data@images %>%
    names() %>%
    unique() %>%
    map_dfr(~{
      GetTissueCoordinates(
        .data,
        image = .x,
        cols = c("row", "col"),
        scale = NULL
      ) %>%
        rownames_to_column(var = "cellid")
    })
}

celltype_isoheight_plot <- function(
    .data, 
    gn = NULL, 
    col_bg = "gray70",
    col_top = "darkred",
    col_isoheight = "white",
    col_white_ratio = .25,
    cols_fill_isoheight = c(
      rep("white", round(100 * .25)),
      colorRampPalette(brewer.pal(5, "YlOrRd")[2:5])(round(100 * (1 - .25)))
    ),
    size_bg = 0.2,
    size_top = size_bg
) {
  
  df <- .data@meta.data %>%
    rownames_to_column("cellid") %>%
    inner_join(
      GetAllCoordinates(.data)
    ) %>%
    as_tibble() %>% mutate(row = -1 * row)
  
  xy_r <- max(df$col, df$row)
  
  p <- ggplot(mapping = aes(x = row, y = col)) +
    ggplot2::geom_point(
      data = df, show.legend = T, 
      color = col_bg, alpha = .8, size = size_bg
    )
  p <- p +
    ggnewscale::new_scale_fill() +
    ggplot2::stat_density_2d_filled(
      data = df %>% dplyr::filter(top_n),
      mapping = aes(fill = ..ndensity.., alpha = ..ndensity.. ),
      geom = "raster", contour = F
    ) +
    scale_fill_gradientn(colours = cols_fill_isoheight) + 
    ggnewscale::new_scale_fill() +
    geom_density_2d(
      data = df %>% dplyr::filter(top_n),
      color = col_isoheight, 
      contour_var	= "ndensity", alpha = .5, 
      show.legend = T, linewidth = .2
    )
  
  p <- p +
    geom_point(
      data = df %>% dplyr::filter(top_n), show.legend = T, 
      color = col_top, alpha = .8, size = size_top
    ) + labs(color = 'top 10%')
  
  x_min <- min(df$row)
  x_max <- max(df$row)
  y_min = min(df$col)
  y_max = max(df$col)
  p <- p +
    ggplot2::coord_fixed(ratio = 1, xlim = c(1.05 * x_min, 0.9 * x_max), 
                         ylim = c(0.9 * y_min, 1.05 * y_max), expand = F) +
    theme_classic(base_size = 8) + 
    theme(
      legend.key.size = unit(3, "mm"), 
      axis.line = element_blank(), 
      axis.title = element_blank(), 
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      plot.title = element_text(hjust = .5), 
      panel.border = element_rect(linewidth = .5, fill = NA)
    ) + labs(title = gn)
}

sel_genes <- list(
  gly = c("Pgam2", "Pkm", "Eno3"), 
  aa = c("Shmt1", "Gatm"), 
  nucleotide = c("Paics", "Ampd1"), 
  fa = c("Hsd17b12", "Isyna1")
  )
feats_lst <- list()
for(gn in unlist(sel_genes)) {
  val <- LayerData(seu1_somite2[["Spatial"]], "counts")[gn, ]
  seu1_somite2$top_n <- case_when(val > quantile(val, probs = .85) ~ TRUE, TRUE ~ FALSE)
  feats_lst[[gn]] <- 
    celltype_isoheight_plot(.data = seu1_somite2, col_bg = "gray70", col_top = "darkred", 
                            col_isoheight = "white", gn = gn, col_white_ratio = .25,
                            cols_fill_isoheight = c(
                              rep("white", round(100 * .25)),
                              colorRampPalette(brewer.pal(5, "YlOrRd")[2:5])(round(100 * (1 - .25)))
                            ), size_bg = 0.2, size_top = .2)
}
pdf("sfig10gj_e11.5_somite_mtb_genes.pdf", width = 2, height = 3)
print(feats_lst)
dev.off()







