.libPaths(new = c(.libPaths(), "~/sbin/R/R-4.3.0"))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat",
          "paletteer", "cowplot", "ComplexHeatmap", "circlize", "parallel", "FELLA")
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/fig5")
workdir |> fs::dir_create() |> setwd()

yaml_fn <- "~/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn = yaml_fn, "stage")[c("E9.5", "E11.5", "E13.5")]

my_theme1 <- theme_classic(base_size = 8) +
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = .5))
my_theme2 <- theme_classic(base_size = 8) + theme(legend.key.size = unit(3, 'mm')) +
  theme(axis.line = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = .5),
        axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_rect(linewidth = .5, fill = NA))

## fig5b-f: m/z changes of somite among the trends from anterior to posterior ------
rds_fn1 <- "~/projects/collabrators/analysis/wangwenjie/align_new/mtb_new/pseudotime/cds_vld.rds"
cds_vld = read_rds(rds_fn1)
p3 <- monocle::plot_cell_trajectory(cds_vld, cell_size= .2, color_by = "Pseudotime") +
  theme_classic(base_size= 8) +
  theme( legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black")) +
  coord_fixed() + viridis::scale_color_viridis()
cds_vld$subarea = as.character(as.numeric(cds_vld$snn_res.0.6) + 1)
p4 <- monocle::plot_cell_trajectory(cds_vld, cell_size= .2, color_by = "subarea") +
  theme_classic(base_size= 8) +
  theme( legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black")) +
  coord_fixed()
ggsave("fig5b_ddrtree_dimplot.pdf", p4 / p3, width = 4, height = 3)

rds_fn1 <-
  "~/projects/collabrators/analysis/wangwenjie/align_new/mtb_new/pseudotime/seu_mz_somite.rds"
seu_somite = read_rds(rds_fn1)
seu_somite$subarea = seu_somite$snn_res.0.6 |> as.numeric() |> as.character() |> fct() |>
  fct_recode("1" = '5', "2" = "2", "3" = "4", "4" = "3", "5" = "1", "6" = "6") |> as.character()
write_rds(seu_somite, "../rds/fig5b_mouse_e11.5_mz_seu_somite.rds")
sp1 = Seurat::SpatialDimPlot(seu_somite, group.by = 'subarea', pt.size.factor = 2,
                             label = T, label.size = 2, label.color = "black", label.box = F,
                             repel = T) &
  theme_classic(base_size= 8) & labs(title = "The subarea of E11.5 somite") &
  theme( legend.key.size = unit(3, "mm")) & coord_fixed()
ggsave("fig5c_spatial_dimplot.pdf", sp1, width = 4, height = 3)
rds_fn2 <-
  "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig5d_mouse_e11.5_gene_seu1_somite2.rds"
seu1_somite2 <- read_rds(rds_fn2)
seu1_somite2$Pseudotime = cds_vld$Pseudotime

p1 = Seurat::DimPlot(seu_somite, reduction= "pca", group.by= "subarea", pt.size = .01,
                     label = T, label.size = 2, label.color = "black", label.box = F,
                     repel = T) +
  theme_classic(base_size= 8) + labs(title = "subarea: m/z") +
  theme( legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"))
ggsave("tst_p1.pdf", p1, width = 3, height = 3)
p2 = Seurat::DimPlot(seu1_somite2, reduction= "pca", group.by= "subarea",
                     pt.size = .01, label = T, label.size = 2, label.color = "black", label.box = F,
                     repel = T) +
  theme_classic(base_size= 8) + labs(title = "subarea: visium") +
  theme( legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"))
p11 <- Seurat::FeaturePlot(seu_somite, features = "Pseudotime", pt.size = .01, reduction = "pca") +
  theme_classic(base_size= 8) + labs(title = "Pseudotime: m/z") +
  theme( legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black")) +
  viridis::scale_color_viridis()
p21 <- Seurat::FeaturePlot(seu1_somite2, features = "Pseudotime", pt.size = .01, reduction = "pca") +
  theme_classic(base_size= 8) + labs(title = "Pseudotime: visium") +
  theme( legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black")) +
  viridis::scale_color_viridis()
p_comb = p2 + p1 + p21 + p11 + patchwork::plot_layout(nrow = 2)
ggsave("fig5d_pca_gene_mz.pdf", p_comb, width = 6, height = 6)

## fig5e: cluster of m/z trend -----
pseudotime <- cds$Pseudotime |> setNames(nm = colnames(cds)) |> sort()
mz_dat <- cds@assayData$exprs[, names(pseudotime)]

mz_idx <- rowSums(mz_dat) > 0
tmp <- mz_dat[mz_idx, ]

pt.matrix <- tmp
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

cm <- clusterData(exp = pt.matrix,
                  cluster.method = "mfuzz",
                  cluster.num = 6,
                  seed = 42)
sel_mz <- cm$long.res |> .$gene |> unique()
mz_info <- tibble(mode = str_sub(sel_mz, end = 3), mz = str_sub(sel_mz, start = 5)) |>
  mutate(mz_id = paste0(mode, "-", sprintf(as.numeric(mz), fmt = "%.5f"))) |>
  left_join(x = ., y = anot_all, by = c("mz_id" = "mz_id", "mode" = "mode", "mz" = "mz"))
write_tsv(mz_info, "trends_mz_info.tsv")

pdf("fig5e_pseudotime_mz_heatmap.pdf", width = 8, height = 10, onefile = FALSE)
visCluster(object = cm,
           plot.type = "both",
           show_row_dend = F,
           column_names_rot = 45,
           annnoblock.text = F,
           ctAnno.col = unname(cols_tissue[1:6]),
           sample.order = colnames(pt.matrix)
)
dev.off()

## fig5f: cluster1 of m/z trends enrichment -----
xlsx_fn5 <-
  "~/projects/collabrators/analysis/wangwenjie/align_new/mtb_new/pseudotime_250116/mz_6trends/fella_enrich_up_res_moran_cluster 1.xlsx"
enrich_res <- readxl::read_xlsx(xlsx_fn5, sheet = 3)

bar1 <- enrich_res |> head(n = 10) |> .[-c(3, 4, 6, 9), ] |>
  dplyr::arrange(p.value) |> mutate(KEGG.name = fct(KEGG.name)) |>
  ggplot(aes(x = -log10(p.value), y = KEGG.name)) +
  geom_bar(fill = cols_tissue[1], stat = "identity") +
  my_theme1 + ggplot2::scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 30)) +
  labs(title = "Cluster 1", y = "") +
  theme(legend.position = "none")
ggsave("fig5f_mz_cluster1_enrichment.pdf", bar1, width = 8, height = 4, unit = "cm")

## fig5h: spatial co-occurrence of signaling pathways and metabolic pathways -----
rds_fn8 <- "~/projects/collabrators/analysis/wangwenjie/mouse/figures/rds/fig4gi_obj_lst_v2.rds"
obj_lst <- read_rds(rds_fn8)
obj1 <- obj_lst[["E115"]]
cord_df = Seurat::GetTissueCoordinates(obj1) |> .[colnames(obj1), ] |> tibble() |>
  mutate(x = imagecol, y = -1 * imagerow)
cord1 <- cbind(cord_df, obj1@meta.data)
comp_df <- tibble(signaling = c("Fgf-Erk signaling pathway"),
                  metabolic = c("Glycolysis / Gluconeogenesis", "Purine metabolism",
                                "Pyrimidine metabolism"))
plst3 <- lapply(1:nrow(comp_df), \(idx) {
  feat1 <- comp_df[["signaling"]][idx] |> as.character() #|> gsub(" / ", "_", x = .)
  feat2 <- comp_df[["metabolic"]][idx] |> as.character() #|> gsub(" / ", "_", x = .)
  df4p <- cord1 |> dplyr::select(all_of(c("x", "y", feat1, feat2))) |>
    dplyr::rename("feat1" = feat1, "feat2" = feat2) |>
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
      data = df4p |> dplyr::filter(top == "top 15%"),
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
      data = df4p |> dplyr::filter(top == "top 15%"),
      contour_var	= "ndensity", alpha = .8,
      show.legend = T, linewidth = .2
    ) +
    coord_fixed() + labs(title = glue("{feat1}\n&\n{feat2}")) +
    theme(plot.title = element_text(hjust = .5))
 return(p2)

}) |> patchwork::wrap_plots(ncol = 1)
pdf(glue("fig5h_e115_signal_mtb_path_co_expr.pdf"), width = 3, height = 10)
  print(plst3)
dev.off()



