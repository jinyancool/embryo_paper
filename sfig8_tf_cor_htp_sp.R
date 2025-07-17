pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "Seurat", 
          "viridis")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/sfig8")
workdir |> fs::dir_create() |> setwd()

yaml_fn <- "~/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn= yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn = yaml_fn, "stage")[c("CS12", "CS14", "CS18")]

my_theme1 <- theme_classic(base_size = 8) + 
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(color = "black"), 
        axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"))

## figS8a-c: mouse TF clustering -----
samples <- c("ME9.5", "ME11.5x1", "ME13.5")
for(samp in samples) {
  csv_fn2 <- 
    glue("/cluster/home/ztao_jh/projects/embryo/analysis/zhangjing/human/rnaseq/pyscenic/mouse/{samp}/{samp}_SCENIC.csv")
  tf_score <- read_csv(csv_fn2) |> mutate(barcode = str_sub(Cell, end = -3))
  
  tf_sds <- tf_score |> as.data.frame() |> dplyr::select(-c("barcode")) |> 
    column_to_rownames("Cell") |> t() |> matrixStats::rowSds()
  cor_mtx <- tf_score |> as.data.frame() |> dplyr::select(-c("barcode")) |> 
    column_to_rownames("Cell") |> as.matrix() |> .[, names(tf_sds[tf_sds > 0])] |> cor() 
  htp1 <- ComplexHeatmap::Heatmap(cor_mtx, show_row_names = F, show_column_names = F, 
                                  clustering_method_rows = "complete", clustering_method_columns = "complete", 
                                  row_split = 7, column_split = 7, name = "cor", raster_by_magick = T)
  pdf(glue("{samp}_cor_htp1.pdf"), width = 5, height = 5)
  print(htp1)
  dev.off()
}

## figS8d-f ----
samples <- c("yao1", "yao2", "yao5")
for(samp in samples) {
  csv_fn2 <- 
    glue("/cluster/home/ztao_jh/projects/embryo/analysis/zhangjing/human/rnaseq/pyscenic/human/{samp}/{samp}_SCENIC.csv")
  tf_score <- read_csv(csv_fn2) |> mutate(barcode = str_sub(Cell, end = -3))
  
  tf_sds <- tf_score |> as.data.frame() |> dplyr::select(-c("barcode")) |> 
    column_to_rownames("Cell") |> t() |> matrixStats::rowSds()
  cor_mtx <- tf_score |> as.data.frame() |> dplyr::select(-c("barcode")) |> 
    column_to_rownames("Cell") |> as.matrix() |> .[, names(tf_sds[tf_sds > 0])] |> cor() 
  htp1 <- ComplexHeatmap::Heatmap(cor_mtx, show_row_names = F, show_column_names = F, 
                                  clustering_method_rows = "complete", clustering_method_columns = "complete", 
                          row_split = 7, column_split = 7, name = "cor", raster_by_magick = T)
  pdf(glue("{samp}_cor_htp1.pdf"), width = 5, height = 5)
  print(htp1)
  dev.off()
  
}
