pkgs <- c(
  "fs",
  "futile.logger",
  "configr",
  "stringr",
  "ggpubr",
  "ggthemes",
  "jhtools",
  "glue",
  "ggsci",
  "patchwork",
  "tidyverse",
  "dplyr",
  "Seurat",
  "viridis",
  "FELLA"
)
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "collabrators"
dataset <- "wangwenjie"
species <- "mouse"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/sfig5")
workdir |> fs::dir_create() |> setwd()

yaml_fn <- "~/projects/collabrators/code/wangwenjie/mouse/figures/configs.yaml"
cols_tissue <- jhtools::show_me_the_colors(config_fn = yaml_fn, "tissue")
stg_cols <- jhtools::show_me_the_colors(config_fn = yaml_fn, "stage")[c("CS12", "CS14", "CS18")]

my_theme1 <- theme_classic(base_size = 8) +
  theme(
    legend.key.size = unit(3, "mm"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

## figS5a: enrichment of each selected tissues in each stage ----
xlsx_fn1 =
  "/cluster/home/jhuang/projects/collabrators/data/wangwenjie/mouse/metabolism/mouse_adjusted/DZLM2023110146_DZLM2024030584-b1-张进-王文杰-空间代谢组-项目报告/2.定性结果/Qualitative.xlsx"
mz_anot_neg = readxl::read_excel(xlsx_fn1, sheet = "neg-all") |> mutate(mode = "neg")
mz_anot_pos = readxl::read_excel(xlsx_fn1, sheet = "pos-all") |> mutate(mode = "pos")
mz_anot = rbind(mz_anot_neg, mz_anot_pos) |> mutate(mz_id = paste0(mode, "-", mz))

csv_fn1 <- "../rds/fig2g_mouse_mz_liver_long_fc.csv"
mz_mmu = read_csv(csv_fn1) |> .[, 1:2] |> dplyr::distinct()

mz_trnd_lst = mz_mmu |> split(f = .$group)
tmpdir <- "~/ref/fella/kegg/mmu"
fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = c("diffusion", "pagerank", "hypergeom")
)


for (trnd in names(mz_trnd_lst)) {
  keg_id <- mz_anot |> dplyr::filter(mz_id %in% mz_trnd_lst[[trnd]][[1]]) |> .$KEGG |>
    na.omit() |> unique()
  anls <- enrich(
    compounds = keg_id,
    data = fella.data,
    method = c("diffusion", "pagerank", "hypergeom"),
    approx = "normality"
  )
  getExcluded(anls)
  enrich_res_lst = pbmcapply::pbmclapply(c("diffusion", "pagerank", "hypergeom"), \(mth) {
    res_tbl <- generateResultsTable(
      method = mth,
      threshold = 1,
      capPscores = 0,
      object = anls,
      data = fella.data
    )
    pth_nms <- vector(mode = "character", length = nrow(res_tbl))
    for (idx in 1:nrow(res_tbl)) {
      pth_nms[idx] <- fella.data@keggdata@id2name[[res_tbl$`KEGG.id`[idx]]]
    }
    res_tbl[["KEGG.name"]] <- stringr::str_split(pth_nms, " - Mus mus", simplify = T)[, 1]
    res_tbl <- res_tbl |> mutate(method = mth, trend = trnd)
  }, mc.cores = 3)
  writexl::write_xlsx(
    enrich_res_lst,
    path = glue("figS5f_enrich_{trnd}.xlsx"),
    col_names = T,
    format_headers = T
  )
}
# figS5f: enrichment of each selected tissues in each stage ----
kegg_mtb_pth = read_csv("~/ref/kegg/mmu")
enrich_res1 = lapply(names(mz_trnd_lst), \(trnd) {
  enrich_res = readxl::read_xlsx(glue("figS5f_enrich_{trnd}.xlsx"), sheet = 1)
}) |> bind_rows()
enrich_res2 = lapply(names(mz_trnd_lst), \(trnd) {
  enrich_res = readxl::read_xlsx(glue("figS5f_enrich_{trnd}.xlsx"), sheet = 2)
}) |> bind_rows()

enrich_res = lapply(names(mz_trnd_lst), \(trnd) {
  enrich_res = readxl::read_xlsx(glue("figS5f_enrich_{trnd}.xlsx"), sheet = 3)
}) |> bind_rows()
mmu_mtb_pth = read_rds("~/ref/kegg/mouse/mmu_mtb_pth_cpd.rds")

p_tbl1 = enrich_res |> dplyr::filter(KEGG.id %in% unique(mmu_mtb_pth[[1]])) |>
  group_by(trend) |> dplyr::slice_min(order_by = p.value, n = 4) |>
  dplyr::arrange(trend, desc(p.value)) |> dplyr::mutate(KEGG.name = as.character(KEGG.name) |> fct())
bar1 = ggplot2::ggplot(p_tbl1, aes(
  x = -log10(p.value),
  y = KEGG.name,
  fill = trend
)) +
  geom_bar(stat = "identity", position = "dodge") + my_theme1 +
  ggplot2::facet_grid(trend ~ ., space = "free", scales = "free_y") +
  scale_fill_manual(values = unname(cols_tissue)[c(1, 3, 4, 5, 6)]) + labs(y = "") +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = .5)
  )
ggsave(
  "mmu_liver_trend_mz_enrich_bar1.pdf",
  bar1,
  width = 8,
  height = 10,
  unit = "cm"
)
