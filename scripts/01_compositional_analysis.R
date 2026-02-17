library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)

set.seed(123)

data_dir <- "data"   

samples_v <- paste0("Tumor", 1:20, "_singlets")

for (s in samples_v) {
  load(file.path(data_dir, paste0(s, ".RData")))
}

objs <- mget(samples_v, inherits = TRUE)


prep_for_merge <- function(o){
  md <- o@meta.data
  
  for (cn in colnames(md)) {
    x <- md[[cn]]
    if (is.matrix(x) && ncol(x) == 1) {
      md[[cn]] <- as.vector(x)
    }
  }
  o@meta.data <- md
  
  if (!"RNA" %in% Assays(o)) {
    old <- DefaultAssay(o)
    o <- RenameAssay(o, old, "RNA")
  }
  DefaultAssay(o) <- "RNA"
  
  if ("SCT" %in% Assays(o)) {
    o[["SCT"]] <- NULL
  }
  
  o <- DietSeurat(
    o,
    assays = "RNA",
    layers = c("counts", "data"),
    dimreducs = NULL,
    graphs = NULL,
    misc = FALSE
  )
  
  return(o)
}

objs <- lapply(objs, prep_for_merge)

NST_combined <- merge(
  x = objs[[1]],
  y = objs[-1],
  add.cell.ids = samples_v,
  project = "NST_combined"
)

NST_combined$Sample_ID <- gsub("^Tumor", "NST", NST_combined$orig.ident)
NST_combined$Sample_ID <- factor(
  NST_combined$Sample_ID,
  levels = paste0("NST", 1:20)
)

id_num <- suppressWarnings(
  as.integer(sub("^NST", "", as.character(NST_combined$Sample_ID)))
)

NST_combined$HER2_status <- case_when(
  id_num <= 6  ~ "0",
  id_num <= 14 ~ "1+",
  id_num <= 20 ~ "2+",
  TRUE ~ NA_character_
)

NST_combined$HER2_status <- factor(
  NST_combined$HER2_status,
  levels = c("0","1+","2+")
)

split_data <- SplitObject(NST_combined, split.by = "Sample_ID")

split_data <- lapply(split_data, function(x) {
  DefaultAssay(x) <- "RNA"
  SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

NST_combined <- merge(split_data[[1]], y = split_data[-1])

features <- SelectIntegrationFeatures(split_data, nfeatures = 3000)

DefaultAssay(NST_combined) <- "SCT"
VariableFeatures(NST_combined) <- features

NST_combined <- RunPCA(NST_combined, verbose = FALSE)
NST_combined <- RunUMAP(NST_combined, dims = 1:15, verbose = FALSE)
NST_combined <- FindNeighbors(NST_combined, dims = 1:15, verbose = FALSE)
NST_combined <- FindClusters(NST_combined, resolution = 0.5, verbose = FALSE)

NST_combined$seurat_clusters <- Idents(NST_combined)

NST_combined$Cell_type2 <- as.character(NST_combined$seurat_clusters)

NST_combined$Category <- case_when(
  grepl("^T-Epi", NST_combined$Cell_type2) ~ "T-Epi",
  NST_combined$Cell_type2 == "Basal/Myo" ~ "Basal/Myo",
  NST_combined$Cell_type2 %in% c("CAFs 1","CAFs 2","Pericytes","Endo","Adipo") ~ "Str cells",
  NST_combined$Cell_type2 %in% c("TAMs (SPP1-)","TAMs (SPP1+)",
                                 "T cells","T cells (Cycling)",
                                 "Plasma/B cells","DCs","Mast cells") ~ "Immune cells",
  TRUE ~ "Other"
)

meta <- NST_combined@meta.data %>%
  dplyr::select(Sample_ID, HER2_status, Category) %>%
  dplyr::filter(!is.na(Sample_ID),
                !is.na(HER2_status),
                !is.na(Category),
                Category != "Other")

tab_prop <- meta %>%
  count(Sample_ID, HER2_status, Category, name = "n") %>%
  group_by(Sample_ID) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

res_comp <- tab_prop %>%
  group_by(Category) %>%
  summarise(
    n_0  = n_distinct(Sample_ID[HER2_status == "0"]),
    n_1  = n_distinct(Sample_ID[HER2_status == "1+"]),
    n_2  = n_distinct(Sample_ID[HER2_status == "2+"]),
    p_kw = kruskal.test(prop ~ HER2_status)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_kw, method = "BH")) %>%
  arrange(p_adj)

write.csv(res_comp,
          "results_Resumen_Kruskal_por_Category.csv",
          row.names = FALSE)

pair_comparisons <- list(
  c("0","1+"),
  c("0","2+"),
  c("1+","2+")
)

p_all <- ggplot(tab_prop,
                aes(x = HER2_status, y = prop, fill = HER2_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(aes(color = HER2_status),
              width = 0.15, size = 1.8, alpha = 0.9) +
  facet_wrap(~ Category, scales = "fixed", ncol = 4) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = c(
    "0"  = "#648FFF",
    "1+" = "#FFB000",
    "2+" = "#DC267F"
  )) +
  scale_color_manual(values = c(
    "0"  = "#355E9E",
    "1+" = "#B37D00",
    "2+" = "#A61D63"
  )) +
  labs(
    title = "Patient Category Proportions by HER2 Status",
    x = "HER2 status",
    y = "Proportion"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  stat_compare_means(method = "kruskal.test",
                     label = "p.format") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = pair_comparisons,
                     p.adjust.method = "BH",
                     label = "p.signif",
                     hide.ns = FALSE)

ggsave("results_Boxplots_Category_HER2_0_1_2.png",
       p_all,
       width = 12,
       height = 8,
       dpi = 600,
       bg = "white")
