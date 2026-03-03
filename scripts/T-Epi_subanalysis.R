library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(future)
library(clustree)
library(openxlsx)

set.seed(123)
options(future.globals.maxSize = 8 * 1024^3)

Idents(NST_combined) <- "Cell_type"

idents_keep <- grep("^T-Epi", levels(NST_combined), value = TRUE)
obj <- subset(NST_combined, idents = idents_keep)

obj_list <- SplitObject(obj, split.by = "Sample_ID")

obj_list <- lapply(obj_list, function(x) {
  DefaultAssay(x) <- "RNA"
  SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
})

obj <- merge(obj_list[[1]], y = obj_list[-1])
DefaultAssay(obj) <- "SCT"

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
VariableFeatures(obj) <- features

obj <- RunPCA(obj, assay = "SCT")
obj <- RunUMAP(obj, dims = 1:15)
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, verbose = FALSE, resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4))

clustree(obj, prefix = "SCT_snn_res.")

res_id <- "SCT_snn_res.0.3"
Idents(obj) <- res_id
obj$seurat_clusters <- obj[[res_id, drop = TRUE]]

DefaultAssay(obj) <- "SCT"
options(future.globals.maxSize = 64 * 1024^3)
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)

mk <- FindAllMarkers(obj, only.pos = TRUE, assay = "SCT")
save(mk, file = "markers_obj_res03.RData")

top100 <- mk %>%
  group_by(cluster) %>%
  slice_head(n = 100)

write.xlsx(top100, file = "top100_markers_obj_res03.xlsx", rowNames = FALSE)

map_ids <- c(
  "0" = "T-Epi 1",
  "1" = "T-Epi 2",
  "2" = "T-Epi 3",
  "3" = "T-Epi 4 (Cycling)",
  "4" = "T-Epi 5",
  "5" = "T-Epi 6",
  "6" = "T-Epi 7",
  "7" = "T-Epi 8 (Cili)",
  "8" = "T-Epi 9",
  "9" = "T-Epi 10"
)

map_ids <- map_ids[names(map_ids) %in% levels(obj)]
obj <- RenameIdents(obj, map_ids)
obj$Cell_type <- as.character(Idents(obj))

md <- obj@meta.data %>%
  mutate(HER2_group = recode(as.character(HER2_group), "Null" = "0")) %>%
  filter(!is.na(Sample_ID), !is.na(HER2_group), !is.na(Cell_type)) %>%
  mutate(HER2_group = factor(HER2_group, levels = c("0", "Low")))

tab <- md %>%
  count(Sample_ID, HER2_group, Cell_type, name = "n") %>%
  group_by(Sample_ID) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  group_by(Cell_type) %>%
  filter(n_distinct(HER2_group) == 2) %>%
  ungroup()

res_comp <- tab %>%
  group_by(Cell_type) %>%
  summarise(
    n_0   = n_distinct(Sample_ID[HER2_group == "0"]),
    n_low = n_distinct(Sample_ID[HER2_group == "Low"]),
    p_wilcox = suppressWarnings(wilcox.test(prop ~ HER2_group)$p.value),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_wilcox, method = "BH")) %>%
  arrange(p_adj)

write.csv(res_comp, "results_Wilcox_CellType_HER2group.csv", row.names = FALSE)

p_all <- ggplot(tab, aes(x = HER2_group, y = prop, fill = HER2_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.65, colour = "grey20") +
  geom_jitter(aes(color = HER2_group), width = 0.15, height = 0, size = 1.8, alpha = 0.9) +
  facet_wrap(~ Cell_type, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("0" = "#648FFF", "Low" = "#DC267F")) +
  scale_color_manual(values = c("0" = "#355E9E", "Low" = "#A61D63")) +
  labs(title = "Patient-level epithelial population proportions", x = "HER2 group", y = "Proportion") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0, size = 14)
  ) +
  coord_cartesian(clip = "off") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("0","Low")),
    label = "p.signif",
    hide.ns = FALSE,
    size = 5,
    fontface = "bold",
    label.y.npc = 0.85,
    tip.length = 0.003,
    bracket.size = 0.4
  )

ggsave("results_Boxplots_Tepi_HER2_0_vs_Low.png", p_all, width = 12, height = 8, dpi = 600, bg = "white")

res_var <- tab %>%
  group_by(Cell_type) %>%
  summarise(mean_prop = mean(prop), sd_prop = sd(prop), CV = sd_prop / mean_prop, .groups = "drop")

p_cv <- ggplot(res_var, aes(x = reorder(Cell_type, CV), y = CV, fill = CV)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_gradient(low = "#C3D6E5", high = "#4D6B8A") +
  theme_classic(base_size = 20) +
  labs(title = "Inter-tumoral variability of epithelial populations", x = "Epithelial population", y = "Coefficient of variation (CV)") +
  theme(legend.position = "none")

ggsave("results_CV_Tepi.png", p_cv, width = 10, height = 6, dpi = 600, bg = "white")

save(obj, file = "Epi_object_res03.RData")

