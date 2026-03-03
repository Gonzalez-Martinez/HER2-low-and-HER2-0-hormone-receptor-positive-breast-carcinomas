suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(future)
  library(clustree)
  library(openxlsx)
})

set.seed(123)
options(future.globals.maxSize = 64 * 1024^3)

stopifnot(exists("NST_combined"))

#GENERAL

Idents(NST_combined) <- "Cell_type"

idents_keep <- grep("^(T cells|TAMs|B|Plasma|Mast|DCs|Osteoclast)", levels(NST_combined), value = TRUE)
Immune_subset <- subset(NST_combined, idents = idents_keep)

DefaultAssay(Immune_subset) <- "SCT"

Immune_subset <- RunPCA(Immune_subset, assay = "SCT", verbose = FALSE)
Immune_subset <- RunUMAP(Immune_subset, dims = 1:13, verbose = FALSE)
Immune_subset <- FindNeighbors(Immune_subset, dims = 1:13, verbose = FALSE)
Immune_subset <- FindClusters(
  Immune_subset,
  resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4),
  verbose = FALSE
)

clustree(Immune_subset, prefix = "SCT_snn_res.")

res_immune <- "SCT_snn_res.0.2"
Idents(Immune_subset) <- res_immune
Immune_subset$seurat_clusters <- Immune_subset[[res_immune, drop = TRUE]]

Immune_subset <- PrepSCTFindMarkers(Immune_subset, assay = "SCT", verbose = FALSE)
markers_Immune_res02 <- FindAllMarkers(Immune_subset, assay = "SCT", only.pos = TRUE, recorrect_umi = FALSE)
save(markers_Immune_res02, file = "results/markers_Immune_res02.RData")

top100_Immune_res02 <- markers_Immune_res02 %>%
  group_by(cluster) %>%
  slice_head(n = 100)

write.xlsx(top100_Immune_res02, file = "results/top100_markers_Immune_res02.xlsx", rowNames = FALSE)

map_ids_immune <- c(
  "0" = "T cells",
  "1" = "TAMs (SPP1-)",
  "2" = "TAMs (SPP1+)",
  "3" = "Plasma cells",
  "4" = "Mast cells",
  "5" = "B cells",
  "6" = "DCs (CD1A+)",
  "7" = "DCs (LAMP3+)",
  "8" = "Osteoclasts"
)

map_ids_immune <- map_ids_immune[names(map_ids_immune) %in% levels(Immune_subset)]
Immune_subset <- RenameIdents(Immune_subset, map_ids_immune)
Immune_subset$Cell_type <- as.character(Idents(Immune_subset))

tab_immune <- Immune_subset@meta.data %>%
  count(Sample_ID, HER2_group, Cell_type) %>%
  group_by(Sample_ID) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

res_comp_immune <- tab_immune %>%
  group_by(Cell_type) %>%
  summarise(
    n_0     = sum(HER2_group == "0"),
    n_low   = sum(HER2_group == "Low"),
    p_wilcox = tryCatch(wilcox.test(prop ~ HER2_group)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_wilcox, method = "BH")) %>%
  arrange(p_adj)

write.xlsx(res_comp_immune, file = "results/wilcox_Immune_res02.xlsx", rowNames = FALSE)

cv_immune <- tab_immune %>%
  group_by(Cell_type) %>%
  summarise(mean_prop = mean(prop), sd_prop = sd(prop), CV = sd_prop / mean_prop, .groups = "drop")

p_cv_immune <- ggplot(cv_immune, aes(x = reorder(Cell_type, CV), y = CV, fill = CV)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = "Inter-tumoral variability of immune populations",
    x = "Immune population",
    y = "Coefficient of variation (CV)"
  ) +
  theme(legend.position = "none")

# T CELLS

Tcells <- subset(Immune_subset, subset = Cell_type == "T cells")
DefaultAssay(Tcells) <- "SCT"

Tcells@meta.data <- Tcells@meta.data[, !grepl("^SCT_snn_res", colnames(Tcells@meta.data)), drop = FALSE]

Tcells <- RunPCA(Tcells, assay = "SCT", verbose = FALSE)
Tcells <- FindNeighbors(Tcells, dims = 1:11, verbose = FALSE)
Tcells <- FindClusters(Tcells, resolution = c(0.2,0.3,0.4), verbose = FALSE)
Tcells <- RunUMAP(Tcells, dims = 1:11, verbose = FALSE)

res_t <- "SCT_snn_res.0.3"
Idents(Tcells) <- res_t
Tcells$seurat_clusters <- Tcells[[res_t, drop = TRUE]]

Tcells <- PrepSCTFindMarkers(Tcells, assay = "SCT", verbose = FALSE)
markers_Tcells_res03 <- FindAllMarkers(Tcells, only.pos = TRUE, assay = "SCT", recorrect_umi = FALSE)
save(markers_Tcells_res03, file = "results/markers_Tcells_res03.RData")

top100_Tcells_res03 <- markers_Tcells_res03 %>%
  group_by(cluster) %>%
  slice_head(n = 100)

write.xlsx(top100_Tcells_res03, file = "results/top100_markers_Tcells_res03.xlsx", rowNames = FALSE)

annot_t <- c(
  "0" = "CD4+ Treg",
  "1" = "CD4+ Naive",
  "2" = "CD8+ Cytotoxic",
  "3" = "CD8+ TRM",
  "4" = "NKT",
  "5" = "ISG-high"
)

annot_t <- annot_t[names(annot_t) %in% levels(Tcells)]
Tcells <- RenameIdents(Tcells, annot_t)
Tcells$Cell_type_T <- as.character(Idents(Tcells))

tab_t <- Tcells@meta.data %>%
  count(Sample_ID, HER2_group, Cell_type_T) %>%
  group_by(Sample_ID) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(HER2_group = factor(HER2_group, levels = c("0","Low")))

res_comp_t <- tab_t %>%
  group_by(Cell_type_T) %>%
  summarise(
    n_0     = sum(HER2_group == "0"),
    n_low   = sum(HER2_group == "Low"),
    p_wilcox = tryCatch(wilcox.test(prop ~ HER2_group)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_wilcox, method = "BH")) %>%
  arrange(p_adj)

write.xlsx(res_comp_t, file = "results/wilcox_Tcells_res03.xlsx", rowNames = FALSE)

p_t_box <- ggplot(tab_t, aes(x = HER2_group, y = prop, fill = HER2_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.65, colour = "grey20") +
  geom_jitter(aes(color = HER2_group), width = 0.15, height = 0, size = 1.8, alpha = 0.9) +
  facet_wrap(~ Cell_type_T, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Patient-level proportions of T cell populations",
    x = "HER2 status",
    y = "Proportion"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text       = element_text(face = "bold", size = 12),
    plot.title       = element_text(face = "bold", hjust = 0, size = 14)
  ) +
  coord_cartesian(clip = "off") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("0","Low")),
    label = "p.signif",
    hide.ns = FALSE,
    size = 4,
    fontface = "bold",
    label.y.npc = 0.85,
    tip.length = 0.003,
    bracket.size = 0.4
  )


# MYELOIDS

Myeloid_subset <- subset(
  Immune_subset,
  subset = Cell_type %in% c(
    "TAMs (SPP1-)",
    "TAMs (SPP1+)",
    "DCs (CD1A+)",
    "DCs (LAMP3+)",
    "Osteoclasts"
  )
)

DefaultAssay(Myeloid_subset) <- "SCT"

Myeloid_subset@meta.data <- Myeloid_subset@meta.data[
  , !grepl("^SCT_snn_res", colnames(Myeloid_subset@meta.data)),
  drop = FALSE
]

Myeloid_subset <- RunPCA(Myeloid_subset, assay = "SCT", verbose = FALSE)
Myeloid_subset <- FindNeighbors(Myeloid_subset, dims = 1:10, verbose = FALSE)
Myeloid_subset <- FindClusters(Myeloid_subset, resolution = c(0.2,0.3,0.4), verbose = FALSE)
Myeloid_subset <- RunUMAP(Myeloid_subset, dims = 1:10, verbose = FALSE)

res_my <- "SCT_snn_res.0.2"
Idents(Myeloid_subset) <- res_my
Myeloid_subset$seurat_clusters <- Myeloid_subset[[res_my, drop = TRUE]]

Myeloid_subset <- PrepSCTFindMarkers(Myeloid_subset, assay = "SCT", verbose = FALSE)
markers_Myeloid_res02 <- FindAllMarkers(Myeloid_subset, only.pos = TRUE, assay = "SCT", recorrect_umi = FALSE)
save(markers_Myeloid_res02, file = "results/markers_Myeloid_res02.RData")

top100_Myeloid_res02 <- markers_Myeloid_res02 %>%
  group_by(cluster) %>%
  slice_head(n = 100)

write.xlsx(top100_Myeloid_res02, file = "results/top100_markers_Myeloid_res02.xlsx", rowNames = FALSE)

annot_my <- c(
  "0" = "TAMs (IFN+)",
  "1" = "TAMs (C1Q+)",
  "2" = "TAMs (M2-like)",
  "3" = "cDCs",
  "4" = "DCs-Lang",
  "5" = "DCs (LAMP3+)",
  "6" = "Osteoclasts"
)

annot_my <- annot_my[names(annot_my) %in% levels(Myeloid_subset)]
Myeloid_subset <- RenameIdents(Myeloid_subset, annot_my)
Myeloid_subset$Cell_type_Myeloid <- as.character(Idents(Myeloid_subset))

tab_my <- Myeloid_subset@meta.data %>%
  count(Sample_ID, HER2_group, Cell_type_Myeloid) %>%
  group_by(Sample_ID) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(HER2_group = factor(HER2_group, levels = c("0","Low")))

res_comp_my <- tab_my %>%
  group_by(Cell_type_Myeloid) %>%
  summarise(
    n_0     = sum(HER2_group == "0"),
    n_low   = sum(HER2_group == "Low"),
    p_wilcox = tryCatch(wilcox.test(prop ~ HER2_group)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_wilcox, method = "BH")) %>%
  arrange(p_adj)

write.xlsx(res_comp_my, file = "results/wilcox_Myeloids_res02.xlsx", rowNames = FALSE)

p_my_box <- ggplot(tab_my, aes(x = HER2_group, y = prop, fill = HER2_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.65, colour = "grey20") +
  geom_jitter(aes(color = HER2_group), width = 0.15, height = 0, size = 1.8, alpha = 0.9) +
  facet_wrap(~ Cell_type_Myeloid, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Patient-level proportions of myeloid populations",
    x = "HER2 status",
    y = "Proportion"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text       = element_text(face = "bold", size = 12),
    plot.title       = element_text(face = "bold", hjust = 0, size = 14)
  ) +
  coord_cartesian(clip = "off") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("0","Low")),
    label = "p.signif",
    hide.ns = FALSE,
    size = 4,
    fontface = "bold",
    label.y.npc = 0.85,
    tip.length = 0.003,
    bracket.size = 0.4
  )

  
    suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(infercnv)
})

set.seed(123)

options(future.globals.maxSize = 64 * 1024^3)

stopifnot(exists("Epi_obj"))
stopifnot(exists("Tcells_obj"))

out_dir <- "results/infercnv_epi_vs_tcells"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(Epi_obj) <- "RNA"
DefaultAssay(Tcells_obj) <- "RNA"

epi_counts <- GetAssayData(Epi_obj, assay = "RNA", layer = "counts")
t_counts   <- GetAssayData(Tcells_obj, assay = "RNA", layer = "counts")

common_genes <- intersect(rownames(epi_counts), rownames(t_counts))

epi_counts <- epi_counts[common_genes, , drop = FALSE]
t_counts   <- t_counts[common_genes, , drop = FALSE]

counts_mat <- Matrix::cBind(epi_counts, t_counts)

epi_cells <- colnames(epi_counts)
t_cells   <- colnames(t_counts)

cell_annot <- tibble(
  cell = c(epi_cells, t_cells),
  group = c(rep("Epithelial", length(epi_cells)), rep("T_cells", length(t_cells)))
)

annotation_file <- file.path(out_dir, "cell_annotations.tsv")
write.table(
  cell_annot,
  file = annotation_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

gene_order_file <- "resources/gencode_gene_pos_hg38.tsv"
stopifnot(file.exists(gene_order_file))

infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts_mat,
  annotations_file = annotation_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = c("T_cells")
)

infercnv_res_dir <- file.path(out_dir, "infercnv_out")
dir.create(infercnv_res_dir, recursive = TRUE, showWarnings = FALSE)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = infercnv_res_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE
)
