library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(clustree)
library(openxlsx)

set.seed(123)
options(future.globals.maxSize = 64 * 1024^3)

Idents(NST_combined_VF) <- "Cell_type2"

idents_keep <- c("CAFs 1", "CAFs 2")
obj <- subset(NST_combined_VF, idents = idents_keep)

DefaultAssay(obj) <- "SCT"
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)

mk_base <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE, recorrect_umi = FALSE)
save(mk_base, file = "markers_CAFs_base.RData")

top100_base <- mk_base %>%
  group_by(cluster) %>%
  slice_head(n = 100)

write.xlsx(top100_base, file = "top100_markers_CAFs_base.xlsx", rowNames = FALSE)

obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:14, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:14, verbose = FALSE)
obj <- FindClusters(
  obj,
  verbose = FALSE,
  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4)
)

clustree(obj, prefix = "SCT_snn_res.")

res_id <- "SCT_snn_res.0.1"
Idents(obj) <- res_id
obj$seurat_clusters <- obj[[res_id, drop = TRUE]]

DimPlot2(obj, group.by = "seurat_clusters", label = TRUE)
DimPlot2(obj, group.by = "Cell_type2", label = TRUE)
DimPlot2(obj, group.by = "seurat_clusters", label = TRUE, split.by = "Sample_ID")
DimPlot2(obj, group.by = "Sample_ID")

DefaultAssay(obj) <- "SCT"
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)

mk <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE, recorrect_umi = FALSE)
save(mk, file = "markers_CAFs_res01.RData")

top100 <- mk %>%
  group_by(cluster) %>%
  slice_head(n = 100)

write.xlsx(top100, file = "top100_markers_CAFs_res01.xlsx", rowNames = FALSE)

save(obj, file = "CAFs_object_res01.RData")

map_ids <- c(
  "0" = "CAFs 1",
  "1" = "CAFs 2",
  "2" = "CAFs 3"
)

map_ids <- map_ids[names(map_ids) %in% levels(obj)]
obj <- RenameIdents(obj, map_ids)
obj$Cell_type <- as.character(Idents(obj))

Idents(obj) <- "Cell_type"
