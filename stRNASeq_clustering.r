library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Hmisc)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(reshape2)

a1 <- Load10X_Spatial(data.dir = 'd://data_backup/sc_embryo/sp_a/', filename = 'filtered_feature_bc_matrix.h5', assay = "Spatial", slice = 'a1')
a1$orig.ident <- 'a1'

a1 <- PercentageFeatureSet(a1, pattern = "^mt-", col.name = "percent.mt")
a1 <- CellCycleScoring(object = a1, s.features = capitalize(tolower(cc.genes.updated.2019$s.genes)), g2m.features = capitalize(tolower(cc.genes.updated.2019$g2m.genes)))
a1 <- SCTransform(a1, assay = "Spatial", verbose = FALSE, return.only.var.genes = FALSE, vars.to.regress = c("percent.mt", 'G2M.Score', 'S.Score')) 

## E155 for example 
##### adding lables
lable1 <- read.csv('./155_type.csv', header = T, stringsAsFactors = F, row.names = 1)
lable2 <- read.csv('./155_site.csv', header = T, stringsAsFactors = F, row.names = 1)

cells <- intersect(rownames(lable1), rownames(lable2))
lable <- data.frame(row.names = cells, type = lable1[cells, 1], site = lable2[cells, 1], 
                        site_type = paste0(lable2[cells, 1], "_", lable1[cells, 1]))

sp_155 <- subset(a1, cells = cells)
sp_155$site <- factor(lable[rownames(sp_155@meta.data), 'site'], levels = lev)
sp_155$type <- factor(lable[rownames(sp_155@meta.data), 'type'], levels = c(0,1))
sp_155$site_type <- factor(lable[rownames(sp_155@meta.data), 'site_type'])

sp_155 <- FindSpatiallyVariableFeatures(sp_155, assay = "SCT", features = VariableFeatures(sp_155)[1:1000], 
                                       selection.method = "markvariogram")

sp_155 <- RunPCA(sp_155, assay = "SCT", verbose = FALSE)
sp_155 <- FindNeighbors(sp_155, reduction = "pca", dims = 1:30)
sp_155 <- FindClusters(sp_155, verbose = FALSE)
sp_155 <- RunUMAP(sp_155, reduction = "pca", dims = 1:30)

############# combine with scRNASeq data
DefaultAssay(sp_155) <- 'SCT'

obj_seurat <- readRDS('./union.rds')
anchors2 <- FindTransferAnchors(reference = sp_155, query = obj_seurat, normalization.method = "SCT")
predictions.assay2 <- TransferData(anchorset = anchors2, refdata = sp_155$site, prediction.assay = TRUE, weight.reduction = obj_seurat[["pca"]], dims = 1:28)
obj_seurat[["predictions"]] <- predictions.assay2

DefaultAssay(obj_seurat) <- "predictions"

obj_seurat$pred_score <- obj_seurat@assays$predictions@data['max',rownames(obj_seurat@meta.data)]
obj_seurat$site <- 'none'
tmp <- obj_seurat@assays$predictions@data[,rownames(obj_seurat@meta.data)]
for (i in 1:num){
  obj_seurat$site[i] <- names(which.max(tmp[,i]))
}