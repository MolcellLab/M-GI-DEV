setwd('c://workstudio/embryo/union_14/')

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(stringr)

fea_num = 2000
mt = 10
dim_num = 30

## e95 for example
tmp <- Read10X(data.dir = ".../e95/filtered_feature_bc_matrix/")
tmp@Dimnames[[2]] <- paste0(tmp@Dimnames[[2]], "-1")
e95 <- CreateSeuratObject(counts = tmp, project = "e95", min.cells = 3, min.features = 200)
e95 <- PercentageFeatureSet(e95, pattern = "^mt-", col.name = "percent.mt")
e95 <- subset(e95, subset = percent.mt < mt)

## merge samples at same time point 
e105 <- merge(x = e105_st, y = e105_i)
e115 <- merge(x = e115_st, y = merge(x = e115_si, y = e115_li))
e135 <- merge(x = e135_st, y = merge(x = e135_s1, y = merge(x = e135_s2, y = e135_l)))
e155 <- merge(x = e155_st, y = merge(x = e155_s1, y = merge(x = e155_s2, y = e155_l)))

tmp <- list(e95 = e95, e105 = e105, e115 = e115, e135 = e135, e155 = e155)
for (n in 1:length(tmp)){
  tmp[[n]] <- NormalizeData(tmp[[n]], verbose = T)
  tmp[[n]] <- FindVariableFeatures(tmp[[n]], selection.method = "vst", nfeatures = fea_num, verbose = FALSE)
}

## integtating
tmp.anchors <- FindIntegrationAnchors(object.list = tmp, dims = 1:dim_num, anchor.features = fea_num)
union <- IntegrateData(anchorset = tmp.anchors, dims = 1:dim_num)

DefaultAssay(union) <- "integrated"
union <- CellCycleScoring(object = union, s.features = capitalize(tolower(cc.genes.updated.2019$s.genes)), g2m.features = capitalize(tolower(cc.genes.updated.2019$g2m.genes)), set.ident = F)
union <- ScaleData(union, verbose = T, vars.to.regress = c("percent.mt", 'G2M.Score', 'S.Score'))
union <- RunPCA(union, npcs = dim_num, verbose = T)
union <- RunUMAP(union, reduction = "pca", dims = 1:dim_num)
union <- RunTSNE(union, reduction = "pca", dims = 1:dim_num)

DimPlot(union, reduction = "umap", group.by = "orig.ident", pt.size = 0.6)
VlnPlot(union, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
FeaturePlot(union, features = "percent.mt")

################ clustering
union <- FindNeighbors(union, dims = 1:dim_num)
union <- FindClusters(object = union, resolution = 2)