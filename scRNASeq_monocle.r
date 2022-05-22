library(monocle)
library(Seurat)
library(ggplot2)

## epithelial cells of stomach for example
st_epi <- readRDS('./st_epi.rds')

data <- as(as.matrix(st_epi@assays$RNA@counts), 'sparseMatrix')
pd <- st_epi@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds_st <- newCellDataSet(data,
                         phenoData = new("AnnotatedDataFrame", data = pd),
                         featureData = new("AnnotatedDataFrame", data = fData),
                         expressionFamily=negbinomial.size())

## Estimate size factors and dispersions
cds_st <- estimateSizeFactors(cds_st)
cds_st <- estimateDispersions(cds_st)

#s1: choose genes that define a cell's progress
disp_table <- dispersionTable(cds_st)
high_disp_genes <- subset(disp_table, dispersion_empirical >= dispersion_fit)
high_disp_genes_filter <- high_disp_genes[-grep("^mt- *|^Rp *", as.character(high_disp_genes$gene_id)),]
cds_st <- setOrderingFilter(cds_st, high_disp_genes_filter$gene_id)

#s2: reduce data dimensionality
cds_st <- reduceDimension(cds_st, method = 'DDRTree')
#s3: order cells along the trajectory
cds_st <- orderCells(cds_st)
#cds_st <- orderCells(cds_st, root_state = stem_state(cds_st))

tmp <- pData(cds_st)
st_epi$pseudotime <- tmp[rownames(st_epi@meta.data), 'Pseudotime']
FeaturePlot(object = st_epi, features = 'pseudotime', pt.size = .9, cols = c('red', 'blue'))
ggsave(filename = 'pseudotime_umap.png', dpi = 300, width = 9, height = 7)

plot_cell_trajectory(cds_st, color_by = "seurat_clusters", show_branch_points = F, cell_size = 1.5, cell_link_size = 1) #+ facet_wrap(facets = "seurat_clusters")

time_cols <- brewer.pal(n = 5, name = 'Spectral')
names(time_cols) <- names(table(pData(cds_st)$time))[1:5]
plot_cell_trajectory(cds_st, color_by = "time", show_branch_points = F, cell_size = 1.5, cell_link_size = 1) + scale_color_manual(values = time_cols, name = 'time')
ggsave(filename = 'traj_time.png', dpi = 300, width = 6, height = 4)


############################### diff exp genes ################################
## pseudotime only 
diff_pseudo <- differentialGeneTest(cds_st[high_disp_genes_filter$gene_id,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
## pseudotime + celltype
diff_pseudo_celltype <- differentialGeneTest(cds_st[high_disp_genes_filter$gene_id,],fullModelFormulaStr = "~seurat_clusters + Pseudotime",reducedModelFormulaStr = "~Pseudotime")

####### branched
BEAM_si <- BEAM(cds_st[high_disp_genes_filter$gene_id,], branch_point = 1, cores = 1)
BEAM_si <- BEAM_si[order(BEAM_si$qval),]

plot_genes_branched_heatmap(cds_st[BEAM_si$gene_short_name[1:100],],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)