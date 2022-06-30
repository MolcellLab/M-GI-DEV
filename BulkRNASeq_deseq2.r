library(DESeq2)
library(ggplot2)
library(corrplot)
library("rgl")
library("RColorBrewer")
library(clusterProfiler)
library(enrichplot)
library(Hmisc)
library(pheatmap)
library(corrplot)
library(org.Mm.eg.db)
library(stringr)

setwd('opt/bulk_seq/deseq2/')

files <- dir('../data_v2/')
samples <- c()

for (i in files){
  df <- read.table(file = paste0('../data_v2/', i), sep = '\t', 
                   header = F, stringsAsFactors = F, row.names = 1)
  name = strsplit(x = i, split = 'Reads')[[1]][1]
  samples <- c(samples,name)
  assign(x = name, value = df[-c(1:4),])
  
}
samples <- samples[c(2,3,1,4,9,12,10,11,5,7,6,8)]

exp_all <- data.frame(row.names = rownames(df)[-c(1:4)])
for (i in samples){
  df <- get(i)
  exp_all <- cbind(exp_all, df$V2[-c(1:4)])
}
colnames(exp_all) <- samples
coldata_all <- data.frame(row.names = colnames(exp_all), 
                          condition = factor(c('s1', 's1', 's2', 's2', 's3', 's3', 's4', 's4', 's5', 's5', 's6', 's6')))

dds_all = DESeqDataSetFromMatrix(countData = exp_all, colData = coldata_all, design = ~ condition)
keep = rowMeans(counts(dds_all))>10
dds_all=dds_all[keep,]
dds_all = DESeq(dds_all)
counts_normalized = as.data.frame(counts(dds_all, normalized = T))

#######################  diff between paires
id_symbol <- read.table(opt/resource/mm_id_name.txt', header = F, sep = '\t', stringsAsFactors = F, row.names = 1)

deseq_fun(ct1 = "Fore_P3_HST", ct2 = "Fore_P6_HST", ex1 = "F4_RAF4", ex2 = "Fore_P7_RAF4", 
           pvl = 0.01, lfc = 1, file = "comp1")
deseq_fun(ct1 = "H4_RAF4", ct2 = "Hind_P7_RAF4", ex1 = "Hind_P3_HST", ex2 = "Hind_P6_HST", 
           pvl = 0.01, lfc = 1, file = "comp2")

########### gsea
gsea_targ <- data.frame(set = c(rep('Foregut', 18), rep('hindgut', 37)), genes = g3)

res1 <- read.table('../comp1/results.txt', sep = '\t', stringsAsFactors = F)
gene_list1 = res1$stat
names(gene_list1) = res1$genes
gene_list1 = sort(gene_list1, decreasing = T)
gsea1 <- GSEA(geneList = gene_list1, TERM2GENE = gsea_targ)
gseaplot2(gsea1, geneSetID = 1:2)
ggsave(filename = 'comp1_gsea.png', dpi = 300, width = 5, height = 6)


########################## functions #############################
deseq_fun <- function(ct1, ct2, ex1, ex2, pvl, lfc, file, num){
  dir.create(paste0('./', file))
  setwd(paste0('./', file))
  ct1.exp <- get(x = ct1)
  ct2.exp <- get(x = ct2)
  ko1.exp <- get(x = ex1)
  ko2.exp <- get(x = ex2)
  ##
  counts <- data.frame(row.names = rownames(ct1.exp), CT_1 = ct1.exp$V2, CT_2 = ct2.exp$V2, 
                       KO_1 = ko1.exp$V2, KO_2 = ko2.exp$V2)
  coldata <- data.frame(row.names = colnames(counts), condition = factor(c('Control', 'Control', 'KO', 'KO')))
  ##
  dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
  dds = DESeq(dds)
  res <- na.omit(results(dds))
  ##
  counts_normalized = counts(dds, normalized = T)
  res_od =  res[order(res$padj),]
  counts_od <- counts_normalized[rownames(res_od),]
  results = as.data.frame(cbind(as.matrix(res_od), counts_od))
  results$genes <- id_symbol[rownames(results), 1]
  write.table(results, file = "results.txt", quote=F,sep="\t", col.names = T, row.names = T)
  ## plots
  vc_plot(result = results, fc = lfc, pv = pvl, file = 'vc', num = 30)
  hm_plot(result = results, num = 100, file = 'hm')
  go_fun(result = results, fc = lfc, pv = pvl, targ = 'mm_sig_gobp', file = 'gobp')
  go_fun(result = results, fc = lfc, pv = pvl, targ = 'mm_sig_kegg', file = 'kegg')
  setwd('../')
}
