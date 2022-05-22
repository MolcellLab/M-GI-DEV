library(Seurat)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(cowplot)
options(stringsAsFactors = FALSE)

## cellchat: st for example
st <- readRDS(object = st, file = './st.rds')

for (time in c('E9.5', 'E10.5', 'E11.5', 'E13.5','E15.5')){
  assign(x = paste0('st_', time, "_cc"), 
         value = cc_song(seurat = st, group = 'celltype_seg', t = time))
}


####################### function #############################
path_plot <- function(path, object, groupSize, receiver){
  for (pathways.show in path){
    ggsave(filename = paste0("./paths/", pathways.show, ".png"), dpi = 300, plot = netVisual_aggregate(object, signaling = pathways.show,  vertex.receiver = receiver, vertex.size = groupSize, pt.title = 15, title.space = 2), width = 16, height = 7)
    ggsave(filename = paste0("./paths/", pathways.show, "_circ.png"), dpi = 300, plot = netVisual_aggregate(object, signaling = pathways.show, layout = "circle", vertex.size = groupSize, pt.title = 15), width = 7, height = 6)
    ggsave(filename = paste0("./paths/", pathways.show, "_cont.png"), dpi = 300, plot = netAnalysis_contribution(object, signaling = pathways.show), width = 8, height = 4)
    ggsave(filename = paste0("./paths/", pathways.show, "_role.png"), dpi = 300, plot = netVisual_signalingRole(object, signaling = pathways.show, width = 16, height = 5, font.size = 13.5, font.size.title = 15), width = 10, height = 6)
  }}

cc_fun <- function(seurat, group, t){
  dir.create(paste0('./', t))
  setwd(paste0('./', t))
  
  print('data prepare')
  obj <- subset(seurat, subset = time == t)
  data.input <- GetAssayData(obj, assay = "RNA", slot = "data") # normalized data matrix
  labels <- obj@meta.data[,group]
  identity <- data.frame(group = labels, row.names = names(labels))
  
  print('cellchat_s1')
  cc <- createCellChat(data = data.input)
  cc <- addMeta(cc, meta = identity, meta.name = "labels")
  cc <- setIdent(cc, ident.use = "labels")
  cc@DB <- CellChatDB.mouse 
  cc <- subsetData(cc) 
  future::plan("multiprocess", workers = 6) # do parallel
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- projectData(cc, PPI.mouse)
  
  print('cellchat_s2')
  ##Inference of cell-cell communication network 
  cc <- computeCommunProb(cc)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_signalingRole(cc, slot.name = "netP") 
  num <- round(length(table(labels))/2)
  vertex.receiver = seq(1,num) 
  signal_path = cc@netP$pathways
  groupsize <- as.numeric(table(cc@idents))
  dir.create('./paths')
  path_plot(object = cc, path = signal_path, groupSize = groupsize, receiver = vertex.receiver)
  
  print('cellchat_s3')
  nPatterns = 5 
  cc <- identifyCommunicationPatterns(cc, pattern = "outgoing", k = nPatterns)
  netAnalysis_river(cc, pattern = "outgoing")
  ggsave(filename = 'outgoing_river.png', width = 16, height = 8, dpi = 300)
  netAnalysis_dot(cc, pattern = "outgoing", font.size.title = 15, font.size = 12)
  ggsave(filename = 'outgoing_dot.png', width = 15, height = 5, dpi = 300)
  cc <- identifyCommunicationPatterns(cc, pattern = "incoming", k = nPatterns)
  netAnalysis_river(cc, pattern = "incoming")
  ggsave(filename = 'incoming_river.png', width = 16, height = 8, dpi = 300)
  netAnalysis_dot(cc, pattern = "incoming", font.size.title = 15, font.size = 12)
  ggsave(filename = 'incoming_dot.png', width = 15, height = 5, dpi = 300)
  
  print('cellchat_s4')
  cc <- computeNetSimilarity(cc, type = "functional", thresh = 0.25)
  cc <- netEmbedding(cc, type = "functional")
  cc <- netClustering(cc, type = "functional", k = 4)
  p1 <- netVisual_embedding(cc, type = "functional", pathway.remove.show = F, label.size = 3.5, title = "Functional classification of network similariity ")
  p2 <- netVisual_embeddingZoomIn(cc, type = "functional", nCol = 2)
  cc <- computeNetSimilarity(cc, type = "structural", thresh = 0.25)
  cc <- netEmbedding(cc, type = "structural")
  cc <- netClustering(cc, type = "structural")
  p3 <- netVisual_embedding(cc, type = "structural", label.size = 3.5, title = "Structural classification of network similariity")
  p4 <- netVisual_embeddingZoomIn(cc, type = "structural", nCol = 2)
  plot_grid(p1, p2)
  ggsave(filename = "netsim_functional.png", width = 20, height = 9, dpi = 300)
  plot_grid(p3, p4)
  ggsave(filename = "netsim_structural.png", width = 20, height = 9, dpi = 300)
  
  print('cellchat_s5')
  cc <- rankNetPairwise(cc)
  saveRDS(object = cc, paste0(t, "_cc.rds"))
  setwd('../')
  return(cc)
}