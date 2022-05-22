# **Mesenchymal-epithelial interaction regulates gastrointestinal tract development in mouse embryos**
This repository contains the scripts of the manuscript **Mesenchymal-epithelial interaction regulates gastrointestinal tract development in mouse embryos**.



## Abstract

After gut tube patterning in early embryos, the cellular and molecular changes of developing stomach and intestine remain largely unknown. Here, combining single-cell RNA-sequencing and spatial RNA-sequencing, we constructed a spatiotemporal transcriptomic landscape of the mouse stomach and intestine during embryonic day E9.5-E15.5. Several novel subpopulations were identified, including Lox+ stomach mesenchyme, Aldh1a3+ small intestinal mesenchyme and Adamdec1+ large intestinal mesenchyme. The regionalization and heterogeneity of both the epithelium and mesenchyme could be traced back to E9.5. The spatiotemporal distributions of cell clusters and the mesenchymal-epithelial interaction analysis indicate that a coordinated development of the epithelium and mesenchyme contribute to the stomach regionalization, intestine segmentation and villus formation. Using the gut tube-derived organoids, we found that the cell fate of the foregut and hindgut could be switched by the regional niche factors including FGFs and RA. Together, this work demonstrates the important function of the mesenchymal-epithelial interactions in early gastrointestinal tract development, laying a foundation for further dissection of the mechanisms governing this process.



## Content

- process_cellranger.sh : script of process of scRNASeq data using cellranger
- process_spaceranger.sh : script of process of spatial RNASeq data using spaceranger
- scRNASeq_clustering.r : script of union clustering of all 14 samples
- scRNASeq_monocle.r : script of trajectory inferencing analysis using monocle
- scRNASeq_cellchat.r : script of cell-cell interaction analysis using cellchat
- scRNASeq_cellphoneDB.sh : script of cell-cell interaction analysis using cellphoneDB
- stRNASeq_clustering.r : script of clustering and integration analysis of spatial RNASeq data
- Bulkseq_deseq2.r : script of bulk RNA-Seq data



## Data availability

The raw data have been deposited in the GEO under accession no. GSE186525.

