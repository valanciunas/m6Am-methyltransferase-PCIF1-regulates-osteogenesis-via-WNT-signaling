library(readr)
library(CytoTRACE)
library(tibble)
library(Seurat)
# Run CytoTrace
Osteo14subj <- readRDS("c:/CytoTRACE/Osteo14subj.RDS")
MSC<- subset(Osteo14subj, idents=c("MSC_Fibro", "MSC", "OB"))

# Prepare CytoTRACE input
cti <- GetAssayData(MSC, slot = 'data', assay = 'RNA')
cti <- as.data.frame(cti) 
cytotrace <- CytoTRACE(cti) 

MSC <- AddMetaData(MSC, cytotrace$CytoTRACE, col.name = "CytoTRACE_score")
MSC$original.seurat.clusters <- Idents(MSC)

ctype.sele = subset(MSC@meta.data, select = c('original.seurat.clusters'))
phe.vec = as.character(ctype.sele$original.seurat.clusters)
names(phe.vec) = rownames(ctype.sele)

plotCytoTRACE(cytotrace,
              emb = MSC@reductions$umap@cell.embeddings)

plotCytoGenes(cytotrace)
