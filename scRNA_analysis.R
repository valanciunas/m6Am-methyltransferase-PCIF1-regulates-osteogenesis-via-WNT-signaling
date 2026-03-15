library(dplyr)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(DoubletFinder)
library(monocle)

# Load Doublet Finder ------
FindDoublets <- function(seurat.rna, PCs = 1:50, exp_rate = 0.02, sct = FALSE){  
  ## pK identification
  sweep.res.list <- paramSweep(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25,
                              pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                              sct = sct)
  
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25, 
                              pK = 0.09, nExp = nExp_poi.adj,
                              reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                              sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  #seurat.rna = subset(seurat.rna, Doublet_Singlet == 'Singlet')
  return(seurat.rna)
}


# Load the cellranger output files
Osteo29 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo29"), project = "Osteo29")
Osteo24 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo24"), project = "Osteo24")
Osteo32 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo32"), project = "Osteo32")
Osteo36 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo36"), project = "Osteo36")
Osteo43 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo43"), project = "Osteo43")
Osteo65 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo65"), project = "Osteo65")
Osteo63 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo63"), project = "Osteo63")
Osteo52 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo67"), project = "Osteo52")
Osteo71 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo71"), project = "Osteo71")
Osteo30 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo30"), project = "Osteo30")
Osteo73 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo73"), project = "Osteo73")
Osteo74 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo74"), project = "Osteo74")
Osteo77 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo77"), project = "Osteo77")
Osteo25 <- CreateSeuratObject(counts = Read10X(data.dir = "/work/zpan/osteo_data/practice/Osteo25"), project = "Osteo78")

# Assign sample identifiers to each Seurat object
Osteo29$sample <- "Osteo29"
Osteo24$sample <- "Osteo24"
Osteo32$sample <- "Osteo32"
Osteo36$sample <- "Osteo36"
Osteo43$sample <- "Osteo43"
Osteo65$sample <- "Osteo65"
Osteo63$sample <- "Osteo63"
Osteo52$sample <- "Osteo52"
Osteo71$sample <- "Osteo71"
Osteo30$sample <- "Osteo30"
Osteo73$sample <- "Osteo73"
Osteo74$sample <- "Osteo74"
Osteo77$sample <- "Osteo77"
Osteo25$sample <- "Osteo25"

# Preprocess samples individually to remove doublets -----
ob.list <- list(
  Osteo24, Osteo29, Osteo32, Osteo36, Osteo43, 
  Osteo52, Osteo65, Osteo63, Osteo71, Osteo30, Osteo73, Osteo74, 
  Osteo77, Osteo25)
for (i in 1:length(ob.list)) {
  ob.list[[i]][["percent.mt"]] <- PercentageFeatureSet(ob.list[[i]], pattern = "^MT-")
  ob.list[[i]] <- subset(ob.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 10)
  # ob.list[[i]] <- subset(ob.list[[i]], cells = cells.use)
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = VariableFeatures(object = ob.list[[i]]))
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:30, reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:50, reduction.name = "UMAP_dim50", reduction.key = "UMAP_dim50_")
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], dims = 1:30)
  ob.list[[i]] <- FindClusters(ob.list[[i]], resolution = 1)
  ob.list[[i]] <- FindDoublets(ob.list[[i]], PCs = 1:30, sct = FALSE, exp_rate = (length(colnames(ob.list[[i]]))/125000))
}

# Preprocess samples individually to remove doublets -----

Osteo24 <- ob.list[[1]]
Osteo29 <- ob.list[[2]]
Osteo32 <- ob.list[[3]]
Osteo36 <- ob.list[[4]]
Osteo43 <- ob.list[[5]]
Osteo52 <- ob.list[[6]]
Osteo65 <- ob.list[[7]]
Osteo63 <- ob.list[[8]]
Osteo71 <- ob.list[[9]]
Osteo30 <- ob.list[[10]]
Osteo73 <- ob.list[[11]]
Osteo74 <- ob.list[[12]]
Osteo77 <- ob.list[[13]]
Osteo25 <- ob.list[[14]]

ob.list_2 <- list(
  Osteo24, Osteo29, Osteo32, Osteo36, Osteo43, Osteo52, Osteo65, Osteo63, Osteo71, Osteo30, 
  Osteo73, Osteo74, Osteo77, Osteo25)
for (i in 1:length(ob.list_2)) {
  ob.list_2[[i]][["percent.mt"]] <- PercentageFeatureSet(ob.list_2[[i]], pattern = "^MT-")
  ob.list_2[[i]] <- subset(ob.list_2[[i]], subset = nFeature_RNA > 300 & nCount_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
  ob.list_2[[i]] <- subset(ob.list_2[[i]], subset = Doublet_Singlet == "Singlet")
}

# Reassign the quality-controlled and doublet-filtered Seurat objects from the list back to individual sample objects
Osteo24_2 <- ob.list_2[[1]]
Osteo29_2 <- ob.list_2[[2]]
Osteo32_2 <- ob.list_2[[3]]
Osteo36_2 <- ob.list_2[[4]]
Osteo43_2 <- ob.list_2[[5]]
Osteo52_2 <- ob.list_2[[6]]
Osteo65_2 <- ob.list_2[[7]]
Osteo63_2 <- ob.list_2[[8]]
Osteo71_2 <- ob.list_2[[9]]
Osteo30_2 <- ob.list_2[[10]]
Osteo73_2 <- ob.list_2[[11]]
Osteo74_2 <- ob.list_2[[12]]
Osteo77_2 <- ob.list_2[[13]]
Osteo25_2 <- ob.list_2[[14]]

Osteo24_2 <- ob.list_2[[1]]
Osteo29_2 <- ob.list_2[[2]]
Osteo32_2 <- ob.list_2[[3]]
Osteo36_2 <- ob.list_2[[4]]
Osteo43_2 <- ob.list_2[[5]]
Osteo52_2 <- ob.list_2[[6]]
Osteo65_2 <- ob.list_2[[7]]
Osteo63_2 <- ob.list_2[[8]]
Osteo71_2 <- ob.list_2[[9]]
Osteo30_2 <- ob.list_2[[10]]
Osteo73_2 <- ob.list_2[[11]]
Osteo74_2 <- ob.list_2[[12]]
Osteo77_2 <- ob.list_2[[13]]
Osteo25_2 <- ob.list_2[[14]]


# Create new seurat obejct after quality control
Osteo29 <- CreateSeuratObject(counts = Osteo29_3, project = "Osteo29_2")
Osteo24 <- CreateSeuratObject(counts = Osteo24_3, project = "Osteo24_2")
Osteo32 <- CreateSeuratObject(counts = Osteo32_3, project = "Osteo32_2")
Osteo36 <- CreateSeuratObject(counts = Osteo36_3, project = "Osteo36_2")
Osteo43 <- CreateSeuratObject(counts = Osteo43_3, project = "Osteo43_2")
Osteo65 <- CreateSeuratObject(counts = Osteo65_3, project = "Osteo65_2")
Osteo63 <- CreateSeuratObject(counts = Osteo63_3, project = "Osteo63_2")
Osteo52 <- CreateSeuratObject(counts = Osteo52_3, project = "Osteo52_2")
Osteo71 <- CreateSeuratObject(counts = Osteo71_3, project = "Osteo71_2")
Osteo30 <- CreateSeuratObject(counts = Osteo30_3, project = "Osteo30_2")
Osteo73 <- CreateSeuratObject(counts = Osteo73_3, project = "Osteo73_2")
Osteo74 <- CreateSeuratObject(counts = Osteo74_3, project = "Osteo74_2")
Osteo77 <- CreateSeuratObject(counts = Osteo77_3, project = "Osteo77_2")
Osteo25 <- CreateSeuratObject(counts = Osteo25_3, project = "Osteo25_2")

# Assign sample identifiers to each Seurat object
Osteo29$sample <- "Osteo29"
Osteo24$sample <- "Osteo24"
Osteo32$sample <- "Osteo32"
Osteo36$sample <- "Osteo36"
Osteo43$sample <- "Osteo43"
Osteo65$sample <- "Osteo65"
Osteo63$sample <- "Osteo63"
Osteo52$sample <- "Osteo52"
Osteo71$sample <- "Osteo71"
Osteo30$sample <- "Osteo30"
Osteo73$sample <- "Osteo73"
Osteo74$sample <- "Osteo74"
Osteo77$sample <- "Osteo77"
Osteo25$sample <- "Osteo25"


# Merge all seurat objects together 

Osteo14subj <- merge(Osteo29,
                     y = c(
                       Osteo24, Osteo32, Osteo36, Osteo43, Osteo65, Osteo63, Osteo52, 
                       Osteo71, Osteo30, Osteo73, Osteo74, Osteo77, Osteo25),
                     add.cell.ids = c(
                       "Osteo29", "Osteo24", "Osteo32", "Osteo36", "Osteo43", "Osteo65", "Osteo63", "Osteo52", 
                       "Osteo71", "Osteo30", "Osteo73", "Osteo74", "Osteo77", "Osteo25"),
                     project = "scOsteo"
)

# Perform normalization, variable feature selection, scaling, Harmony-based batch integration, clustering, and UMAP visualization for the 14 dataset
Osteo14subj <- NormalizeData(Osteo14subj, verbose=F)
Osteo14subj <- FindVariableFeatures(Osteo14subj)
Osteo14subj <- ScaleData(Osteo14subj)

Osteo14subj <- RunPCA(Osteo14subj, features = VariableFeatures(object = Osteo14subj))
Osteo14subj <- IntegrateLayers(
  object = Osteo14subj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
Osteo14subj <- JoinLayers(Osteo14subj)
Osteo14subj <- FindNeighbors(Osteo14subj, reduction = "harmony",dims = 1:30)
Osteo14subj <- FindClusters(Osteo14subj, resolution = 0.5)
head(Idents(Osteo14subj))
Osteo14subj <- RunUMAP(Osteo14subj, reduction = "harmony", dims = 1:50, reduction.name = "umap", verbose = FALSE)
DimPlot(Osteo14subj, reduction = "umap", label = TRUE) 

#find markers for every cluster compared to all remaining cells, report only the positive
Osteo14subj.markers <- FindAllMarkers(Osteo14subj, only.pos = TRUE)
Osteo14subj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj< 0.05)

#Annotate each cluster based on the marker gene expression
new.cluster.ids_2 <- c("Endothelial", "OB", "MSC_Fibro", "Neutrophil", "B", "RBC", "RBC", "Endothelial", "VSMCs", "Macrophage", 
                       "Macrophage", "NK/T", "MSC")
names(new.cluster.ids_2) <- levels(Osteo14subj)
Osteo14subj <- RenameIdents(Osteo14subj, new.cluster.ids_2)
DimPlot(Osteo14subj, reduction = "umap", label = TRUE) 
Osteo14subj$celltype<-Idents(Osteo14subj)

#visualizing gene expression
VlnPlot(Osteo14subj, feature= "PCIF1")
FeaturePlot(Osteo14subj, feature= "PCIF1")

saveRDS(Osteo14subj, "Osteo14subj.rds") 
