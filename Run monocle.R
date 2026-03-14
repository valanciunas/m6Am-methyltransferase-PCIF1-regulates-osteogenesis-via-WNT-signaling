library(monocle)
# Run monocle
Osteo14subj <- readRDS("c:/CytoTRACE/Osteo14subj.RDS")
MSC<- subset(Osteo14subj, idents=c("MSC_Fibro", "MSC", "OB"))

scRNA<- MSC
ct <- scRNA@assays$RNA$counts
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
pd <- new("AnnotatedDataFrame",
          data=scRNA@meta.data)
sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds
sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)
fdif = "diff_test_res.Rdata"
if(!file.exists(fdif)){
  diff_test_res <- differentialGeneTest(sc_cds,
                                        fullModelFormulaStr = "~celltype",
                                        cores = 4)
  save(diff_test_res,file = fdif)
}
load(fdif)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
head(ordering_genes)
length(ordering_genes)
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
plot_ordering_genes(sc_cds)

sc_cds <- reduceDimension(sc_cds, max_components = 2, method = 'DDRTree')
sc_cds <- orderCells(sc_cds)
plot_cell_trajectory(sc_cds, color_by = "State")
plot_cell_trajectory(sc_cds, color_by = 'Pseudotime') 
plot_cell_trajectory(sc_cds, color_by = 'celltype')