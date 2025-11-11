#--- loading ---#
suppressPackageStartupMessages(library(scRNAseq))
sce.zeisel <- ZeiselBrainData()
suppressPackageStartupMessages(library(scater))
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel,
                                      id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))
#--- gene-annotation ---#
suppressPackageStartupMessages(library(org.Mm.eg.db))
rowData(sce.zeisel)$Ensembl <- mapIds(org.Mm.eg.db,
                                      keys=rownames(sce.zeisel), keytype="SYMBOL", column="ENSEMBL")
#> 'select()' returned 1:many mapping between keys and columns
#--- quality-control ---#
stats <- perCellQCMetrics(sce.zeisel, subsets=list(
  Mt=rowData(sce.zeisel)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent",
                                              "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]
#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.zeisel)
sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clusters)
sce.zeisel <- logNormCounts(sce.zeisel)
#--- variance-modelling ---#
dec.zeisel <- modelGeneVarWithSpikes(sce.zeisel, "ERCC")
top.hvgs <- getTopHVGs(dec.zeisel, prop=0.1)
library(scran)
top.zeisel <- getTopHVGs(dec.zeisel, n=2000)
set.seed(100)
sce.zeisel <- fixedPCA(sce.zeisel, subset.row=top.zeisel)

set.seed(123)
library(densvis)
dt <- densne(reducedDim(sce.zeisel, "PCA"), dens_frac = 0.5, dens_lambda = 0.5)
reducedDim(sce.zeisel, "dens-SNE") <- dt
dm <- densmap(reducedDim(sce.zeisel, "PCA"), dens_frac = 0.5, dens_lambda = 0.5)
reducedDim(sce.zeisel, "densMAP") <- dm
sce.zeisel <- runUMAP(sce.zeisel) # for comparison
sce.zeisel <- runTSNE(sce.zeisel) # for comparison

ts <- reducedDim(sce.zeisel, "TSNE")
tu <- reducedDim(sce.zeisel, "UMAP")
ds <- scale(dt)
ts <- scale(ts)
du <- scale(dm)
tu <- scale(tu)

# for (class in unique(sce.zeisel$level1class)) {
#   print(class)
#   vars_ds <- colVars(ds[sce.zeisel$level1class == class, ])
#   vars_ts <- colVars(ts[sce.zeisel$level1class == class, ])
#   print(vars_ds)
#   print(vars_ts)

#   vars_du <- colVars(du[sce.zeisel$level1class == class, ])
#   vars_tu <- colVars(tu[sce.zeisel$level1class == class, ])

#   print(vars_du)
#   print(vars_tu)
# }

vars_ds <- colVars(ds[sce.zeisel$level1class == "endothelial-mural", ])
vars_ts <- colVars(ts[sce.zeisel$level1class == "endothelial-mural", ])
print(vars_ds)
#> [1] 0.5667786 0.2703578
print(vars_ts)
#>     TSNE1     TSNE2 
#> 0.5781887 0.1228506


vars_du <- colVars(du[sce.zeisel$level1class == "endothelial-mural", ])
vars_tu <- colVars(tu[sce.zeisel$level1class == "endothelial-mural", ])

#> Error: mean(vars_du) > mean(vars_tu) is not TRUE
print(vars_du)
#> [1] 0.1654369 0.2237998
print(vars_tu)
#>     UMAP1     UMAP2 
#> 0.2775908 0.5723159
stopifnot(mean(vars_ds) > mean(vars_ts))
stopifnot(mean(vars_du) > mean(vars_tu))
