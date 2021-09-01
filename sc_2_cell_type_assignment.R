library(Seurat)
library(patchwork)
library(ggplot2)
library(scales)
library(harmony)
library(viridis)
library(monocle3)
library(dplyr)
cols = c("gray","yellow","orange","red")
MOE = readRDS("data/MOE.rds")

#plot markers
marker_genes2 = c("Krt5","Ascl1","Neurod1","Gap43","Omp","Ascl3","Cyp2g1","Muc5b","S100b","Tmem212","C1qb","Ccr2","Igkc","Cd3g","Ngp")
plots = list()
for (i in 1:length(marker_genes2)) {
  plots[[i]] = FeaturePlot(MOE, features = marker_genes2[i], order = T, cols = cols)+NoLegend()+ggtitle(marker_genes2[i])
}
wrap_plots(plots, ncol = 5)


#clustering
EM = as.sparse(GetAssayData(MOE, assay = "RNA"))
gene_metadata = as.data.frame(rownames(EM))                  
rownames(gene_metadata) = gene_metadata$`rownames(EM)`
colnames(gene_metadata) = "gene_short_name"
cds <- new_cell_data_set(as(EM, "sparseMatrix"),
                         cell_metadata = MOE@meta.data, gene_metadata = gene_metadata)
cds@int_colData@listData[["reducedDims"]][["UMAP"]] <-MOE@reductions[["umap"]]@cell.embeddings
cds = cluster_cells(cds,reduction_method = "UMAP",k = 5)
MOE$UMAP_cluster =clusters(cds, reduction_method = "UMAP")

DimPlot(MOE, label = T, group.by = "UMAP_cluster")+NoLegend()

EM = GetAssayData(MOE, slot = "counts")
EM = EM[grep("Olfr", rownames(EM)),]
EM[EM>0] = 1
ORnum = colSums(as.matrix(EM))
MOE$ORnum = ORnum
OR_multi = ORnum>1
MOE$OR_multi = OR_multi
FeaturePlot(MOE, features = "ORnum", order = T)
DimPlot(MOE, group.by = "OR_multi")

#cell type assigment
metadf = MOE@meta.data
metadf$cell_type = "OSN"
metadf[metadf$UMAP_cluster%in%c(63,33,65,72),]$cell_type = "HBC"
metadf[metadf$UMAP_cluster%in%c(85,66,86),]$cell_type = "GBC"
metadf[metadf$UMAP_cluster%in%c(78,69),]$cell_type = "INP"
metadf[metadf$UMAP_cluster%in%c(51,82,11,56,53,38,58),]$cell_type = "iOSN"
metadf[metadf$UMAP_cluster%in%c(54,30,77),]$cell_type = "MV"
metadf[metadf$UMAP_cluster%in%c(20,31,18,17,15,3,7,36),]$cell_type = "SUS"
metadf[metadf$UMAP_cluster%in%c(87,60),]$cell_type = "BG"
metadf[metadf$UMAP_cluster%in%c(75,90),]$cell_type = "OEC"
metadf[metadf$UMAP_cluster%in%c(62),]$cell_type = "EC"
metadf[metadf$UMAP_cluster%in%c(74,49,81),]$cell_type = "MC"
metadf[metadf$UMAP_cluster%in%c(68),]$cell_type = "MP"
metadf[metadf$UMAP_cluster%in%c(71),]$cell_type = "BC"
metadf[metadf$UMAP_cluster%in%c(83),]$cell_type = "TC"
metadf[metadf$UMAP_cluster%in%c(67),]$cell_type = "NP"
metadf[metadf$UMAP_cluster%in%c(89),]$cell_type = "NKN"
MOE@meta.data = metadf
my_level = rev(c("HBC","GBC","INP","iOSN","OSN","MV","SUS","BG","OEC","EC","MC","MP","BC","TC","NP","NKN"))
MOE$cell_type <- factor(x = MOE$cell_type, levels = my_level)
DimPlot(MOE, group.by = "cell_type", label = T)

saveRDS(MOE, "MOE.rds")

Idents(MOE) = MOE$cell_type
MOE_markers = FindAllMarkers(MOE)
saveRDS(MOE_markers, "data/MOE_markers.rds")
#isolate OSN cluster
Idents(MOE) = MOE$cell_type
OSN = subset(MOE, ident = "OSN")
OSN = GetAssayData(OSN, "counts")
OSN = CreateSeuratObject(OSN)
OSN = NormalizeData(OSN)
OSN = FindVariableFeatures(OSN)
OSN = ScaleData(OSN)
OSN = RunPCA(OSN)
ElbowPlot(OSN, ndims = 50)
OSN = RunHarmony(OSN, group.by.vars = "orig.ident")
OSN = RunUMAP(OSN, reduction = "harmony", dims = 1:30)
DimPlot(OSN, group.by = "orig.ident", split.by = "orig.ident", ncol = 3)
saveRDS(OSN, "data/OSN.rds")


