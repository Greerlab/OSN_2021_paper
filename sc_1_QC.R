# the mis-mapped counts of Olfr genes are adjusted with OR_deconvolution script (https://github.com/elisadonnard/ORdeconvolution) 
#ran the count.txt with Scrublet.ipynb
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(monocle3)
cols = c("gray","yellow","orange","red")
mtx = readRDS("data/raw_counts.rds")

#add the hb %, mt % and doublet score to the meta.data
mtx[["percent.hb"]] <- PercentageFeatureSet(mtx, pattern = "^Hb")
mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^mt-")
files = list.files(pattern = "db_score")
db_list = list()
for (i in 1:length(files)) {
  db_list[[i]] = read.table(paste0("",files[i]), header = T, row.names = 1)
}
db_df = do.call(rbind, db_list)
db_df = db_df[colnames(mtx),]
mtx$log_db = log1p(db_df$doublet_scores)
mtx$log_umi = log1p(mtx$nCount_RNA)

# plots for dignosis and cutoff
VlnPlot(mtx, features = c("nFeature_RNA", "nCount_RNA", "percent.hb", "percent.mt"), ncol = 4, pt.size = 0)
p1 = ggplot(mtx@meta.data, aes(x = nCount_RNA))+geom_histogram(binwidth =100)+xlim(c(500,10000))
p2 = ggplot(mtx@meta.data,aes(log_umi))+
  stat_bin(aes(y=cumsum(..count..)),geom="line",color="blue")+ylab("cumulative cell count")
p1+p2

#filter
clean_mtx = subset(mtx, subset = percent.hb < 1 
                   & percent.mt < 10 
                   & nFeature_RNA > 500 
                   & nCount_RNA >=1500
                   & nCount_RNA <=50000)
ncol(clean_mtx) #48162

# reanalyze the filtered matrix
clean_mtx = NormalizeData(clean_mtx)
clean_mtx = FindVariableFeatures(clean_mtx) 
clean_mtx = ScaleData(clean_mtx)
clean_mtx = RunPCA(clean_mtx)
  
ElbowPlot(clean_mtx, ndims = 50)
clean_mtx = RunUMAP(clean_mtx, dims = 1:30)

DimPlot(clean_mtx,split.by = "orig.ident",ncol = 3)
DimPlot(clean_mtx, label = T)+NoLegend()

p1 = FeaturePlot(clean_mtx, features = "log_db", cols = cols,order = T)
p2 = FeaturePlot(clean_mtx, features = "log_umi", cols = cols,order = T)
p1|p2

#check marker expression
marker_genes2 = c("Krt5","Ascl1","Neurod1","Gap43","Omp","Ascl3","Cyp2g1","Muc5b","S100b","Tmem212","C1qb","Ccr2","Igkc","Cd3g","Ngp")
plots = list()
for (i in 1:length(marker_genes2)) {
  plots[[i]] = FeaturePlot(clean_mtx, features = marker_genes2[i], order = T, cols = cols)+NoLegend()+ggtitle(marker_genes2[i])
}
wrap_plots(plots, ncol = 5)

#check OR expression
EM = GetAssayData(clean_mtx, slot = "counts")
EM = EM[grep("Olfr", value = T, rownames(EM)),]
EM[EM>0] = 1
OR_num = colSums(as.matrix(EM))
clean_mtx$OR_num = OR_num
FeaturePlot(clean_mtx, features = "OR_num", order = T)

#clustering
EM = as.sparse(GetAssayData(clean_mtx, assay = "RNA"))
gene_metadata = as.data.frame(rownames(EM))                  
rownames(gene_metadata) = gene_metadata$`rownames(EM)`
colnames(gene_metadata) = "gene_short_name"
cds <- new_cell_data_set(as(EM, "sparseMatrix"),
                         cell_metadata = clean_mtx@meta.data, gene_metadata = gene_metadata)
cds@int_colData@listData[["reducedDims"]][["UMAP"]] <-clean_mtx@reductions[["umap"]]@cell.embeddings
cds = cluster_cells(cds,reduction_method = "UMAP",k = 8)
clean_mtx$UMAP_cluster =clusters(cds, reduction_method = "UMAP")
DimPlot(clean_mtx, label = T, group.by = "UMAP_cluster")+NoLegend()

#remove the potrntial doublets with high dbscore and highly expressing multiple cell markers
cluster_to_remove = c(32,46,44,38,24,53)
Idents(clean_mtx) = clean_mtx$UMAP_cluster
clean_mtx2 = subset(clean_mtx, ident = cluster_to_remove, invert = T) #45317
clean_mtx2 = GetAssayData(clean_mtx2, slot = "counts")
clean_mtx2 = CreateSeuratObject(clean_mtx2)
clean_mtx2 = NormalizeData(clean_mtx2)
clean_mtx2 <- FindVariableFeatures(clean_mtx2)
clean_mtx2 <- ScaleData(clean_mtx2)
clean_mtx2 <- RunPCA(clean_mtx2)
ElbowPlot(clean_mtx2, ndims = 50)
clean_mtx2 = RunHarmony(clean_mtx2, group.by.vars = "orig.ident")
clean_mtx2 <- RunUMAP(clean_mtx2, reduction = "harmony", dims = 1:30)

saveRDS(clean_mtx2, "data/MOE.rds")









