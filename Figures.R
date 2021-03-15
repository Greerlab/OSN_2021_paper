library(Seurat)
library(DESeq2)
library(BioVenn)
library(patchwork)
library(harmony)
library(dplyr)
library(ggplot2)
library(scales)
library(ggtext)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(rstatix)
library(EnhancedVolcano)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(GO.db)
library(Hmisc)
library(forcats)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(oce)
library(jpeg)
library(MLmetrics)
library(ggrepel)
library(reshape2)
library(spdep)
library(reticulate)
cols = c("gray","yellow","orange","red")
dds = readRDS("data/dds.rds") # DEseq2 objecut for bulk RNAseq
rds = readRDS("data/vsd.rds") # DEseq2 objecut for bulk RNAseq, transformed
MOE = readRDS("data/MOE.rds") # all cell types detected in MOE
OSN = readRDS("data/OR_OSN.rds") #single OR expressing OSN 
Simi = readRDS("data/OR_similarity_score_defult.rds") # protein similarity from pairwis alignment
OR_info = read.table("data/OR_loci.txt", header = T) # the file containing the info of OR: chr, genome position, cluster, class.
map_pre = readRDS("data/map_pre.rds") # reconstructed map
img = readJPEG("data/projection2.jpg") # the background of glomeruli map
source("DEseq2_fun.R") #for DESseq2 customized function
source("pair_analysis_harmony.R") # for silhouette coefficient analysis
source("plot_features.R") # for ploting the reconstructed map
source("Balanced_score.py") # for calculating balanced score

##### for the analysis with bulk RNAseq, MOR28 = Olfr1507, MOR174-9 = Olfr73, M72 = Olfr160 #####
##### 1b #####
plotPCA(vsd , intgroup = "genotype")
PCA_df = plotPCA(vsd , intgroup = "genotype", returnData = T)
p = ggplot(PCA_df,aes(x=PC1,y=PC2,color=group ))+geom_point(size = 5)+
  xlab(paste0("PC1:79% variance"))+ylab(paste0("PC2:15% variance"))+
  theme_bw()+theme(legend.position="bottom")
ggsave(plot=p,height=5.5,width=5,dpi=300, filename="Figure1b.pdf", useDingbats=FALSE)

##### 1d #####
res1 <- results(dds, contrast=c("genotype","M72","MOR174_9"))
res2 <- results(dds, contrast=c("genotype","MOR28","MOR174_9"))
res3 <- results(dds, contrast=c("genotype","MOR28","M72"))
res1_sig = get_sig(res1, thresh = 0.01)
res2_sig = get_sig(res2, thresh = 0.01)
res3_sig = get_sig(res3, thresh = 0.01)
genes = unique(c(rownames(res1_sig),rownames(res2_sig),rownames(res3_sig)))

coldata = as.data.frame(colData(dds))[,1,drop =F]

jpeg("Figure1d.jpeg", height = 8, width = 8, unit ="in", res = 300)
cols <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(15))
print(pheatmap(assay(vsd)[genes,], cluster_rows=TRUE, show_rownames=FALSE, color = cols,
               cluster_cols=TRUE, scale="row", border_color=NA, treeheight_row = 0, annotation_col = coldata))
dev.off()

##### 1c #####
set1 = rownames(res1_sig)
set2 = rownames(res2_sig)
set3 = rownames(res3_sig)
pdf("Figure1c.pdf", width = 5, height = 5)
draw.venn(set1,set2,set3, title = "", xtitle = "M72 vs. MOR174_9", ytitle = "OR28 vs. MOR174_9", ztitle = "MOR28 vs. M72", subtitle = "")
dev.off()

##### 1e #####
##Figure1e
ensembl=useMart("ensembl")
listDatasets(ensembl)  # find the dataset for your species 
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#listAttributes(ensembl)
foo <- getBM(attributes=c('entrezgene_id',
                          'mgi_symbol'),mart = ensembl)

backgroundgenes = rownames(assay(dds))
backgroundid = as.character(foo$entrezgene_id[foo$mgi_symbol%in%backgroundgenes])

genes = rownames(res1_sig)
genes_id = as.character(foo$entrezgene_id[foo$mgi_symbol%in%genes])
BPgenes1 = enrichGO(gene = genes_id, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff=0.05, pAdjustMethod="BH",readable = T, universe = backgroundid)

genes = rownames(res2_sig)
genes_id = as.character(foo$entrezgene_id[foo$mgi_symbol%in%genes])
BPgenes2 = enrichGO(gene = genes_id, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff=0.05, pAdjustMethod="BH",readable = T, universe = backgroundid)

genes = rownames(res3_sig)
genes_id = as.character(foo$entrezgene_id[foo$mgi_symbol%in%genes])
BPgenes3 = enrichGO(gene = genes_id, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff=0.05, pAdjustMethod="BH",readable = T, universe = backgroundid)

BPgenes1@result$Sample = "M72_MOR174_9"
BPgenes2@result$Sample = "MOR28_MOR174_9"
BPgenes3@result$Sample = "MOR28_M72"

results=data.frame(rbind(BPgenes1@result,BPgenes2@result,BPgenes3@result))
results = results[results$p.adjust<0.05,]
tmp = tmp = as.data.frame(table(results$Description))
tmp = tmp[tmp$Freq==3,]
results = results[results$Description%in%tmp$Var1,]
Item = c("sensory organ development",
         "Rho protein signal transduction",
         "regulation of cell-substrate adhesion",
         "regulation of actin filament-based process",
         "regulation of actin cytoskeleton organization",
         "neuron projection guidance",
         "negative regulation of cell adhesion",
         "chemotaxis",
         "cell-substrate adhesion",
         "cell projection assembly",
         "axonogenesis",
         "axon guidance")
results = results[results$Description%in%Item,]

p =ggplot(results, 
          aes(x = Sample, y = Description)) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), high = "blue" ,low="red") +
  ylab(NULL) +xlab(NULL)+
  ggtitle("GO biological pathway enrichment")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot=p,height=6,width=6,dpi=300, filename="Figure1e.pdf", useDingbats=FALSE)

##### 1f #####
genes = grep("Sema|Robo|Kirrel|Eph|Efn|Nrp|Plxn", rownames(assay(vsd)), value = T)
pdf("Figure1f.pdf", height = 10, width = 8)
pheatmap(assay(vsd)[genes,], cluster_rows=TRUE, show_rownames=T, color = cols,
         cluster_cols=T, scale="row", border_color=NA, treeheight_row = 0, annotation_col = coldata)
dev.off()

##### S1b #####
pl = list()
genes = c("Olfr73","Olfr160","Olfr1507")
genes = genes[order(nchar(genes), genes)]
for (i in 1:length(genes)) {
  pl[[i]] = plot_gene(dds, gene=genes[i], intgroup="genotype")+ theme(axis.text.x = element_blank(),axis.title.x=element_blank())+ggtitle(paste0(genes[i]))
}
p = wrap_plots(pl,nrow = 1)+ plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.key.size = unit(3,"line"), legend.text = element_text(size=10))
ggsave(plot=p,height=10,width=15,dpi=300, filename="FigureS1b.pdf", useDingbats=FALSE)

##### S1c #####
tpm = read.table("data/genes_expression_tpm.tsv", header = T, sep = "\t", row.names = 1)
tpm = tpm[,-1]
tpm = tpm[rowSums(tpm)>0,]
genes = c("Omp","Cnga2", "Gap43")
df = tpm[genes,]
df = melt(as.matrix(df))
colnames(df) = c("gene","sample","value")
df$intgroup = gsub("_.*","",df$sample)
p = ggplot(df, aes(x = gene, y = value, fill = intgroup)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~ sample, scales = "free_y")  +
  labs(y = 'TPM')+theme_bw()
ggsave(plot=p,height=10,width=12,dpi=300, filename="FigureS1c.pdf", useDingbats=FALSE)

##### S2b #####
jpeg("FigureS2b.jpeg", units = "in", width = 9, height = 9,res = 300)
DimPlot(MOE, group.by = "cell_type")+ NoLegend() +NoAxes()
dev.off()

##### S2c #####
marker_genes = c("Krt5","Ascl1","Neurod1","Gap43","Cnga2","Ascl3","Cyp2g1","Muc5b","S100b","Tmem212","C1qb","Ccr2","Igkc","Cd3g","Ngp","H1f5")
plots <- VlnPlot(MOE, features = marker_genes, group.by = "cell_type", pt.size = 0, ncol = 1,combine = F)
plots <- lapply(X = plots, FUN = function(x) x + theme(axis.text.y.left = element_text(size = 7),legend.position = "none"))
pdf("FigureS2c.pdf", width = 12, height = 25)
wrap_plots(plots = plots, ncol = 1)
dev.off()

##### S2d #####
df = readRDS("data/MOE_markers.rds")
top10 <- df %>% purrr::map_df(rev) %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
p = DotPlot(object = MOE, features = unique(as.character(top10$gene)),dot.scale = 5)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_gradient2(low = "#0571B0", mid = "white", high = "#CA0020", midpoint = 0)
ggsave(plot=p,height=5,width=20,dpi=300, filename="FigureS2a.pdf", useDingbats=FALSE)

##### 2a #####
ORs = sample(unique(OSN$OR_identity), 10)
#the OR random selected was "Olfr728","Olfr741","Olfr1260","Olfr390","Olfr1440","Olfr1030","Olfr571","Olfr1507","Olfr806", and"Olfr536"
my_color_palette <- hue_pal()(length(ORs))
Idents(OSN) = OSN$OR_identity
plots = list()
for (i in 1:length(ORs)) {
  cells = WhichCells(OSN, idents = c(paste0(ORs[i])))
  plots[[i]] = DimPlot(OSN, cells.highlight = cells, cols.highlight = my_color_palette[i])+NoAxes()+NoLegend()+ggtitle(ORs[i])
}
jpeg("Figure2a.jpeg", width = 15, height = 6, res = 300, units = "in")
wrap_plots(plots, ncol = 5)
dev.off()

##### S3a #####
markers = c("Acsm4","Nrp2","Plxna1","Nrp1")
plots = list()
for (i in 1:length(markers)) {
  plots[[i]] = FeaturePlot(OSN, features = markers[i], order = T, cols = cols)+NoAxes()+NoLegend()+ggtitle("")
}
jpeg("FigureS3a.jpeg", width = 6, height = 6, res = 300, units = "in")
wrap_plots(plots, ncol = 2)
dev.off()

##### S3b #####
metadf = OSN@meta.data
OSNrm = OSN@assays$RNA@counts
OSNrm = OSNrm[grep("Olfr",rownames(OSNrm), value = T,invert = T ),]
OSNrm = CreateSeuratObject(OSNrm)
OSNrm = NormalizeData(OSNrm)
OSNrm = FindVariableFeatures(OSNrm)
OSNrm = ScaleData(OSNrm)
OSNrm = RunPCA(OSNrm)
OSNrm@meta.data = metadf
OSNrm = RunHarmony(OSNrm, group.by.vars = "orig.ident")
OSNrm <- RunUMAP(OSNrm, reduction = "harmony", dims = 1:12)
saveRDS(OSNrm,"OR_OSN_rm.rds")
my_color_palette <- hue_pal()(length(ORs))
Idents(OSNrm) = OSN$OR_identity
plots = list()
for (i in 1:length(ORs)) {
  cells = WhichCells(OSNrm, idents = c(paste0(ORs[i])))
  plots[[i]] = DimPlot(OSNrm, cells.highlight = cells, cols.highlight = my_color_palette[i])+NoAxes()+NoLegend()
}
jpeg("FigureS3b.jpeg", width = 15, height = 6, res = 300, units = "in")
wrap_plots(plots, ncol = 5)
dev.off()

##### 2b #####
Freq = as.data.frame(table(OSN$OR_identity))
Freq = Freq[order(Freq$Freq, decreasing = T),]
pair_input = as.character(Freq$Var1[1:10])

matrix_input = OSN
Ident_input = "OR_identity"
Idents(matrix_input) = matrix_input@meta.data[,Ident_input]
Pairs = combn(unique(pair_input),2)
UMAP_plot = list()
for (i in 1:ncol(Pairs)) {
  target_cell_A = WhichCells(matrix_input, idents = Pairs[1,i])
  target_cell_B = WhichCells(matrix_input, idents = Pairs[2,i])
  ncell = min(length(target_cell_A),length(target_cell_B))
  Pair_mtx = subset(matrix_input, cells = c(sample(target_cell_A,ncell),sample(target_cell_B,ncell)))
  Pair_mtx = NormalizeData(Pair_mtx, verbose = FALSE)
  Pair_mtx = FindVariableFeatures(Pair_mtx, verbose = FALSE)
  Pair_mtx = ScaleData(Pair_mtx, verbose = FALSE)
  Pair_mtx <- RunPCA(Pair_mtx, verbose = FALSE)
  Pair_mtx <- RunHarmony(Pair_mtx, group.by.vars = "orig.ident")
  Pair_mtx <- RunUMAP(Pair_mtx, reduction = "harmony",verbose = FALSE, dims = 1:50)
  UMAP_plot[[i]] = DimPlot(Pair_mtx, reduction = "umap", group.by = Ident_input, pt.size = .5, cols = c("#35A2FF","#FF6A98"))+
    NoLegend()+NoAxes()+theme(plot.background = element_rect(colour = "snow3", size = .5))+ggtitle("")
}
jpeg("Figure2b.jpeg", width = 55, height = 55, units = "in", res = 300)
my_layout <- rbind(c(NA,1:9), c(rep(NA,2),10:17), c(rep(NA,3),18:24), c(rep(NA,4),25:30),c(rep(NA,5),31:35),
                     c(rep(NA,6),36:39),c(rep(NA,7),40:42),c(rep(NA,8),43:44),c(rep(NA,9),45))
grid.arrange(grobs = UMAP_plot, layout_matrix = my_layout)
dev.off()

##### S3c #####
matrix_input = OSNrm
Ident_input = "OR_identity"
Idents(matrix_input) = matrix_input@meta.data[,Ident_input]
Pairs = combn(unique(pair_input),2)
UMAP_plot = list()
for (i in 1:ncol(Pairs)) {
  target_cell_A = WhichCells(matrix_input, idents = Pairs[1,i])
  target_cell_B = WhichCells(matrix_input, idents = Pairs[2,i])
  ncell = min(length(target_cell_A),length(target_cell_B))
  Pair_mtx = subset(matrix_input, cells = c(sample(target_cell_A,ncell),sample(target_cell_B,ncell)))
  Pair_mtx = NormalizeData(Pair_mtx, verbose = FALSE)
  Pair_mtx = FindVariableFeatures(Pair_mtx, verbose = FALSE)
  Pair_mtx = ScaleData(Pair_mtx, verbose = FALSE)
  Pair_mtx <- RunPCA(Pair_mtx, verbose = FALSE)
  Pair_mtx <- RunHarmony(Pair_mtx, group.by.vars = "orig.ident")
  Pair_mtx <- RunUMAP(Pair_mtx, reduction = "harmony",verbose = FALSE, dims = 1:50)
  UMAP_plot[[i]] = DimPlot(Pair_mtx, reduction = "umap", group.by = Ident_input, pt.size = .5, cols = c("#35A2FF","#FF6A98"))+
    NoLegend()+NoAxes()+theme(plot.background = element_rect(colour = "snow3", size = .5))+ggtitle("")
}
jpeg("FigureS3c.jpeg", width = 55, height = 55, units = "in", res = 300)
my_layout <- rbind(c(NA,1:9), c(rep(NA,2),10:17), c(rep(NA,3),18:24), c(rep(NA,4),25:30),c(rep(NA,5),31:35),
                   c(rep(NA,6),36:39),c(rep(NA,7),40:42),c(rep(NA,8),43:44),c(rep(NA,9),45))
grid.arrange(grobs = UMAP_plot, layout_matrix = my_layout)
dev.off()

##### S3d #####
OSNrm$OR_identity = sample(OSNrm$OR_identity)
matrix_input = OSNrm
Ident_input = "OR_identity"
Idents(matrix_input) = matrix_input@meta.data[,Ident_input]
Pairs = combn(unique(pair_input),2)
UMAP_plot = list()
for (i in 1:ncol(Pairs)) {
  target_cell_A = WhichCells(matrix_input, idents = Pairs[1,i])
  target_cell_B = WhichCells(matrix_input, idents = Pairs[2,i])
  ncell = min(length(target_cell_A),length(target_cell_B))
  Pair_mtx = subset(matrix_input, cells = c(sample(target_cell_A,ncell),sample(target_cell_B,ncell)))
  Pair_mtx = NormalizeData(Pair_mtx, verbose = FALSE)
  Pair_mtx = FindVariableFeatures(Pair_mtx, verbose = FALSE)
  Pair_mtx = ScaleData(Pair_mtx, verbose = FALSE)
  Pair_mtx <- RunPCA(Pair_mtx, verbose = FALSE)
  Pair_mtx <- RunHarmony(Pair_mtx, group.by.vars = "orig.ident")
  Pair_mtx <- RunUMAP(Pair_mtx, reduction = "harmony",verbose = FALSE, dims = 1:50)
  UMAP_plot[[i]] = DimPlot(Pair_mtx, reduction = "umap", group.by = Ident_input, pt.size = .5, cols = c("#35A2FF","#FF6A98"))+
    NoLegend()+NoAxes()+theme(plot.background = element_rect(colour = "snow3", size = .5))+ggtitle("")
}
jpeg("FigureS3d.jpeg", width = 55, height = 55, units = "in", res = 300)
my_layout <- rbind(c(NA,1:9), c(rep(NA,2),10:17), c(rep(NA,3),18:24), c(rep(NA,4),25:30),c(rep(NA,5),31:35),
                   c(rep(NA,6),36:39),c(rep(NA,7),40:42),c(rep(NA,8),43:44),c(rep(NA,9),45))
grid.arrange(grobs = UMAP_plot, layout_matrix = my_layout)
dev.off()

##### 2c #####
OR654 = as.character(unique(Freq$Var1[Freq$Freq>=7]))
OSN654 = subset(OSN, ident = OR654)
OR654 = CreateSeuratObject(OR654) %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
saveRDS(OR654, "data/OR654.rds")

OR654rm = subset(OSNrm, ident = OR654)
OR654rm = CreateSeuratObject(OR654rm) %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
saveRDS(OR654rm, "data/OR654rm.rds")

OR654rmS = OR654rm
OR654rmS$OR_identity = sample(OR654rmS$OR_identity)
saveRDS(OR654rmS, "data/OR654rm.rds")
#Each of above 3 inputs was used to run the Sil.R on hpcc for parallelizing 213531 combinations from 654 groups of OR
#The result silhouette coefficent was stored as "OR_",start_n,"_",end_n,".rds", start_n and end_n refer to the subset of combinations.
#all the data were downloaded to data/

setwd("data/sil/")
OR = list.files("OR/",pattern = "0.rds|1.rds",full.names = T)
RM = list.files("Removed/",pattern = "0.rds|1.rds",full.names = T)
SH = list.files("shuffle/",pattern = "0.rds|1.rds",full.names = T)

OR_v = c()
for (i in 1:length(OR)) {
  OR_v = c(OR_v, readRDS(OR[i]))
}

RM_v = c()
for (i in 1:length(RM)) {
  RM_v = c(RM_v, readRDS(RM[i]))
}

SH_v = c()
for (i in 1:length(SH)) {
  SH_v = c(SH_v, readRDS(SH[i]))
}

df = as.data.frame(cbind(c(OR_v,RM_v,SH_v), c(rep("OR", length(OR_v)), rep("RM", length(RM_v)), rep("SH", length(SH_v)))))
df$V1 = as.numeric(df$V1)
df$V2 = as.factor(df$V2)
pairwise_t_test(data = df, V1 ~ V2, p.adjust.method = "BH") #check individual combination for p value
pdf("Figures2c.pdf", width = 4, height = 6)
ggboxplot(df, x = "V2", y="V1",
          fill = "V2", palette = c("#00AFBB", "#E7B800", "#FC4E07"),xlab = " ", ylab = "Silhouette Coefficient")+NoLegend()+
  stat_compare_means(method = "anova", label.y = 1)+      # check global p-value
  scale_x_discrete(labels=c("With OR", "Without OR","Random"))
dev.off()
write.table(df, "F2c_result.txt")

##### 2d #####
OR10 = Freq$Var1[1:10]
Idents(OSN) = OSN$OR_identity
OR10 = as.character(OR10)
OSN10 = subset(OSN, idents = OR10)
subdf = OSN10@meta.data
OSN10 = OSN10@assays$RNA@counts
OSN10 = OSN10[grep("Olfr",rownames(OSN10), invert = T),]
OSN10 = CreateSeuratObject(OSN10)
OSN10 = NormalizeData(OSN10)
OSN10 = FindVariableFeatures(OSN10, nfeatures = 3000)
OSN10 = ScaleData(OSN10)
OSN10 = RunPCA(OSN10, npcs = 50)
OSN10@meta.data = subdf
OSN10 = RunHarmony(OSN10, group.by.vars = "orig.ident")
OSN10 = FindNeighbors(OSN10, reduction = "harmony")
OSN10 = RunSPCA(OSN10, graph = "RNA_snn")
df = as.data.frame(OSN10@reductions$spca@cell.embeddings)
df$observed = OSN10$OR_identity
df$barcodes = rownames(df)
write.table(df, "data/svm_10_rm_all_50SPCs_3000G.csv", sep = ",")
saveRDS(OSN10,"data/OR10_rm_all_50SPCs_3000G.rds")

OR654 = as.character(unique(Freq$Var1[Freq$Freq>=7]))
Idents(OSN) = OSN$OR_identity
OSN654 = subset(OSN, idents = OR654)
subdf = OSN654@meta.data
OSN654 = OSN654@assays$RNA@counts
OSN654 = OSN654[grep("Olfr",rownames(OSN654), invert = T),]
OSN654 = CreateSeuratObject(OSN654)
OSN654 = NormalizeData(OSN654)
OSN654 = FindVariableFeatures(OSN654, nfeatures = 3000)
OSN654 = ScaleData(OSN654)
OSN654 = RunPCA(OSN654, npcs = 50)
OSN654@meta.data = subdf
OSN654 = RunHarmony(OSN654, group.by.vars = "orig.ident")
OSN654 = FindNeighbors(OSN654, reduction = "harmony")
OSN654 = RunSPCA(OSN654, graph = "RNA_snn")
df = as.data.frame(OSN654@reductions$spca@cell.embeddings)
df$observed = OSN654$OR_identity
df$barcodes = rownames(df)
write.table(df, "data/svm_654_rm_all_50SPCs_3000G.csv", sep = ",")
saveRDS(OSN654,"data/OR654_rm_all_50SPCs_3000G.rds")
# above 2 csv file was used as the input of Balanced_SVC.py for 100 repeats, ran on the hpcc, the results were downlaoded to data/
#10
OR = list.files("data/clf_10/", pattern = "results", full.names = T)
Sh = list.files("data/clf_10_shuffled/", pattern = "results", full.names = T)

OR_list = list()
for (i in 1:length(OR)) {
  OR_list[[i]] = read.csv(OR[i])
}
OR_list = do.call(rbind.data.frame,OR_list)


Sh_list = list()
for (i in 1:length(Sh)) {
  Sh_list[[i]] = read.csv(Sh[i])
}
Sh_list = do.call(rbind.data.frame,Sh_list)
OR_list$group = "OR"
Sh_list$group = "sh"

OR_list$test = factor(OR_list$test, levels = unique(OR_list$test))
OR_list$pred = factor(OR_list$pred, levels = unique(OR_list$test))
Sh_list$test = factor(Sh_list$test, levels = unique(Sh_list$test))
Sh_list$pred = factor(Sh_list$pred, levels = unique(Sh_list$test))

OSNdf = OSN@meta.data

OR_list$UMI = OSNdf[OR_list$barcodes,"nCount_RNA"]
Sh_list$UMI = OSNdf[Sh_list$barcodes,"nCount_RNA"]

grouped_acc <- function(df, groups) {
  acc_v = c()
  for (i in unique(df[,groups])) {
    tmp = df[df[,groups]==i,]
    acc_v = c(acc_v,mean(balanced_score(tmp$test,tmp$pred)))
  }
  return(acc_v)
}

merged = rbind.data.frame(OR_list, Sh_list)
write.table(merged,"F2d_10_raw_result.txt")

th = c(2600, 4000, 6000, 8000)
OR_acc_list = list()
Sh_acc_list = list()
for (z in 1:length(th)) {
  th_OR_list = OR_list[OR_list$UMI>=th[z],]
  th_Sh_list = Sh_list[Sh_list$UMI>=th[z],]
  OR_acc_list[[z]] = grouped_acc(th_OR_list,"rep")
  Sh_acc_list[[z]] = grouped_acc(th_Sh_list,"rep")
}

OR_acc = cbind(unlist(OR_acc_list),c(rep(th[1],length(OR_acc_list[[1]])),
                                     rep(th[2],length(OR_acc_list[[2]])),
                                     rep(th[3],length(OR_acc_list[[3]])),
                                     rep(th[4],length(OR_acc_list[[4]]))),rep("OR",length(unlist(Sh_acc_list))))

Sh_acc = cbind(unlist(Sh_acc_list),c(rep(th[1],length(Sh_acc_list[[1]])),
                                     rep(th[2],length(Sh_acc_list[[2]])),
                                     rep(th[3],length(Sh_acc_list[[3]])),
                                     rep(th[4],length(Sh_acc_list[[4]]))),rep("Shuffle",length(unlist(Sh_acc_list))))
df = rbind.data.frame(OR_acc,Sh_acc)

colnames(df) = c("Balanced_accuracy", "Threshold", "Group")
write.table(df, "data/F2d_10_result.txt")

df$Balanced_accuracy = as.numeric(df$Balanced_accuracy)
df$Balanced_accuracy = df$Balanced_accuracy*100
df = read.table("data/F2d_10_result.txt")
df$Threshold = as.character(df$Threshold)
df$Group = factor(df$Group, levels = c("Shuffle", "OR"))

pdf("Figure2d_left.pdf", width = 5, height = 7)
ggplot(df, aes(x = Group, y = Balanced_accuracy, fill = Threshold, color = Threshold))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+ylab("Balanced accuracy")+stat_summary(fun=mean, colour="red", geom="text",show.legend = FALSE, aes(x = Group, y = Balanced_accuracy, group = Threshold,label=round(..y.., digits=2)),
                                                    position = position_dodge(width = .75), vjust=-4)+ theme_classic()+
  scale_color_manual(values = c("#EB4136", "#009044", "#1C75B8", "#7F3F94"))+
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = expansion(mult = c(0.15, 0.15)))
dev.off()

for (i in 1:length(unique(df$Threshold))) {
  print(compare_means(Balanced_accuracy ~ Group, data = df[df$Threshold==unique(df$Threshold)[i],], method = "t.test"))
}

#654
OR = list.files("data/clf/", pattern = "results", full.names = T)
Sh = list.files("data/clf_shuffled/", pattern = "results", full.names = T)

OR_list = list()
for (i in 1:length(OR)) {
  OR_list[[i]] = read.csv(OR[i])
}
OR_list = do.call(rbind.data.frame,OR_list)


Sh_list = list()
for (i in 1:length(Sh)) {
  Sh_list[[i]] = read.csv(Sh[i])
}
Sh_list = do.call(rbind.data.frame,Sh_list)
OR_list$group = "OR"
Sh_list$group = "sh"

OR_list$test = factor(OR_list$test, levels = unique(OR_list$test))
OR_list$pred = factor(OR_list$pred, levels = unique(OR_list$test))
Sh_list$test = factor(Sh_list$test, levels = unique(Sh_list$test))
Sh_list$pred = factor(Sh_list$pred, levels = unique(Sh_list$test))


OR_list$UMI = OSNdf[OR_list$barcodes,"nCount_RNA"]
Sh_list$UMI = OSNdf[Sh_list$barcodes,"nCount_RNA"]

grouped_acc <- function(df, groups) {
  acc_v = c()
  for (i in unique(df[,groups])) {
    tmp = df[df[,groups]==i,]
    acc_v = c(acc_v,mean(balanced_score(tmp$test,tmp$pred)))
  }
  return(acc_v)
}

merged = rbind.data.frame(OR_list, Sh_list)
write.table(merged,"data/F2d_654_raw_result.txt")

th = c(2600, 4000, 6000, 8000)
OR_acc_list = list()
Sh_acc_list = list()
for (z in 1:length(th)) {
  th_OR_list = OR_list[OR_list$UMI>=th[z],]
  th_Sh_list = Sh_list[Sh_list$UMI>=th[z],]
  OR_acc_list[[z]] = grouped_acc(th_OR_list,"rep")
  Sh_acc_list[[z]] = grouped_acc(th_Sh_list,"rep")
}

OR_acc = cbind(unlist(OR_acc_list),c(rep(th[1],length(OR_acc_list[[1]])),
                                     rep(th[2],length(OR_acc_list[[2]])),
                                     rep(th[3],length(OR_acc_list[[3]])),
                                     rep(th[4],length(OR_acc_list[[4]]))),rep("OR",length(unlist(Sh_acc_list))))

Sh_acc = cbind(unlist(Sh_acc_list),c(rep(th[1],length(Sh_acc_list[[1]])),
                                     rep(th[2],length(Sh_acc_list[[2]])),
                                     rep(th[3],length(Sh_acc_list[[3]])),
                                     rep(th[4],length(Sh_acc_list[[4]]))),rep("Shuffle",length(unlist(Sh_acc_list))))
df = rbind.data.frame(OR_acc,Sh_acc)

colnames(df) = c("Balanced_accuracy", "Threshold", "Group")
df$Balanced_accuracy = as.numeric(df$Balanced_accuracy)
df$Balanced_accuracy = df$Balanced_accuracy*100
write.table(df, "data/F2d_654_result.txt")
df = read.table("data/F2d_654_result.txt")
df$Threshold = as.character(df$Threshold)
df$Group = factor(df$Group, levels = c("Shuffle", "OR"))

pdf("Figure2d_right.pdf", width = 5, height = 7)
ggplot(df, aes(x = Group, y = Balanced_accuracy, fill = Threshold, color = Threshold))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+ylab("Balanced accuracy")+stat_summary(fun=mean, colour="red", geom="text",show.legend = FALSE, aes(x = Group, y = Balanced_accuracy, group = Threshold,label=round(..y.., digits=2)),
                                                    position = position_dodge(width = .75), vjust=-4)+ theme_classic()+
  scale_color_manual(values = c("#EB4136", "#009044", "#1C75B8", "#7F3F94"))+
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = expansion(mult = c(0.15, 0.15)))
dev.off()

for (i in 1:length(unique(df$Threshold))) {
  print(compare_means(Balanced_accuracy ~ Group, data = df[df$Threshold==unique(df$Threshold)[i],], method = "t.test"))
}

##### S4a #####
# Get the biological process GO term with at least 150 genes
k = keys(org.Mm.eg.db, keytype="GO")
df = select(org.Mm.eg.db, keys=k, columns=c("ONTOLOGY"), keytype="GO")
df = df[df$ONTOLOGY=="BP",]
BPGO = df$GO
k = unique(BPGO)
df = select(org.Mm.eg.db, keys=k, columns=c("SYMBOL"), keytype="GO")
df = df[,c("GO","SYMBOL")]
df = df[df$SYMBOL%in%rownames(OSN),]
x = df %>%
  group_by(GO) %>%
  summarise(count = n_distinct(SYMBOL))
x = x[x$count>=150,]
ref = x$GO
saveRDS(ref, "Go_ID.rds")
# get 10 OR goups
OR10 = Freq$Var1[1:10]
Idents(OSN) = OSN$OR_identity
OR10 = as.character(OR10)
OSN10 = subset(OSN, idents = OR10)
subdf = OSN10@meta.data
OSN10 = OSN10@assays$RNA@counts
OSN10 = OSN10[grep("Olfr",rownames(OSN10), invert = T),]
OSN10 = CreateSeuratObject(OSN10)
OSN10 = NormalizeData(OSN10)
OSN10 = FindVariableFeatures(OSN10, nfeatures = 3000)
saveRDS(OSN10,"data/OSN10.rds")

# The OSN10.rds and Go_ID.rds was used as the input of GO_150_input_generation.R to generate the input of classifier
# Each GO, 100 different input was generated by ramdom sampling the genes in that GO term as the variable genes for PCA
# The Balanced_SVC_GO.py was used to run the prediction result of each GO. the results were downloaded to the data/ 
goterms <- as.data.frame(Term(GOTERM))
gonames = goterms[ref,]
ref = gsub(":","_",ref)

file_list = list()
for (i in 1:length(ref)) {
  if (length(list.files(path = "data/clf_10_GO_150/results/",
                        pattern = ref[i]))==0) {
    next
  }
  myfiles = lapply(list.files(path = "data/clf_10_GO_150/results/",
                              pattern = ref[i]), read.csv)
  myfiles = do.call(rbind.data.frame, myfiles)
  myfiles$X = gonames[i]
  file_list[[i]] = myfiles
}
file_list = do.call(rbind.data.frame,file_list)
file_list = file_list[!(file_list$X=="biological_process"),]
file_list$X = capitalize(as.character(file_list$X))

pdf("FigureS4a.pdf", width = 20, height = 8)
ggplot(file_list, aes(x = fct_reorder(X, balanced_score, .desc = T), y = balanced_score))+geom_boxplot()+
  xlab(NULL)+ylab("Balanced Accuracy")+ylim(c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))+
  theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), legend.position = "none")
dev.off()

##### 2e and S4b #####
#reload Simi and OR_info
merged = read.table("data/F2d_654_raw_result.txt")
OR_info = OR_info[!duplicated(OR_info$gene.name),]
rownames(OR_info) = OR_info$gene.name
OSN$OR_cluster = OR_info[OSN$OR_identity,"OR_cluster"]
OSN$chr = OR_info[OSN$OR_identity,"chr"]
OSN$class = OR_info[OSN$OR_identity,"class"]
saveRDS(OSN,"data/OR_OSN.rds")
OSN = OSN[OSN$OR_identity%in%unique(merged$test),]
Simi = Simi[unique(merged$test),unique(merged$test)]
wrong = merged[!(merged$test==merged$pred),]

chr_vac = c()
cluster_vac = c()
similarity_vac = c()
for (i in 1:nrow(wrong)) {
  print(i)
  OR1 = wrong$test[i]
  OR2 = wrong$pred[i]
  chr1 = OR_info[OR1,"chr"]
  chr2 = OR_info[OR2,"chr"]
  cl1 = OR_info[OR1,"OR_cluster"]
  cl2 = OR_info[OR2,"OR_cluster"]
  similarity_vac[i] = Simi[OR1,OR2]
  chr_vac[i] = as.numeric(chr1==chr2)
  cluster_vac[i] = as.numeric(cl1==cl2)
}
wrong$chr = chr_vac
wrong$OR_cluster = cluster_vac
wrong$similarity = similarity_vac
wrong = wrong[!(is.na(wrong$OR_cluster)),]
wrong$above_0.75 = wrong$similarity>=0.75
wrong$above_0.65 = wrong$similarity>=0.65
wrong$above_0.55 = wrong$similarity>=0.55
wrong$above_0.45 = wrong$similarity>=0.45
wrong$above_0.35 = wrong$similarity>=0.35
wrong$above_0.25 = wrong$similarity>=0.25
wrong$above_0.15 = wrong$similarity>=0.15
write.table(wrong, "data/prediction_error_info.txt", sep = "\t", col.names = NA)
wrong = read.table("data/prediction_error_info.txt", header = T)
wrong$Bin0_1 = wrong$similarity>=0&wrong$similarity<0.1
wrong$Bin1_2 = wrong$similarity>=0.1&wrong$similarity<0.2
wrong$Bin2_3 = wrong$similarity>=0.2&wrong$similarity<0.3
wrong$Bin3_4 = wrong$similarity>=0.3&wrong$similarity<0.4
wrong$Bin4_5 = wrong$similarity>=0.4&wrong$similarity<0.5
wrong$Bin5_6 = wrong$similarity>=0.5&wrong$similarity<0.6
wrong$Bin6_7 = wrong$similarity>=0.6&wrong$similarity<0.7
wrong$Bin7_8 = wrong$similarity>=0.7&wrong$similarity<0.8
wrong$Bin8_9 = wrong$similarity>=0.8&wrong$similarity<0.9
wrong$Bin9_0 = wrong$similarity>=0.9&wrong$similarity<1
wrong$above_0.75 = wrong$similarity>=0.75

summary_df = wrong %>% group_by(group, rep) %>% summarise(chr=mean(chr),
                                                          OR_cluster=mean(OR_cluster),
                                                          above_0.75 = mean(above_0.75),
                                                          Bin0_1 = mean(Bin0_1),
                                                          Bin1_2 = mean(Bin1_2),
                                                          Bin2_3 = mean(Bin2_3),
                                                          Bin3_4 = mean(Bin3_4),
                                                          Bin4_5 = mean(Bin4_5),
                                                          Bin5_6 = mean(Bin5_6),
                                                          Bin6_7 = mean(Bin6_7),
                                                          Bin7_8 = mean(Bin7_8),
                                                          Bin8_9 = mean(Bin8_9),
                                                          Bin9_0 = mean(Bin9_0))
summary_df_s = summary_df[summary_df$group=="sh",]
summary_df_t = summary_df[summary_df$group=="OR",]
bins = c("Bin0_1","Bin1_2","Bin2_3","Bin3_4","Bin4_5","Bin5_6","Bin6_7","Bin7_8","Bin8_9","Bin9_0")
bin_list = list()
for (i in 1:length(bins)) {
  tmp = as.data.frame(summary_df)
  reps = unique(summary_df$rep[which(summary_df[,paste0(bins[i])]==0)])
  tmp = tmp[!(tmp$rep%in%reps),]
  tmps = tmp[tmp$group=="sh",]
  tmpt = tmp[tmp$group=="OR",]
  bin_list[[i]] = tmpt[,paste0(bins[i])]/tmps[,paste0(bins[i])]
}

len = unlist(lapply(bin_list, length))
namev = c()
for (i in 1:length(len)) {
  namev = c(namev,rep(bins[i],len[i]))
}

new_df = as.data.frame(cbind(fold_change = c((summary_df_t$chr/summary_df_s$chr),(summary_df_t$OR_cluster/summary_df_s$OR_cluster), (summary_df_t$above_0.75/summary_df_s$above_0.75)),
                             group = c(rep("chr",1000),rep("OR_cluster",1000),rep("above_0.75",1000))))
new_df$fold_change = as.numeric(new_df$fold_change)
pairwise_t_test(data = new_df, fold_change ~ group, p.adjust.method = "BH")
new_df$fold_change = log2(as.numeric(new_df$fold_change))
pdf("Figure2e.pdf", width = 5, height = 7)
ggboxplot(as.data.frame(new_df), x = "group", y = "fold_change", color = "group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "none")+
  ylab("fold change")+scale_color_manual(values = c("#EB4136", "#009044", "#1C75B8"))
dev.off()


new_df = as.data.frame(cbind(fold_change = c(unlist(bin_list)),
                             group = c(namev)))
new_df$fold_change = log2(as.numeric(new_df$fold_change))
pdf("FigureS4b.pdf", width = 5, height = 7)
ggboxplot(as.data.frame(new_df), x = "group", y = "fold_change", color = "group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "none")+
  ylab("fold change")
dev.off()

#reload Simi
Simi = as.data.frame(as.matrix(Simi))

OSN654 = readRDS("data/OR654_rm_all_50SPCs_3000G.rds")
df = as.data.frame(OSN654@reductions$pca@cell.embeddings)
df_mean = aggregate(df[,1:50], list(OSN654$OR_identity), mean)
rownames(df_mean) = df_mean$Group.1
df_mean = df_mean[,-1]
cor_df = cor(t(df_mean))
cor_df = as.data.frame(cor_df)
rownames(cor_df) = rownames(df_mean)
colnames(cor_df) = rownames(df_mean)
cor_df = melt(as.matrix(cor_df))
cor_df2 = cor_df[!duplicated(t(apply(cor_df, 1, sort))),]
cor_df2 = cor_df2[!(cor_df2$Var1==cor_df2$Var2),]
cor_df2 = as.data.frame(cor_df2)

Simi = Simi[rownames(Simi)%in%rownames(df_mean),colnames(Simi)%in%rownames(df_mean)]
Simi = Simi[match(rownames(df_mean), rownames(Simi)),match(rownames(df_mean), colnames(Simi))]
Simi2 = melt(as.matrix(Simi))
Simi2 = Simi2[!duplicated(t(apply(Simi2, 1, sort))),]
Simi2 = Simi2[!(Simi2$Var1==Simi2$Var2),]

all(cor_df2$Var1 == Simi2$Var1)
all(cor_df2$Var2 == Simi2$Var2)

Simi2$corelation = cor_df2$value
colnames(Simi2) = c("OR1","OR2","protein_sequence_similarity","correlation")


MinMeanSEMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

Simi2$groups = "0.0~0.1"
Simi2[Simi2$protein_sequence_similarity<0.2&Simi2$protein_sequence_similarity>=0.1,]$groups = "0.1~0.2"
Simi2[Simi2$protein_sequence_similarity<0.3&Simi2$protein_sequence_similarity>=0.2,]$groups = "0.2~0.3"
Simi2[Simi2$protein_sequence_similarity<0.4&Simi2$protein_sequence_similarity>=0.3,]$groups = "0.3~0.4"
Simi2[Simi2$protein_sequence_similarity<0.5&Simi2$protein_sequence_similarity>=0.4,]$groups = "0.4~0.5"
Simi2[Simi2$protein_sequence_similarity<0.6&Simi2$protein_sequence_similarity>=0.5,]$groups = "0.5~0.6"
Simi2[Simi2$protein_sequence_similarity<0.7&Simi2$protein_sequence_similarity>=0.6,]$groups = "0.6~0.7"
Simi2[Simi2$protein_sequence_similarity<0.8&Simi2$protein_sequence_similarity>=0.7,]$groups = "0.7~0.8"
Simi2[Simi2$protein_sequence_similarity<0.9&Simi2$protein_sequence_similarity>=0.8,]$groups = "0.8~0.9"
Simi2[Simi2$protein_sequence_similarity>=0.9,]$groups = "0.9~1.0"

pdf("FigureS4c.pdf", width = 4, height = 5)
ggplot(Simi2, aes(x=groups, y=correlation, color = groups))+
  geom_errorbar(stat="summary", fun.data="mean_se",width=.2)+theme(legend.position = "none")+ylab("transcriptome correlation")+xlab("protein sequence similarity score")+theme_classic()+theme(legend.position = "none")
dev.off()

##### 2f #####
#reload OSN, Simi, OR_info
OR_set1 = unique(OSN$OR_identity)
OR_set2 = unique(OR_info$gene.name)
OR_set3 = unique(rownames(as.matrix(Simi)))
common_OR = intersect(OR_set1,intersect(OR_set2,OR_set3))

OR_info = OR_info[OR_info$gene.name%in%common_OR,]
OR_info = OR_info[!duplicated(OR_info$gene.name),]
OR_info = OR_info[!is.na(OR_info$OR_cluster),]
Simi2 = Simi[common_OR,common_OR]

Idents(OSN) = OSN$OR_identity
OSN = subset(OSN, idents = common_OR)
metadf = OSN@meta.data
EM = OSN@assays$RNA@counts
EM = EM[-grep("Olfr", rownames(EM)),]
OSN = CreateSeuratObject(EM)
OSN@meta.data = metadf
OSN = NormalizeData(OSN, verbose = FALSE)
OSN = FindVariableFeatures(OSN, verbose = FALSE)
OSN = ScaleData(OSN, verbose = FALSE)
OSN = RunPCA(OSN, verbose = FALSE)
OSN@meta.data = metadf
OSN = RunHarmony(OSN, group.by.vars = "orig.ident",verbose = F)
df = as.data.frame(OSN@reductions$harmony@cell.embeddings)
df$OR_identity = metadf$OR_identity
df_mean = aggregate(df[,1:50], list(df$OR_identity), mean)
rownames(df_mean) = df_mean$Group.1
df_mean = df_mean[,-1]
cor_df = cor(t(df_mean))

frequency_list = as.data.frame(table(OR_info$OR_cluster))
target_cluster = as.character(frequency_list$Var1[frequency_list$Freq>=10])

result_list = list()
for (i in 1:length(unique(target_cluster))) {
  ORs = unique(OR_info$gene.name[OR_info$OR_cluster==target_cluster[i]])
  subdf = OR_info[OR_info$gene.name%in%ORs,c("gene.name","start")]
  mat = as.matrix(Simi(subdf$start))
  colnames(mat) = subdf$gene.name
  rownames(mat) = subdf$gene.name
  longdf1 = melt(mat)
  longOR_info = melt(Simi2[ORs,ORs])
  longdf3 = melt(cor_df[ORs,ORs])
  merged = cbind(longdf1,longOR_info$value,longdf3$value)
  colnames(merged) = c("OR1","OR2","Simiance","similarity","correlation")
  merged = merged[!(merged$OR1==merged$OR2),]
  merged$cluster = target_cluster[i]
  result_list[[i]] = merged
}
all_result = do.call(rbind.data.frame,result_list)

write.table(all_result,"data/TTS_results.txt", quote = F, sep = "\t")

p1 = ggscatter(all_result, x = "Simiance", y = "similarity", size = 0, alpha = 0)+
  xlab("TSS Simiance")+ylab("Protein sequence similarity")+
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_viridis()+
  stat_cor(method = "pearson", label.y = 1.1)+ theme(legend.position = "right")+labs(fill = "Counts")
x = cor.test(all_result$Simiance, all_result$similarity)
p2 = ggscatter(all_result, x = "Simiance", y = "correlation", size = 0, alpha = 0)+
  xlab("TSS Simiance")+ylab("Trancriptome correlation")+
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_viridis()+
  stat_cor(method = "pearson", label.y = 1.1)+ theme(legend.position = "right")+labs(fill = "Counts")

pdf("Figure2f.pdf", width = 10, height = 5)
p1|p2
dev.off()

##### S4d #####
#reload OSN, Simi, OR_info
OR_set1 = unique(OSN$OR_identity)
OR_set2 = unique(OR_info$gene.name)
OR_set3 = unique(rownames(as.matrix(Simi)))
common_OR = intersect(OR_set1,intersect(OR_set2,OR_set3))

Idents(OSN) = OSN$OR_identity
OSN = subset(OSN, idents = common_OR)
metadf = OSN@meta.data
EM = OSN@assays$RNA@counts
EM = EM[-grep("Olfr", rownames(EM)),]
OSN = CreateSeuratObject(EM)
OSN@meta.data = metadf
OSN = NormalizeData(OSN, verbose = FALSE)
OSN = FindVariableFeatures(OSN, verbose = FALSE)
OSN = ScaleData(OSN, verbose = FALSE)
OSN = RunPCA(OSN, verbose = FALSE)
OSN = RunHarmony(OSN, group.by.vars = "orig.ident",verbose = F)
df = as.data.frame(OSN@reductions$harmony@cell.embeddings)
df$OR_identity = metadf$OR_identity
colnames(df) = gsub("harmony","PC",colnames(df))
df_mean = aggregate(df[,1:50], list(df$OR_identity), mean)
rownames(df_mean) = df_mean$Group.1
df_mean = df_mean[,-1]
cor_df = cor(t(df_mean))
cor_df = as.data.frame(cor_df)

Simi = as.data.frame(as.matrix(Simi))
Simi = Simi[rownames(cor_df),rownames(cor_df)]
all(colnames(Simi) == colnames(cor_df))
all(rownames(Simi) == rownames(cor_df))

OR_info = as.data.frame(OR_info)
OR_info = OR_info[OR_info$gene.name%in%common_OR,]
OR_info = OR_info[!duplicated(OR_info$gene.name),]

breaksList = seq((-1), 1, by = 0.01)

OR_info = OR_info[order(OR_info$chr,OR_info$center, decreasing = F),]
cor_df2 = cor_df[OR_info$gene.name,OR_info$gene.name]
Simi2 = Simi[OR_info$gene.name,OR_info$gene.name]
cor_df2[lower.tri(cor_df2)] <- Simi2[lower.tri(Simi2)]
all(rownames(cor_df)==df$gene.name)
anndf = OR_info[,c("gene.name","chr","OR_cluster")]
rownames(anndf) = anndf$gene.name
anndf = anndf[,-1]
colnames(anndf) = c("chr","OR_cluster")

tmp_x = unique(anndf$chr)
anndf$chr = factor(anndf$chr, levels = tmp_x[order(nchar(tmp_x), tmp_x)])
tmp_y = unique(anndf$OR_cluster)
anndf$OR_cluster = factor(anndf$OR_cluster, levels = tmp_y[order(nchar(tmp_y), tmp_y)])
col1 <- hue_pal()(length(levels(anndf$chr)))
names(col1) <- levels(anndf$chr)
col2 <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(anndf$OR_cluster)))
names(col2) <- levels(anndf$OR_cluster)
ha <- HeatmapAnnotation(df = anndf, col = list(chr = col1, OR_cluster = col2), annotation_name_side = "left")

split_name = OR_info$chr
#jpeg("All_chr_transcriptome_protein_similarity_filtered.jpeg", height = 15, width = 17, units = "in", res = 300)
#Heatmap(cor_df2, cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha, row_split = split_name, column_split = split_name,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
#dev.off()

anndf2 = anndf[grep("chr10",anndf$chr),]
tmp_x = as.character(unique(anndf2$chr))
anndf2$chr = factor(anndf2$chr, levels = tmp_x[order(nchar(tmp_x), tmp_x)])
tmp_y = as.character(unique(anndf2$OR_cluster))
anndf2$OR_cluster = factor(anndf2$OR_cluster, levels = tmp_y[order(nchar(tmp_y), tmp_y)])
col1 <- hue_pal()(length(levels(anndf2$chr)))
names(col1) <- levels(anndf2$chr)
col2 <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(anndf2$OR_cluster)))
names(col2) <- levels(anndf2$OR_cluster)
ha <- HeatmapAnnotation(df = anndf2[,"OR_cluster",drop = F], col = list(OR_cluster = col2), annotation_name_side = "left")
OR_ch10 = cor_df2[which(split_name=="chr10"),which(split_name=="chr10")]
pdf("FigureS4d_chr10.pdf", height = 10, width = 12)
Heatmap(as.data.frame(OR_ch10), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
dev.off()


anndf2 = anndf[grep("chr2",anndf$chr),]
tmp_x = as.character(unique(anndf2$chr))
anndf2$chr = factor(anndf2$chr, levels = tmp_x[order(nchar(tmp_x), tmp_x)])
tmp_y = as.character(unique(anndf2$OR_cluster))
anndf2$OR_cluster = factor(anndf2$OR_cluster, levels = tmp_y[order(nchar(tmp_y), tmp_y)])
col1 <- hue_pal()(length(levels(anndf2$chr)))
names(col1) <- levels(anndf2$chr)
col2 <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(anndf2$OR_cluster)))
names(col2) <- levels(anndf2$OR_cluster)
ha <- HeatmapAnnotation(df = anndf2[,"OR_cluster",drop = F], col = list(OR_cluster = col2), annotation_name_side = "left")
OR_ch2 = cor_df2[which(split_name=="chr2"),which(split_name=="chr2")]
pdf("FigureS4d_chr2.pdf", height = 10, width = 12)
Heatmap(as.data.frame(OR_ch2), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
dev.off()


anndf2 = anndf[grep("chr7",anndf$chr),]
tmp_x = as.character(unique(anndf2$chr))
anndf2$chr = factor(anndf2$chr, levels = tmp_x[order(nchar(tmp_x), tmp_x)])
tmp_y = as.character(unique(anndf2$OR_cluster))
anndf2$OR_cluster = factor(anndf2$OR_cluster, levels = tmp_y[order(nchar(tmp_y), tmp_y)])
col1 <- hue_pal()(length(levels(anndf2$chr)))
names(col1) <- levels(anndf2$chr)
col2 <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(anndf2$OR_cluster)))
names(col2) <- levels(anndf2$OR_cluster)
ha <- HeatmapAnnotation(df = anndf2[,"OR_cluster",drop = F], col = list(OR_cluster = col2), annotation_name_side = "left")
OR_ch7 = cor_df2[which(split_name=="chr7"),which(split_name=="chr7")]
pdf("FigureS4d_chr7.pdf", height = 10, width = 12)
Heatmap(as.data.frame(OR_ch7), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
dev.off()

##### 3b #####
img = readJPEG("data/projection2.jpg")
dt = readRDS(paste0("data/OB2_slideseq/Slide11.rds"))
genes = c("Kctd12","Calb2", "Doc2g", "Cdhr1","Pcp4","Sox11")
for (i in 1:length(genes)) {
  jpeg(paste("Marker",genes[i],".jpeg", sep = "_"), res = 300, units = "in", width = 5, height = 5)
  print(FeaturePlot(dt, genes[i], order = T, reduction = "Spatial", max.cutoff = "q95", pt.size = 0.5)+NoAxes()+NoLegend()&scale_color_viridis(option = "C", direction = 1))
  dev.off()
}

##### 3c #####
set = c(1,1,1,2,2,2)
slide = c(6,9,10,11,12,15)
ORs = c("Olfr791","Olfr16","Olfr1018","Olfr479","Olfr517","Olfr1250")
for (i in 1:6) {
  dt = readRDS(paste0("data/OB",set[i],"_slideseq/Slide",slide[i],".rds"))
  EM = GetAssayData(dt, "counts")
  OR = ORs[i]
  EM = EM[OR,,drop = F]
  cell_list = colnames(EM)[colSums(as.matrix(EM))>0]
  p1 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+NoLegend()+NoAxes()& 
    scale_color_viridis(option = "D", direction = -1)
  p2 = p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =2)
  jpeg(paste0("OB",set[i],slide[i],OR,".jpeg"), res = 300, units = "in", width = 5, height = 5)
  print(p2)
  dev.off()
}

##### 3d #####
df = read.table("data/map_obs.txt", header = T)
pdf("Figure3d.pdf", width = 10, height = 10)
ggplot(df, aes(x = x, y = y, label = OR))+background_image(img)+geom_point()+
  geom_text_repel()+xlab("A---P")+ylab("V---D")+xlim(c(1.00000,21.29802))+ylim(c(-4584.6427,292.5073))+theme_classic()
dev.off()

##### 3e #####
df = read.table("data/map_654.txt")
colnames(df) = c("OR", "cell", "type", "x", "y")
ORs = df[df$type=="slide-seq","OR"]
df2 = df[df$OR%in%ORs,]

OR_candidates = c("Olfr47","Olfr788","Olfr360","Olfr16","Olfr598","Olfr1110")

plots = list()
for (i in 1:length(OR_candidates)) {
  plots[[i]] = ggplot(df2[df2$OR==OR_candidates[i],], aes(x = x, y = y, color = type))+
    background_image(img)+geom_point(size =2)+xlab("A---P")+ylab("V---D")+xlim(c(1.00000,21.29802))+ylim(c(-4584.6427,292.5073))+theme_classic()+NoLegend()
}
pdf("Figure3e.pdf", width = 30, height = 20)
wrap_plots(plots, ncol = 3)
dev.off()

##### 3f #####
df2 = df[df$OR%in%ORs,]
df_obs = df2[df2$type=="slide-seq",]
df_pre = df2[df2$type=="predicted",]

df_obs = df_obs[order(df_obs$OR),]
df_pre = df_pre[order(df_pre$OR),]
my_data = as.data.frame(cbind(predicted = df_pre$x,
                              observed = df_obs$x))

p1 = ggplot(my_data,aes(predicted, observed)) + geom_point(color = "black") + 
  geom_smooth(method=lm) + ggtitle(paste("A-P axis","R^2 =",round(R2_Score(df_pre$x,df_obs$x),2),"RMSE =",round(RMSE(df_pre$x,df_obs$x),2), sep = " ")) +
  xlab("Predecited") + ylab("Observed")+theme_classic()

my_data = as.data.frame(cbind(predicted = df_pre$y,
                              observed = df_obs$y))
p2 = ggplot(my_data,aes(predicted, observed)) + geom_point(color = "black") + 
  geom_smooth(method=lm) + ggtitle(paste("D-V axis","R^2 =",round(R2_Score(df_pre$y,df_obs$y),2),"RMSE =",round(RMSE(df_pre$y,df_obs$y),2), sep = " ")) +
  xlab("Predecited") + ylab("Observed")+theme_classic()

pdf("Figure3f.pdf", width = 8, height = 4)
p1+p2
dev.off()

##### S5a #####
dt = readRDS("data/OB2_slideseq/Slide11.rds")
jpeg("FigureS5a.jpeg", res = 300, units = "in", width = 5, height = 5)
FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
dev.off()

##### S5b #####
#LOOCV Homography pre calculated with python script Common_OR_transform_hold1
files =list.files(pattern = "transformed_coords")
dflist = list()
for (i in 1:length(files)) {
  df = read.table(files[i])
  colnames(df) = c("OR","x","y","rep")
  ORs = gsub("transformed_coords_1to2_|.txt","",files[i])
  df = df[df$OR==ORs,]
  dflist[[i]] = df
}
df = do.call(rbind.data.frame, dflist)
df1 = df[df$rep==1,]
df2 = df[df$rep==2,]
dfx = as.data.frame(cbind("set1" = df1$x, "set2" = df2$x))
dfy = as.data.frame(cbind("set1" = df1$y, "set2" = df2$y))
p1 = ggscatter(dfx, x = "set1", y = "set2", 
               add = "reg.line", conf.int = TRUE, add.params = list(color = "blue", fill = "lightblue"),
               cor.coef = TRUE, cor.method = "pearson", title = "A-P")
p2 = ggscatter(dfy, x = "set1", y = "set2", 
               add = "reg.line", conf.int = TRUE, add.params = list(color = "blue", fill = "lightblue"),
               cor.coef = TRUE, cor.method = "pearson", title = "D-V")


pdf("FigureS5b.pdf", width = 10, height = 5)
p1|p2 
dev.off()

##### S5C #####
dt = readRDS("data/OB1_slideseq/Slide8.rds")

EM = GetAssayData(dt, "counts")
OR = "Olfr414"
EM = EM[OR,,drop = F]
cell_list = colnames(EM)[colSums(as.matrix(EM))>0]
p1 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
p2 = p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =2)

dt = readRDS("~/Dropbox (UMass Medical School)/OB1_slideseq/Analysis/Slide9.rds")
EM = GetAssayData(dt, "counts")
OR = "Olfr16"
EM = EM[OR,,drop = F]
cell_list = colnames(EM)[colSums(as.matrix(EM))>0]
p3 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
p4 = p3+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =2)

jpeg("FigureS5c.jpeg", res = 300, units = "in", width = 10, height = 5)
p2|p4
dev.off()

##### 4a #####
map_pre = readRDS("data/map_pre.rds")
pdf("Figure4a.pdf", width = 5, height = 5)
ggplot(map_pre, aes(x = x, y = y))+background_image(img)+xlim(c(1.00000,21.29802))+ylim(c(-4584.6427,292.5073))+geom_point(color = "snow3")+
  theme(legend.position = "none")+
  geom_point(data=map_pre[map_pre$OR == "Olfr160",],color="magenta",size=2)+
  geom_point(data=map_pre[map_pre$OR == "Olfr1507",],color="green",size=2)+
  geom_point(data=map_pre[map_pre$OR == "Olfr73",],color="green",size=2)+
  theme_classic()+NoAxes()
dev.off()

##### 4b #####
pdf("Figure4c.pdf", width = 5, height = 5)
plot_features(map_pre,"class")+NoAxes()+NoLegend()
dev.off()

##### 4c #####
p1 = plot_features(map_pre,"Acsm4")
p2 = plot_features(map_pre,"Nrp2")
p3 = plot_features(map_pre,"Plxna1")
p4 = plot_features(map_pre,"Nrp1")

pdf("Figure4c.pdf", width = 12, height = 10)
wrap_plots(list(p1,p2,p3,p4), ncol = 2)
dev.off()

##### 4d #####
genes = c("Sema3a","Robo2","Mycbp2","Dpysl2","Ntf3","Ccdc141")
genes = genes[genes%in%colnames(map_pre)]
plots = list()
for (i in 1:length(genes)) {
  plots[[i]] = plot_features(map_pre,genes[i])+NoAxes()+NoLegend()+ggtitle(genes[i])
}
pdf("Figure4d.pdf", width = 15, height = 10)
wrap_plots(plots, ncol = 3)
dev.off()

##### 4e #####
Freq = as.data.frame(table(map_pre$OR_cluster))
OR_cluster = as.character(Freq$Var1[Freq$Freq>10]) # only check the greek with more than 10 glomeruli
OR_cluster = OR_cluster[order(nchar(OR_cluster), OR_cluster)]
g2dist_df = map_pre[,c("x","y")]
g2dist_df$x = scales::rescale(g2dist_df$x, to=c(0,3000))
g2dist_df$y = scales::rescale(g2dist_df$y, to=c(0,3000))
rownames(g2dist_df) = map_pre$OR
g2dist_df = dist(g2dist_df)
g2dist_df = melt(as.matrix(g2dist_df))

qdflist = list()
plot_list = list()
for (x in 1:length(OR_cluster)) {
  ORs = map_pre$OR[map_pre$OR_cluster==OR_cluster[x]]
  subdist_df = g2dist_df[(g2dist_df$Var1%in%ORs&g2dist_df$Var2%in%ORs),]
  rank_list = list()
  for (i in 1:length(ORs)) {
    gdist_df2 = subdist_df[subdist_df$Var1==ORs[i],]
    rownames(gdist_df2) = gdist_df2$Var2
    gdist_df2 = gdist_df2[ORs,]
    tmpdf = cbind(gdist_df2$value, map_pre[map_pre$OR%in%ORs,paste("sim",ORs[i], sep = "_")])
    colnames(tmpdf) = c("glomerulus distance","protein sequence similarity")
    rownames(tmpdf) = ORs
    tmpdf = tmpdf[-grep(paste0("^",ORs[i],"$"), rownames(tmpdf)),]
    tmpdf = as.data.frame(tmpdf)
    tmpdf = tmpdf[order(tmpdf$`glomerulus distance`, decreasing = F),]
    rank_list[[i]] = tmpdf$`protein sequence similarity`
  }
  c2df = do.call(rbind, rank_list)
  rownames(c2df) = ORs
  colnames(c2df) = 1:ncol(c2df)
  long_c2df = melt(c2df)
  qdflist[[x]]  = long_c2df %>%
    mutate(quantile = ntile(Var2, 10)) %>% group_by(quantile) %>% summarise_at(vars(value),.funs =  mean)
}
qdf = do.call(rbind, qdflist)
pdf("Figure4e.pdf", width = 4, height = 4)
ggplot(qdf, aes(x = quantile, y = value))+ 
  stat_summary(geom = "point", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, position = "dodge", alpha = 0.5)+geom_smooth(method='glm')+
  xlab("ranked glomeruli position by qunatile (close to distant)")+ylab("Average protein similarity score")+theme_classic()+ scale_x_continuous(breaks=seq(0,10,1))
dev.off()

##### 4f #####
g2dist_df = map_pre[,c("x","y")]
g2dist_df$x = scales::rescale(g2dist_df$x, to=c(0,3000))
g2dist_df$y = scales::rescale(g2dist_df$y, to=c(0,3000))
rownames(g2dist_df) = map_pre$OR
g2dist_df = dist(g2dist_df)
g2dist_df = melt(as.matrix(g2dist_df))

ORs = map_pre$OR
rank_list = list()
for (i in 1:length(ORs)) {
  gdist_df2 = g2dist_df[g2dist_df$Var1==ORs[i],]
  rownames(gdist_df2) = gdist_df2$Var2
  gdist_df2 = gdist_df2[map_pre$OR,]
  tmpdf = cbind(gdist_df2$value, map_pre[,paste("sim",ORs[i], sep = "_")])
  colnames(tmpdf) = c("glomerulus distance","protein sequence similarity")
  rownames(tmpdf) = map_pre$OR
  tmpdf = tmpdf[-grep(paste0("^",ORs[i],"$"), rownames(tmpdf)),]
  tmpdf = as.data.frame(tmpdf)
  tmpdf = tmpdf[order(tmpdf$`glomerulus distance`, decreasing = F),]
  rank_list[[i]] = tmpdf$`protein sequence similarity`
}
c2df = do.call(rbind, rank_list)

rownames(c2df) = ORs
colnames(c2df) = 1:ncol(c2df)

long_c2df = melt(c2df)
pdf("Figure4f.pdf", width = 5, height = 5)
ggplot(long_c2df, aes(x = Var2, y = value))+ 
  stat_summary(geom = "point", fun = mean, position = "dodge")+
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, position = "dodge", alpha = 0.2)+geom_smooth(method = 'glm')+
  xlab("Ranked glomeruli position (min to max)")+ylab("Average protein similarity score")+theme_classic()
dev.off()

##### 4g ######
index = grep("Ligand",colnames(map_pre))
ketones = c(11,12,17,18,19,20)
FA = c(22,37,45,48,54,58)
targets = c(FA,ketones)
cols = hue_pal()(6)
cols2 = c(rep(cols[1],6),rep(cols[3],6))
ligand_plots = list()
for (i in 1:length(targets)) {
  ligand_plots[[i]] = ggplot(map_pre, aes(x, y))+background_image(img)+xlim(c(1.00000,21.29802))+ylim(c(-4584.6427,292.5073))+geom_point(color = "snow3")+
    geom_point(data=map_pre[map_pre[,index[targets[i]]] == "Yes",],color=cols2[i],size=2)+theme_classic()+ggtitle(paste0(gsub("Ligand: ","",colnames(map_pre)[index[targets[i]]])))+NoAxes()
}
pdf("Figure4g.pdf", width = 30, height = 10)
wrap_plots(ligand_plots, ncol = 6)
dev.off()

##### S6a #####
OSN = readRDS("OR_OSN.rds")
coefficients = read.table("data/coefficients.txt")
coefficients_p = read.table("data/coefficient_p_values.txt")
colnames(coefficients_p) = c("axis","PC","p_value")
coefficients_p$coef = c(coefficients$V1,coefficients$V2)
coefficients_p_sig = coefficients_p[coefficients_p$p_value<0.05,]

loading = Loadings(object = OSN[["harmony"]])
loading = abs(loading)
coefficients_p_sig$coef = abs(coefficients_p_sig$coef)

xtmp = list()
x_coefficients_p_sig = coefficients_p_sig[coefficients_p_sig$axis=="x",]
for (i in 1:nrow(x_coefficients_p_sig)) {
  PCs = x_coefficients_p_sig$PC[i]
  xtmp[[i]] = loading[,PCs]*x_coefficients_p_sig$coef[i]
}
ytmp = list()
y_coefficients_p_sig = coefficients_p_sig[coefficients_p_sig$axis=="y",]
for (i in 1:nrow(y_coefficients_p_sig)) {
  PCs = y_coefficients_p_sig$PC[i]
  ytmp[[i]] = loading[,PCs]*y_coefficients_p_sig$coef[i]
}

x_load = as.data.frame(do.call(cbind, xtmp))
rownames(x_load) = rownames(loading)
y_load = as.data.frame(do.call(cbind, ytmp))
rownames(y_load) = rownames(loading)
x_load$summary = rowSums(x_load)
y_load$summary = rowSums(y_load)
x_load = x_load[order(x_load$summary, decreasing = T),]
y_load = y_load[order(y_load$summary, decreasing = T),]
write.table(x_load,"data/gene_contribution_x.txt", col.names = NA, sep = "\t")
write.table(y_load,"data/gene_contribution_y.txt", col.names = NA, sep = "\t")

x_load = read.table("data/gene_contribution_x.txt", header = T, row.names = 1)
y_load = read.table("data/gene_contribution_y.txt", header = T, row.names = 1)

ensembl=useMart("ensembl")
listDatasets(ensembl) 
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
listAttributes(ensembl)
foo <- getBM(attributes=c('entrezgene_id',
                          'mgi_symbol'),mart = ensembl)
target_genes = (rownames(x_load)[!is.na(rownames(x_load))])[1:150]
backgroundgenes = rownames(OSN_filtered)
backgroundid = as.character(foo$entrezgene_id[foo$mgi_symbol%in%backgroundgenes])
targetid = as.character(foo$entrezgene_id[foo$mgi_symbol%in%target_genes])

xBPgenes = enrichGO(gene = targetid, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff=1, pAdjustMethod="BH", universe = backgroundid, qvalueCutoff=1, minGSSize=5)
xBPgenes_df = xBPgenes@result
xBPgenes_df$GeneRatio = as.numeric(sub('\\/.*', '', xBPgenes_df$GeneRatio))/as.numeric(sub('.*\\/', '', xBPgenes_df$GeneRatio))


target_genes = (rownames(y_load)[!is.na(rownames(y_load))])[1:150]
backgroundgenes = rownames(OSN_filtered)
backgroundid = as.character(foo$entrezgene_id[foo$mgi_symbol%in%backgroundgenes])
targetid = as.character(foo$entrezgene_id[foo$mgi_symbol%in%target_genes])

yBPgenes = enrichGO(gene = targetid, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff=1, pAdjustMethod="BH", universe = backgroundid, qvalueCutoff=1, minGSSize=5)
yBPgenes_df = yBPgenes@result
yBPgenes_df$GeneRatio = as.numeric(sub('\\/.*', '', yBPgenes_df$GeneRatio))/as.numeric(sub('.*\\/', '', yBPgenes_df$GeneRatio))

xBPgenes_df = xBPgenes_df[xBPgenes_df$p.adjust<0.05,]
yBPgenes_df = yBPgenes_df[yBPgenes_df$p.adjust<0.05,]

xBPgenes_df = xBPgenes_df[with(xBPgenes_df, order(-p.adjust,GeneRatio, decreasing = T)),]
xBPgenes_df = xBPgenes_df[1:10,]
xBPgenes_df = xBPgenes_df[xBPgenes_df$Count>=10,]
xBPgenes_df$group  = "A-P axis"
yBPgenes_df = yBPgenes_df[with(yBPgenes_df, order(-p.adjust,GeneRatio, decreasing = T)),]
yBPgenes_df = yBPgenes_df[1:10,]
yBPgenes_df = yBPgenes_df[yBPgenes_df$Count>=10,]
yBPgenes_df$group  = "D-V axis"

df = rbind.data.frame(xBPgenes_df,yBPgenes_df)
df$Description = factor(df$Description, levels = rev(unique(sort(df$Description))) )

pdf("FigureS6a.pdf", width = 10, height = 5)
ggplot(df, 
       aes(x = group, y = Description)) + 
  geom_count(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(high = "blue" ,low="red") +
  ylab(NULL) +xlab(NULL)+
  ggtitle("GO biological pathway enrichment")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size_area()
dev.off()

##### S6b #####
source_python("lr.py")
Freq = as.data.frame(table(OSN$OR_identity))
Freq = Freq[Freq$Freq>=7,]

goterms <- as.data.frame(Term(GOTERM))
k = rownames(goterms)[which(goterms$`Term(GOTERM)`=="axon guidance")]
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("SYMBOL"), keytype="GO")
df = df[,c("GO","SYMBOL")]

genes = df$SYMBOL

ORs = as.character(Freq$Var1)
Idents(OSN) = OSN$OR_identity
ORs = as.character(ORs)
OSN654 = subset(OSN, idents = ORs)
subdf = OSN654@meta.data
OSN654 = OSN654@assays$RNA@counts
OSN654 = OSN654[grep("Olfr",rownames(OSN654), invert = T),]
OSN654 = CreateSeuratObject(OSN654)
OSN654 = NormalizeData(OSN654)
genes = unique(genes[genes%in%rownames(OSN654)])
OSN654@assays$RNA@var.features = genes
OSN654 = ScaleData(OSN654, features = genes)
OSN654 = RunPCA(OSN654, npcs = 50) #159
OSN654@meta.data = subdf
OSN654 = RunHarmony(OSN654, group.by.vars = "orig.ident")
dfs = as.data.frame(OSN654@reductions$harmony@cell.embeddings)
dfs$observed = OSN654$OR_identity
result = lr(dfs)
colnames(result) = c("OR","num","type","x","y")
refOR = result$OR[result$type == "slide-seq"]
result = result[result$OR%in%refOR,]
result_obs = result[result$type=="slide-seq",]
result_pre = result[result$type=="predicted",]
result_obs = result_obs[order(result_obs$OR),]
result_pre = result_pre[order(result_pre$OR),]
dfx = data.frame(cbind("predicted" = result_pre$x,"observed"=result_obs$x))
dfy = data.frame(cbind("predicted" = result_pre$y,"observed"=result_obs$y))

p1 = ggplot(dfx,aes(predicted, observed)) + geom_point(color = "black") + 
  geom_smooth(method=lm) + ggtitle(paste("A-P axis","R^2 =",round(R2_Score(dfx$predicted,dfx$observed),2),"RMSE =",round(RMSE(dfx$predicted,dfx$observed),2), sep = " ")) +
  xlab("predicted") + ylab("Observed")+theme_classic()
p2 = ggplot(dfy,aes(predicted, observed)) + geom_point(color = "black") + 
  geom_smooth(method=lm) + ggtitle(paste("A-P axis","R^2 =",round(R2_Score(dfy$predicted,dfy$observed),2),"RMSE =",round(RMSE(dfy$predicted,dfy$observed),2), sep = " ")) +
  xlab("predicted") + ylab("Observed")+theme_classic()

pdf("FigureS6b.pdf", width = 8, height = 4)
p1|p2|p3|p4
dev.off()

##### S6c #####
genes = c("Tubb3","Dclk1","Efna5","Dlx5","Fezf1","Sema6c","Nexn","Etv1")
genes = genes[genes%in%colnames(map_pre)]
plots = list()
for (i in 1:length(genes)) {
  plots[[i]] = plot_features(map_pre,genes[i])+NoAxes()+NoLegend()+ggtitle(genes[i])
}
pdf("FigureS6c.pdf", width = 40, height = 5)
wrap_plots(plots, ncol = 8)
dev.off()

##### S7 #####
OR_cluster = unique(map_pre$OR_cluster)
OR_cluster = OR_cluster[order(nchar(OR_cluster), OR_cluster)]
OR_cluster = OR_cluster[-1]
cols = hue_pal()(length(OR_cluster))
plots = list()
for (i in 1:length(OR_cluster)) {
  plots[[i]] = ggplot(map_pre, aes(x,y))+background_image(img)+xlim(c(1.00000,21.29802))+ylim(c(-4584.6427,292.5073))+ geom_point(color = "snow3")+
    geom_point(data=map_pre[map_pre$OR_cluster == OR_cluster[i],],color=cols[i],size=2)+theme_classic()+NoLegend()+ggtitle(OR_cluster[i])+NoAxes()
}  
pdf("FigureS7.pdf", width = 50, height = 30)
wrap_plots(plots, ncol = 10)
dev.off()

##### S8a #####
map_preII = map_pre[map_pre$class=="II",]
g2dist_df = map_preII[,c("x","y")]
g2dist_df$x = scales::rescale(g2dist_df$x, to=c(0,3000))
g2dist_df$y = scales::rescale(g2dist_df$y, to=c(0,3000))
rownames(g2dist_df) = map_preII$OR
g2dist_df = dist(g2dist_df)
g2dist_df = melt(as.matrix(g2dist_df))

ORs = map_preII$OR
rank_list = list()
for (i in 1:length(ORs)) {
  gdist_df2 = g2dist_df[g2dist_df$Var1==ORs[i],]
  rownames(gdist_df2) = gdist_df2$Var2
  gdist_df2 = gdist_df2[map_preII$OR,]
  tmpdf = cbind(gdist_df2$value, map_preII[,paste("sim",ORs[i], sep = "_")])
  colnames(tmpdf) = c("glomerulus distance","protein sequence similarity")
  rownames(tmpdf) = map_preII$OR
  tmpdf = tmpdf[-grep(paste0("^",ORs[i],"$"), rownames(tmpdf)),]
  tmpdf = as.data.frame(tmpdf)
  tmpdf = tmpdf[order(tmpdf$`glomerulus distance`, decreasing = F),]
  rank_list[[i]] = tmpdf$`protein sequence similarity`
}
c2df = do.call(rbind, rank_list)

rownames(c2df) = ORs
colnames(c2df) = 1:ncol(c2df)

long_c2df = melt(c2df)
long_c2df$Var2 = as.numeric(long_c2df$Var2)
long_c2df$Var3 = log1p(long_c2df$Var2)
pdf("FigureS8a.pdf", width = 5, height = 5)
ggplot(long_c2df, aes(x = Var2, y = value))+ 
  stat_summary(geom = "point", fun = mean, position = "dodge")+
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, position = "dodge", alpha = 0.2)+#geom_smooth(method = 'glm')+
  geom_smooth(method="lm", formula = y~log1p(x))+
  xlab("Ranked glomeruli position (min to max)")+ylab("Average protein similarity score")+theme_classic()
dev.off()

##### S8b #####
index = grep("Ligand",colnames(map_pre))
ketones = c(11,12,17,18,19,20)
FA = c(22,37,45,48,54,58)
targets = c(FA,ketones)
index = index[-targets]
ligand_plots = list()
for (i in 1:length(index)) {
  if (sum(map_pre[,index[i]]=="Yes")<2) {
    next}
  x = length(ligand_plots)
  ligand_plots[[x+1]] = ggplot(map_pre, aes(x, y))+background_image(img)+xlim(c(1.00000,21.29802))+ylim(c(-4584.6427,292.5073))+geom_point(color = "snow3")+
    geom_point(data=map_pre[map_pre[,index[i]] == "Yes",],color= "blue",size=2)+theme_classic()+ggtitle(paste0(gsub("Ligand: ","",colnames(map_pre)[index[i]])))+NoAxes()
}
pdf("FigureS8b.pdf", width = 15, height = 15)
wrap_plots(ligand_plots, ncol = 6)
dev.off()

##### S9a #####
OSN.integrated = readRDS("data/OSN_654int.rds")
map_pre = readRDS("data/map_pre.rds")
index = grep("Ligand",colnames(map_pre))
ketones = c(11,12,17,18,19,20)
FA = c(22,37,45,48,54,58)
ketones_index = index[ketones]
FA_index = index[FA]
FAdf = map_pre[,FA_index]
ketonesdf = map_pre[,ketones_index]
rownames(FAdf) = map_pre$OR
rownames(ketonesdf) = map_pre$OR

ketones_OR = c()
for (i in 1:ncol(ketonesdf)) {
  ketones_OR = c(ketones_OR,rownames(ketonesdf)[which(ketonesdf[,i]=="Yes")])
}
FA_OR = c()
for (i in 1:ncol(FAdf)) {
  FA_OR = c(FA_OR,rownames(FAdf)[which(FAdf[,i]=="Yes")])
}
OR1 = unique(ketones_OR)[!unique(ketones_OR)%in%unique(FA_OR)]
OR2 = unique(FA_OR)[!unique(FA_OR)%in%unique(ketones_OR)]
OSN$ligands = NA
metadf = OSN@meta.data
metadf[metadf$OR_identity%in%OR1,]$ligands = "Ketones"
metadf[metadf$OR_identity%in%OR2,]$ligands = "Fatty acids"
OSN@meta.data = metadf
Idents(OSN) = OSN$ligands
OSN = NormalizeData(OSN, assay = "RNA")
ligand_markers <- FindMarkers(OSN, ident.1 = "Fatty acids", ident.2 = "Ketones", min.pct = 0, only.pos = F, logfc.threshold = 0, test.use="LR", latent.vars = "orig.ident")
ligand_markers$p_val_adj_log = -log10(ligand_markers$p_val_adj)

k = keys(org.Mm.eg.db, keytype="GO")
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("ONTOLOGY"), keytype="GO")
df = df[df$ONTOLOGY=="BP",]
BPGO = df$GO

k = unique(BPGO)
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("SYMBOL"), keytype="GO")
df = df[,c("GO","SYMBOL")]

goterms <- as.data.frame(Term(GOTERM))
tmp = c()
for (i in 1:nrow(df)) {
  tmp[i] = goterms[df$GO[i],]
}
df$GO_term = tmp

genes = rownames(ligand_markers)[ligand_markers$p_val_adj<10e-4&ligand_markers$avg_logFC>.5|ligand_markers$p_val_adj<10e-4&ligand_markers$avg_logFC<(-.5)]
df2 = df[df$SYMBOL%in%genes,]
df2 = df2[df2$GO_term%in%c("lipid metabolic process","oxidation-reduction process"),]
subgenes = c(unique(df2$SYMBOL),"Anxa5","Glyatl3","Acss2") # genes known to be invovled in lipid metabolic process but haven't updated to the GO
pdf("FigureS9a.pdf", width = 8, height = 12)
EnhancedVolcano(ligand_markers,
                lab = rownames(ligand_markers),
                selectLab = subgenes,
                x = 'avg_logFC',
                y = 'p_val_adj',
                title = "Fatty acid vs Ketone",
                subtitle = "",
                pCutoff = 10e-4,
                FCcutoff = .5,
                pointSize = 2,
                labSize = 6,
                drawConnectors = T,
                labCol = 'black',
                boxedLabels = TRUE,
                colAlpha = 1)
dev.off()

##### S9b #####
map_pre2 = map_pre[,(-grep("sim_|Ligand", colnames(map_pre)))]
map_pre2$x = scales::rescale(map_pre2$x, to=c(0,3000))
map_pre2$y = scales::rescale(map_pre2$y, to=c(0,3000))
map_pre3 = as.data.frame(Rotation(map_pre2[,c('x','y')],45*pi/180))
colnames(map_pre3) = c("x","y")
map_pre3 = cbind(map_pre3,map_pre2[,!colnames(map_pre2)%in%c("x","y")])
px = c()
py = c()
cor_x = c()
cor_y = c()
for (i in 10:ncol(map_pre3)) {
  tmpx = cor.test(map_pre3[,"x"],map_pre3[,i])
  cor_x = c(cor_x,tmpx$estimate)
  px = c(px,tmpx$p.value)
  tmpy = cor.test(map_pre3[,"y"],map_pre3[,i])
  cor_y = c(cor_y,tmpy$estimate)
  py = c(py,tmpy$p.value)
}

df = as.data.frame(cbind(colnames(map_pre3)[10:ncol(map_pre3)],cor_x,cor_y,px,py))
df[,2:5] = apply(df[,2:5], 2, as.numeric)
colnames(df) = c("gene", "cor_x", "cor_y", "px", "py")
df$label = df$V1

write.table(df, "data/rotated_position_correlation.txt", sep = "\t")
cor_df = read.table("data/rotated_position_correlation.txt")

cor_df$label = cor_df$gene
cor_df[!(cor_df$cor_x<(-0.6)|
           cor_df$cor_x>0.6|
           cor_df$cor_y<(-0.6)|
           cor_df$cor_y>0.6),]$label = NA
cor_df$color = "blue"
cor_df[!(cor_df$cor_x<(-0.6)|
           cor_df$cor_x>0.6|
           cor_df$cor_y<(-0.6)|
           cor_df$cor_y>0.6),]$color = "gray"
pdf("FigureS9b.pdf", width = 5, height = 5)
ggplot(cor_df, aes(cor_x,cor_y, label = label,color = color))+geom_point()+
  ggrepel::geom_label_repel()+xlab("correlation coefficient to AD-PV axis")+
  ylab("correlation coefficient to AV-PD axis")+theme_classic()+NoLegend()+scale_color_manual(values = c('black','#999999'))
dev.off()

##### S10c #####
OSN = readRDS("data/OR654_rm_all_50SPCs_3000G.rds")
df = OSN@reductions$harmony@cell.embeddings
df_mean = aggregate(df[,1:50], list(OSN$OR_identity), mean)
rownames(df_mean) = df_mean$Group.1
df_mean = df_mean[,-1]
cor_df = cor(t(df_mean))
cor_df = melt(as.matrix(cor_df))
g2dist_df = map_pre[,c("x","y")]
g2dist_df$x = scales::rescale(g2dist_df$x, to=c(0,3000))
g2dist_df$y = scales::rescale(g2dist_df$y, to=c(0,3000))
rownames(g2dist_df) = map_pre$OR
g2dist_df = dist(g2dist_df)
g2dist_df = melt(as.matrix(g2dist_df))

ORs = map_pre$OR
rank_list = list()
for (i in 1:length(ORs)) {
  gdist_df2 = g2dist_df[g2dist_df$Var1==ORs[i],]
  rownames(gdist_df2) = gdist_df2$Var2
  gdist_df2 = gdist_df2[map_pre$OR,]
  cor_df2 = cor_df[cor_df$Var1==paste0(ORs[i]),]
  rownames(cor_df2) = cor_df2$Var2
  cor_df2 = cor_df2[ORs,]
  tmpdf = cbind(gdist_df2$value, cor_df2$value)
  colnames(tmpdf) = c("glomerulus distance","transciptome correlation")
  rownames(tmpdf) = ORs
  tmpdf = tmpdf[-grep(paste0("^",ORs[i],"$"), rownames(tmpdf)),]
  tmpdf = as.data.frame(tmpdf)
  tmpdf = tmpdf[order(tmpdf$`glomerulus distance`, decreasing = F),]
  rank_list[[i]] = tmpdf$`transciptome correlation`
}
c2df = do.call(rbind, rank_list)

rownames(c2df) = ORs
colnames(c2df) = 1:ncol(c2df)
long_c2df = melt(c2df)
pdf("FigureS9c.pdf", width = 5, height = 5)
ggplot(long_c2df, aes(x = Var2, y = value))+ 
  stat_summary(geom = "point", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, position = "dodge", alpha = 0.2)+geom_smooth(method='glm')+
  xlab("Ranked glomeruli distance (min to max)")+ylab("Average transciptome correlation")+theme_classic()
dev.off()

