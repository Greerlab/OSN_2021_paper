library(reticulate)
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(BioVenn)
library(RColorBrewer)
library(pheatmap)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)
library(reshape2)
library(Seurat)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(caret)
library(scales)
library(harmony)
library(gridExtra)
library(rstatix)
library(reticulate)
library(ggtext)
library(AnnotationDbi)
library(GO.db)
library(viridis)
library(MLmetrics)
library(ComplexHeatmap)
library(jpeg)
library(ape)
library(ggnewscale)
library(ggtree)
library(EnhancedVolcano)


vsd = readRDS("data/vsd.rds") # DEseq2 objecut for bulk RNAseq, transformed
dds = readRDS("data/dds.rds") # DEseq2 objecut for bulk RNAseq
MOE = readRDS("data/MOE.rds") # all cell types detected in MOE
OSN = readRDS("data/OR_OSN.rds") #high confident single OR expressing OSN 
OSN_rm = readRDS("data/OR_OSN_rm.rds") #high confident single OR expressing OSN with OR information removed
Simi = readRDS("data/OR_similarity_score_default.rds") # protein similarity from pairwis alignment
OR_info = read.table("data/OR_cluster_info.txt", header = T) # the file containing the info of OR: chr, genome position, cluster, class...

map_obs = read.table("data/map_obs.txt", header = T) # map observed by slide seq
map_pre = readRDS("data/map_pre.rds") # reconstructed map
# note that the x and y column in the map_pre and map_obs are slide position and pixel
# to transform into um : x*130 (slides are spaced with 130 um) and y*0.65 (beads are spaced with 0.65 um)
img = readJPEG("data/projection3.jpg") # the background of glomeruli map
source("DEseq2_fun.R") # DESseq2 customized function
source("Seruat_fun.R") # Seurat customized function
source("pair_analysis_harmony.R") # for silhouette coefficient analysis
source("plot_features.R") # for ploting the reconstructed map
source_python("Balanced_score.py") # for calculating balanced score

##### for the analysis with bulk RNAseq, MOR28 = Olfr1507, MOR174-9 = Olfr73, M72 = Olfr160 #####
##### 1b #####
plotPCA(vsd , intgroup = "genotype")
PCA_df = plotPCA(vsd , intgroup = "genotype", returnData = T)
p = ggplot(PCA_df,aes(x=PC1,y=PC2,color=group ))+geom_point(size = 5)+
  xlab(paste0("PC1:79% variance"))+ylab(paste0("PC2:15% variance"))+
  theme_bw()+theme(legend.position="bottom")
ggsave(plot=p,height=5.5,width=5,dpi=300, filename="plots/Figure1b.pdf", useDingbats=FALSE)

##### 1c and table S1 #####
res1 <- results(dds, contrast=c("genotype","M72","MOR174_9"))
res2 <- results(dds, contrast=c("genotype","MOR28","MOR174_9"))
res3 <- results(dds, contrast=c("genotype","MOR28","M72"))
res1_sig = get_sig(res1, thresh = 0.01)
res2_sig = get_sig(res2, thresh = 0.01)
res3_sig = get_sig(res3, thresh = 0.01)
set1 = rownames(res1_sig)
set2 = rownames(res2_sig)
set3 = rownames(res3_sig)
pdf("plots/Figure1c.pdf", width = 15, height = 15)
draw.venn(set1,set2,set3, title = "", xtitle = "M72 vs. MOR174_9", ytitle = "MOR28 vs. MOR174_9", ztitle = "MOR28 vs. M72", subtitle = "")
dev.off()

DEG1 = as.data.frame(res1)
DEG1$Genes = rownames(DEG1)
DEG2 = as.data.frame(res2)
DEG2$Genes = rownames(DEG2)
DEG3 = as.data.frame(res3)
DEG3$Genes = rownames(DEG3)
list_of_datasets <- list("M72_M174_9" = DEG1, "MOR28_M174_9" = DEG2, "MOR28_M72" = DEG3)
openxlsx::write.xlsx(list_of_datasets, file = "data/Supplementary table 1.xlsx")

##### 1d #####
res1 <- results(dds, contrast=c("genotype","M72","MOR174_9"))
res2 <- results(dds, contrast=c("genotype","MOR28","MOR174_9"))
res3 <- results(dds, contrast=c("genotype","MOR28","M72"))
res1_sig = get_sig(res1, thresh = 0.01)
res2_sig = get_sig(res2, thresh = 0.01)
res3_sig = get_sig(res3, thresh = 0.01)
genes = unique(c(rownames(res1_sig),rownames(res2_sig),rownames(res3_sig)))
coldata = as.data.frame(colData(dds))[,1,drop =F]

jpeg("plots/Figure1d.jpeg", height = 8, width = 8, unit ="in", res = 300)
cols <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(15))
print(pheatmap(assay(vsd)[genes,], cluster_rows=TRUE, show_rownames=FALSE, color = cols,
               cluster_cols=TRUE, scale="row", border_color=NA, treeheight_row = 0, annotation_col = coldata))
dev.off()

##### 1e #####
#listDatasets(ensembl)
ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
#listAttributes(ensembl)
foo = getBM(attributes=c('entrezgene_id',
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
ggsave(plot=p,height=6,width=6,dpi=300, filename="plots/Figure1e.pdf", useDingbats=FALSE)

##### 1f #####
genes = grep("Sema|Robo|Kirrel|Eph|Efn|Nrp|Plxn", rownames(assay(vsd)), value = T)
pdf("plots/Figure1f.pdf", height = 10, width = 8)
pheatmap(assay(vsd)[genes,], cluster_rows=TRUE, show_rownames=T, color = cols,
         cluster_cols=T, scale="row", border_color=NA, treeheight_row = 0, annotation_col = coldata)
dev.off()

##### E1b #####
pl = list()
genes = c("Olfr73","Olfr160","Olfr1507")
genes = genes[order(nchar(genes), genes)]
for (i in 1:length(genes)) {
  pl[[i]] = plot_gene(dds, gene=genes[i], intgroup="genotype")+ theme(axis.text.x = element_blank(),axis.title.x=element_blank())+ggtitle(paste0(genes[i]))
}
p = wrap_plots(pl,nrow = 1)+ plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.key.size = unit(3,"line"), legend.text = element_text(size=10))
ggsave(plot=p,height=10,width=15,dpi=300, filename="plots/FigureE1b.pdf", useDingbats=FALSE)

##### E1c #####
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
ggsave(plot=p,height=10,width=12,dpi=300, filename="plots/FigureE1c.pdf", useDingbats=FALSE)

##### E2b #####
jpeg("plots/FigureE2b.jpeg", units = "in", width = 9, height = 9,res = 300)
DimPlot(MOE, group.by = "cell_type")+ NoLegend() +NoAxes()
dev.off()

##### E2c #####
marker_genes = c("Krt5","Ascl1","Neurod1","Gap43","Cnga2","Ascl3","Cyp2g1","Muc5b","S100b","Tmem212","C1qb","Ccr2","Igkc","Cd3g","Ngp","H1f5")
plots <- VlnPlot(MOE, features = marker_genes, group.by = "cell_type", pt.size = 0, ncol = 1,combine = F)
plots <- lapply(X = plots, FUN = function(x) x + theme(axis.text.y.left = element_text(size = 7),legend.position = "none"))
pdf("plots/FigureE2c.pdf", width = 12, height = 25)
wrap_plots(plots = plots, ncol = 1)
dev.off()

##### E2d #####
df = readRDS("data/MOE_markers.rds")
top10 <- df %>% purrr::map_df(rev) %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
p = DotPlot(object = MOE, group.by = "cell_type",features = unique(as.character(top10$gene)),dot.scale = 5)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_gradient2(low = "#0571B0", mid = "white", high = "#CA0020", midpoint = 0)
ggsave(plot=p,height=5,width=20,dpi=300, filename="plots/FigureE2d.pdf", useDingbats=FALSE)

##### E3a #####
Idents(MOE) = MOE$cell_type
aOSN = subset(MOE, idents = "OSN")
EM = GetAssayData(aOSN, slot = "counts")
#get rid of none OR expressing OSN
Receptor = grep("^Taar|^Vmn1|^Vmn2|^Fpr|^Ms4a|^Tas1r|^Tas2r", rownames(aOSN), value = T)
REM = EM[Receptor,]
REM[REM==1]=0
REM[REM>0]=1
cells_to_remove = colnames(REM)[colSums(REM)>0]

OR = grep("^Olfr", rownames(aOSN), value = T)
OEM = EM[OR,!colnames(EM)%in%cells_to_remove]
OEM[OEM==1]=0
OEM[OEM>0]=1

ORcount = as.data.frame(colSums(OEM))
ggplot(ORcount, aes(x = `colSums(OEM)`))+geom_histogram(binwidth = 1)+labs(x = "Number of Olfr genes expressed", y = "Cell count")+
  stat_bin(binwidth=1, geom="text", colour="black", size=3.5,
           aes(label=..count..),vjust = -1)+ scale_x_continuous(breaks = seq(1, 15, by = 1))+theme_classic()
ggsave("plots/FigureE3a.pdf", width = 7, height = 6)

##### E3b #####
CleanOSN = OEM[,colSums(OEM)==1] 
EM = GetAssayData(aOSN, slot = "data")
tmp = EM[OR,colnames(CleanOSN)]
tmp = colSums(as.matrix(tmp))
tmp= as.data.frame(tmp)
ggplot(tmp, aes(x = tmp))+geom_histogram(binwidth = 0.01)+labs(title = "Normalized OR expression in all 1 OR OSN", x = "expression")
tmp = tmp[order(tmp$tmp, decreasing = T),,drop=F]
tmp = tmp[1:round(nrow(tmp)*.9,0),,drop = F]# remove the OR OSN with the OR expression at buttom 10%
Clean_OR_OSN_cells = rownames(tmp)#21673
metadf = aOSN@meta.data
metadf = metadf[Clean_OR_OSN_cells,]
nrow(metadf[metadf$nCount_RNA<=15000&metadf$nCount_RNA>=2600,])#20375
ggplot(metadf, aes(x = nCount_RNA))+geom_histogram(binwidth = 150)+
  geom_vline(xintercept = c(2600,15000))+theme_classic()+
  labs(title = "One OR expressing OSN", x = "Total UMI per cell", y = "Cell count")
ggsave("plots/FigureE3b.pdf", width = 7, height = 6)
#This 20375 cells are used for future OSN study, the OR_OSN.rds

##### E3c #####
# reload the OR_OSN.rds
metadf = OSN@meta.data
freq = metadf%>% group_by(orig.ident) %>% count(OR_identity)
freq$OR_identity = with(freq, reorder(OR_identity, n, mean))
freq$label = freq$OR_identity
genes = map_obs$OR
freq[!freq$label%in%genes,]$label = NA

ggplot(freq, aes(OR_identity, n))+stat_summary(geom = "point", fun = mean, alpha = 0.5, color = "red", size = 0.3)+
  stat_summary(geom = "errorbar", fun.data = mean_sd, position = "dodge", alpha = 0.2)+
  stat_summary(aes(label = label),geom = "text_repel", fun = mean, position = "dodge", min.segment.length = 0.1,segment.size = 0.2,
               hjust= -2, vjust =-2,
               box.padding   = 0.35, 
               point.padding = 0.5, nudge_x= -0.1, direction = "y",force_pull = 0,max.overlaps = Inf)+
  theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("OR type")+ylab("Cell number")
ggsave("plots/FigureE3c.pdf", width = 20, height = 10, units = "in")

##### E3d #####
metadf = OSN@meta.data
cell_count = c()
for (i in 1:nrow(OR_info)) {
  cell_count[i] = sum(metadf$OR_identity==OR_info$gene_name[i])
}
OR_info$cells = cell_count
OR_info$gene_name = factor(OR_info$gene_name, levels = OR_info$gene_name)
OR_info$seqnames = as.character(OR_info$seqnames)
OR_info$seqnames = factor(OR_info$seqnames, levels = unique(OR_info$seqnames[order(nchar(OR_info$seqnames), OR_info$seqnames)]))
OR_info$OR_cluster = factor(OR_info$OR_cluster, levels = unique(OR_info$OR_cluster[order(nchar(OR_info$OR_cluster), OR_info$OR_cluster)]))
OR_info$id= 1:nrow(OR_info)

base_data <- OR_info %>% 
  group_by(seqnames) %>% 
  summarize(start=min(id), end=max(id)) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

ggplot(OR_info, aes(gene_name,cells, fill = OR_cluster))+geom_bar(aes(gene_name,cells, fill = OR_cluster),stat='identity')+ylim(-500,400+1)+
  scale_fill_discrete(name = "OR cluster", labels = c(1:68))+theme_void()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  geom_segment(data=OR_info, aes(x = gene_name, y = 400, xend = 0, yend = 400), colour = "grey", alpha=.5, size=0.1 , inherit.aes = FALSE ) +
  geom_segment(data=OR_info, aes(x = gene_name, y = 300, xend = 0, yend = 300), colour = "grey", alpha=.5, size=0.1 , inherit.aes = FALSE ) +
  geom_segment(data=OR_info, aes(x = gene_name, y = 200, xend = 0, yend = 200), colour = "grey", alpha=.5, size=0.1 , inherit.aes = FALSE ) +
  geom_segment(data=OR_info, aes(x = gene_name, y = 100, xend = 0, yend = 100), colour = "grey", alpha=.5, size=0.1 , inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -50, xend = end, yend = -50), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )+
  geom_text(data=base_data, aes(x = title, y = -100, label=seqnames), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)+
  annotate("text", x = 0, y = c(100, 200, 300, 400), label = c("100", "200", "300", "400") , color="grey", size=3 , angle=0, fontface="bold", hjust=1)+
  coord_polar(start = 0) #+theme(plot.background = element_rect(fill = "white"))
ggsave("plots/FigureE3d.pdf", width = 15, height = 12)
#ggsave("plots/FigureE3d.jpeg", width = 15, height = 12, dpi = 300, units = "in")

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
jpeg("plots/Figure2a.jpeg", width = 15, height = 6, res = 300, units = "in")
wrap_plots(plots, ncol = 5)
dev.off()

### following analysis was used to generate the same plot but avoiding the PCs dominated by the most abundant OR OSN types
# reload OR_OSN.rds
Idents(OSN) = OSN$OR_identity
Freq = as.data.frame(table(OSN$OR_identity))
ORs5up = as.character(Freq$Var1[Freq$Freq>5])

OSN5up = subset(OSN, ident = ORs5up)
OSNrest= subset(OSN, ident = ORs5up, invert = T)
metadata = OSN5up@meta.data
#take min number of cells
ncells = 6

# split each OSN group of min number
chunks <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
ncells_list = list()
for (i in 1:length(ORs5up)) {
  sub =  metadata[metadata$OR_identity==ORs5up[i],]
  ncells_list[[i]] = chunks(sample(rownames(sub)), ncells)
}

# average the normalized transcriptome into avergae cells
metadata$group = NA
for (i in 1:ncells) {
  cell = unlist(lapply(ncells_list, `[[`, i))
  metadata$group[rownames(metadata)%in%cell]=i
}
metadata$ave_cell = paste0(metadata$OR_identity,"_",metadata$group)
OSN5up@meta.data = metadata
Idents(OSN5up) = OSN5up$ave_cell
OSN5up_ave = AverageExpression(OSN5up, return.seurat = T)
OSN_balanced = merge(OSN5up_ave, OSNrest)
OSN_balanced = FindVariableFeatures(OSN_balanced, nfeatures = 2000)
genes = VariableFeatures(OSN_balanced)
genes = genes[-grep("Xist|mt-|Rp",genes)]

OSN_balanced_norm = OSN_balanced@assays$RNA@data[genes,]
prin_comp <- prcomp(t(OSN_balanced_norm), center = TRUE,scale. = TRUE)
projected_OSN = scale(t(as.matrix(OSN@assays$RNA@data[genes,])), prin_comp$center, prin_comp$scale) %*% prin_comp$rotation 
OSN[["bpca"]] <- CreateDimReducObject(embeddings = projected_OSN[,1:50], key = "BPC_", assay = DefaultAssay(OSN))
OSN = RunUMAP(OSN, reduction = "bpca", dims = 1:12)

ORs = sample(unique(OSN$OR_identity), 10)
#ORs = c("Olfr728","Olfr741","Olfr1260","Olfr390","Olfr1440","Olfr1030","Olfr571","Olfr1507","Olfr806", "Olfr536")
my_color_palette <- hue_pal()(length(ORs))
Idents(OSN) = OSN$OR_identity
plots = list()
for (i in 1:length(ORs)) {
  cells = WhichCells(OSN, idents = c(paste0(ORs[i])))
  plots[[i]] = DimPlot(OSN, cells.highlight = cells, cols.highlight = my_color_palette[i])+NoAxes()+NoLegend()+ggtitle("")
}
jpeg("plots/Figure2a.jpeg", width = 15, height = 6, res = 300, units = "in")
wrap_plots(plots, ncol = 5)
dev.off()

##### E4a #####
markers = c("Acsm4","Nrp2","Plxna1","Nrp1")
cols = c("gray","yellow","orange","red")
plots = list()
for (i in 1:length(markers)) {
  plots[[i]] = FeaturePlot(OSN, features = markers[i], order = T, cols = cols)+NoAxes()+NoLegend()+ggtitle("")
}
jpeg("plots/FigureE4a.jpeg", width = 6, height = 6, res = 300, units = "in")
wrap_plots(plots, ncol = 2)
dev.off()

##### E4b #####
my_color_palette <- hue_pal()(length(ORs))
Idents(OSN_rm) = OSN_rm$OR_identity
plots = list()
for (i in 1:length(ORs)) {
  cells = WhichCells(OSN_rm, idents = c(paste0(ORs[i])))
  plots[[i]] = DimPlot(OSN_rm, cells.highlight = cells, cols.highlight = my_color_palette[i])+NoAxes()+NoLegend()
}
jpeg("plots/FigureE4b.jpeg", width = 15, height = 6, res = 300, units = "in")
wrap_plots(plots, ncol = 5)
dev.off()

### following analysis was used to generate the same plot but avoiding the PCs dominated by the most abundant OR OSN_rm types
# reload OR_OSN_rm_rm.rds

Idents(OSN_rm) = OSN_rm$OR_identity
Freq = as.data.frame(table(OSN_rm$OR_identity))
ORs5up = as.character(Freq$Var1[Freq$Freq>5])

OSN_rm5up = subset(OSN_rm, ident = ORs5up)
OSN_rmrest= subset(OSN_rm, ident = ORs5up, invert = T)
metadata = OSN_rm5up@meta.data
#take min number of cells
ncells = 6

# split each OSN_rm group of min number
chunks <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
ncells_list = list()
for (i in 1:length(ORs5up)) {
  sub =  metadata[metadata$OR_identity==ORs5up[i],]
  ncells_list[[i]] = chunks(sample(rownames(sub)), ncells)
}

# average the normalized transcriptome into avergae cells
metadata$group = NA
for (i in 1:ncells) {
  cell = unlist(lapply(ncells_list, `[[`, i))
  metadata$group[rownames(metadata)%in%cell]=i
}
metadata$ave_cell = paste0(metadata$OR_identity,"_",metadata$group)
OSN_rm5up@meta.data = metadata
Idents(OSN_rm5up) = OSN_rm5up$ave_cell
OSN_rm5up_ave = AverageExpression(OSN_rm5up, return.seurat = T)
OSN_rm_balanced = merge(OSN_rm5up_ave, OSN_rmrest)
OSN_rm_balanced = FindVariableFeatures(OSN_rm_balanced, nfeatures = 2000)
genes = VariableFeatures(OSN_rm_balanced)
genes = genes[-grep("Xist|mt-|Rp",genes)]

OSN_rm_balanced_norm = OSN_rm_balanced@assays$RNA@data[genes,]
prin_comp <- prcomp(t(OSN_rm_balanced_norm), center = TRUE,scale. = TRUE)
projected_OSN_rm = scale(t(as.matrix(OSN_rm@assays$RNA@data[genes,])), prin_comp$center, prin_comp$scale) %*% prin_comp$rotation 
OSN_rm[["bpca"]] <- CreateDimReducObject(embeddings = projected_OSN_rm[,1:50], key = "BPC_", assay = DefaultAssay(OSN_rm))
OSN_rm = RunUMAP(OSN_rm, reduction = "bpca", dims = 1:25)

ORs = sample(unique(OSN_rm$OR_identity), 10)
#ORs = c("Olfr728","Olfr741","Olfr1260","Olfr390","Olfr1440","Olfr1030","Olfr571","Olfr1507","Olfr806", "Olfr536")
my_color_palette <- hue_pal()(length(ORs))
Idents(OSN_rm) = OSN_rm$OR_identity
plots = list()
for (i in 1:length(ORs)) {
  cells = WhichCells(OSN_rm, idents = c(paste0(ORs[i])))
  plots[[i]] = DimPlot(OSN_rm, cells.highlight = cells, cols.highlight = my_color_palette[i])+NoAxes()+NoLegend()+ggtitle("")
}
jpeg("plots/FigureE4b.jpeg", width = 15, height = 6, res = 300, units = "in")
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
jpeg("plots/Figure2b.jpeg", width = 11, height = 11, units = "in", res = 300) # data might look slightly different since the random selction, but the result and conclusion won't be affected
my_layout <- rbind(c(NA,1:9), c(rep(NA,2),10:17), c(rep(NA,3),18:24), c(rep(NA,4),25:30),c(rep(NA,5),31:35),
                   c(rep(NA,6),36:39),c(rep(NA,7),40:42),c(rep(NA,8),43:44),c(rep(NA,9),45))
grid.arrange(grobs = UMAP_plot, layout_matrix = my_layout)
dev.off()

##### E4c #####
matrix_input = OSN_rm
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
jpeg("plots/FigureE4c.jpeg", width = 11, height = 11, units = "in", res = 300) # data might look slightly different since the random selction, but the result and conclusion won't be affected
my_layout <- rbind(c(NA,1:9), c(rep(NA,2),10:17), c(rep(NA,3),18:24), c(rep(NA,4),25:30),c(rep(NA,5),31:35),
                   c(rep(NA,6),36:39),c(rep(NA,7),40:42),c(rep(NA,8),43:44),c(rep(NA,9),45))
grid.arrange(grobs = UMAP_plot, layout_matrix = my_layout)
dev.off()

##### E4d #####
OSN_rm$OR_identity = as.character(sample(OSN_rm$OR_identity))
matrix_input = OSN_rm
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
jpeg("plots/FigureE4d.jpeg", width = 11, height = 11, units = "in", res = 300) # data might look slightly different since the random selction, but the result and conclusion won't be affected
my_layout <- rbind(c(NA,1:9), c(rep(NA,2),10:17), c(rep(NA,3),18:24), c(rep(NA,4),25:30),c(rep(NA,5),31:35),
                   c(rep(NA,6),36:39),c(rep(NA,7),40:42),c(rep(NA,8),43:44),c(rep(NA,9),45))
grid.arrange(grobs = UMAP_plot, layout_matrix = my_layout)
dev.off()

##### 2c #####
OR654 = as.character(unique(Freq$Var1[Freq$Freq>=7]))
OSN654 = subset(OSN, ident = OR654)
OR654 = OSN654 %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
saveRDS(OR654, "data/OR_OSN_654.rds")

OSN654rm = subset(OSN_rm, ident = OR654)
OSN654rm = OSN654rm %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
saveRDS(OR654rm, "data/OR_OSN_654_rm.rds")

OSN654rmS = OSN654rm
OSN654rmS$OR_identity = sample(OSN654rmS$OR_identity)
saveRDS(OSN654rmS, "data/OR_OSN_654_rms.rds")

# Each of above 3 inputs was used to run the Sil.R on hpcc for parallelizing 213531 combinations from 654 groups of OR
# The result silhouette coefficent was stored as "OR_",start_n,"_",end_n,".rds", start_n and end_n refer to the subset of combinations.
# all the data were downloaded to data/

OR = list.files("data/sil/OR/",pattern = "0.rds|1.rds",full.names = T)
RM = list.files("data/sil/Removed/",pattern = "0.rds|1.rds",full.names = T)
SH = list.files("data/sil/shuffle/",pattern = "0.rds|1.rds",full.names = T)

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
for (i in 1:length(unique(df$V2))) {
  inputs = unique(df$V2)[-i]
  a = df[df$V2==inputs[1],"V1"]
  b = df[df$V2==inputs[2],"V1"]
  print(inputs)
  x = t.test(a,b, alternative = "two.sided")
  print(x)
  print(x$p.value)
}

pdf("plots/Figures2c.pdf", width = 4, height = 6)
ggboxplot(df, x = "V2", y="V1",
          color = "V2", palette = c("#EB4136", "#009044", "#1C75B8"),xlab = " ", ylab = "Silhouette Coefficient")+NoLegend()+
  stat_compare_means(method = "anova", label.y = 1)+      # check global p-value
  scale_x_discrete(labels=c("With OR", "Without OR","Random"))
dev.off()
write.table(df, "data/F2c_result.txt")

##### 2d #####
##left 10 OR
Freq = Freq[order(Freq$Freq, decreasing = T),]
OR10 = as.character(unique(Freq$Var1[1:10]))
OR10 = subset(OSN_rm, ident = OR10)
OR10 = OR10 %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
saveRDS(OR654, "data/OR_10.rds")

# classifier_input_generation_10_OR_hpcc.R was used to generate the train and test split on hpcc
# SVM_10.py was used run the classifer on hpcc and the results were downloaded to data/SVM_10

OR = list.files("data/SVM_10/", pattern = "Result", full.names = T)
Sh = list.files("data/SVM_10_shuffle/", pattern = "Result", full.names = T)

OR_list = list()
for (i in 1:length(OR)) {
  OR_list[[i]] = read.csv(OR[i], row.names = 1)
  OR_list[[i]]$rep = i
}
OR_list = do.call(rbind.data.frame,OR_list)


Sh_list = list()
for (i in 1:length(Sh)) {
  Sh_list[[i]] = read.csv(Sh[i], row.names = 1)
  Sh_list[[i]]$rep = i
}
Sh_list = do.call(rbind.data.frame,Sh_list)
OR_list$group = "OR"
Sh_list$group = "sh"

OR_list$observed = factor(OR_list$observed, levels = unique(OR_list$observed))
OR_list$predicted = factor(OR_list$predicted, levels = unique(OR_list$observed))
Sh_list$observed = factor(Sh_list$observed, levels = unique(Sh_list$observed))
Sh_list$predicted = factor(Sh_list$predicted, levels = unique(Sh_list$observed))

OSNdf = OSN@meta.data

OR_list$UMI = OSNdf[OR_list$barcode,"nCount_RNA"]
Sh_list$UMI = OSNdf[Sh_list$barcode,"nCount_RNA"]

grouped_acc <- function(df, groups) {
  acc_v = c()
  for (i in unique(df[,groups])) {
    tmp = df[df[,groups]==i,]
    acc_v = c(acc_v,mean(balanced_score(tmp$observed,tmp$predicted)))
  }
  return(acc_v)
}

merged = rbind.data.frame(OR_list, Sh_list)

write.table(merged,"data/F2d_10_raw_result.txt")

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

df$Threshold = as.character(df$Threshold)
df$Group = factor(df$Group, levels = c("Shuffle", "OR"))

for (i in 1:length(unique(df$Threshold))) {
  inputs = unique(df$Threshold)[i]
  tmp = df[df$Threshold==inputs,]
  a = tmp[tmp$Group=="OR","Balanced_accuracy"]
  b = tmp[tmp$Group=="Shuffle","Balanced_accuracy"]
  print(inputs)
  x = t.test(a,b, alternative = "two.sided")
  print(x)
  print(x$p.value)
}
write.table(df, "data/F2d_10_result.txt")
pdf("plots/Figure2d_left.pdf", width = 5, height = 7)
ggplot(df, aes(x = Group, y = Balanced_accuracy, fill = Threshold, color = Threshold))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+ylab("Balanced accuracy")+stat_summary(fun=mean, colour="red", geom="text",show.legend = FALSE, aes(x = Group, y = Balanced_accuracy, group = Threshold,label=round(..y.., digits=2)),
                                                    position = position_dodge(width = .75), vjust=-4)+ theme_classic()+
  scale_color_manual(values = c("#EB4136", "#009044", "#1C75B8", "#7F3F94"))+
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = expansion(mult = c(0.15, 0.15)))
dev.off()

for (i in 1:length(unique(df$Threshold))) {
  print(compare_means(Balanced_accuracy ~ Group, data = df[df$Threshold==unique(df$Threshold)[i],], method = "t.test"))
}

##right 654 OR
# classifier_input_generation_654_OR_hpcc.R was used to generate the train and test split on hpcc
# SVM_654.py was used run the classifer on hpcc and the results were downloaded to data/SVM_654

##right 654 OR
OR = list.files("data/SVM_654/", pattern = "Result", full.names = T)
Sh = list.files("data/SVM_654_shuffle/", pattern = "Result", full.names = T)

OR_list = list()
for (i in 1:length(OR)) {
  OR_list[[i]] = read.csv(OR[i], row.names = 1)
  OR_list[[i]]$rep = i
}
OR_list = do.call(rbind.data.frame,OR_list)


Sh_list = list()
for (i in 1:length(Sh)) {
  Sh_list[[i]] = read.csv(Sh[i], row.names = 1)
  Sh_list[[i]]$rep = i
}
Sh_list = do.call(rbind.data.frame,Sh_list)
OR_list$group = "OR"
Sh_list$group = "sh"

OR_list$observed = factor(OR_list$observed, levels = unique(OR_list$observed))
OR_list$predicted = factor(OR_list$predicted, levels = unique(OR_list$observed))
Sh_list$observed = factor(Sh_list$observed, levels = unique(Sh_list$observed))
Sh_list$predicted = factor(Sh_list$predicted, levels = unique(Sh_list$observed))

OSNdf = OSN@meta.data

OR_list$UMI = OSNdf[OR_list$barcode,"nCount_RNA"]
Sh_list$UMI = OSNdf[Sh_list$barcode,"nCount_RNA"]

grouped_acc <- function(df, groups) {
  acc_v = c()
  for (i in unique(df[,groups])) {
    tmp = df[df[,groups]==i,]
    acc_v = c(acc_v,mean(balanced_score(tmp$observed,tmp$predicted)))
  }
  return(acc_v)
}

merged = rbind.data.frame(OR_list, Sh_list)
merged$result = as.numeric(merged$observed==merged$predicted)

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

df$Threshold = as.character(df$Threshold)
df$Group = factor(df$Group, levels = c("Shuffle", "OR"))

for (i in 1:length(unique(df$Threshold))) {
  inputs = unique(df$Threshold)[i]
  tmp = df[df$Threshold==inputs,]
  a = tmp[tmp$Group=="OR","Balanced_accuracy"]
  b = tmp[tmp$Group=="Shuffle","Balanced_accuracy"]
  print(inputs)
  x = t.test(a,b, alternative = "two.sided")
  print(x)
  print(x$p.value)
}
write.table(df, "data/F2d_654_result.txt")
pdf("plots/Figure2d_right.pdf", width = 5, height = 7)
ggplot(df, aes(x = Group, y = Balanced_accuracy, fill = Threshold, color = Threshold))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+ylab("Balanced accuracy")+stat_summary(fun=mean, colour="red", geom="text",show.legend = FALSE, aes(x = Group, y = Balanced_accuracy, group = Threshold,label=round(..y.., digits=2)),
                                                    position = position_dodge(width = .75), vjust=-4)+ theme_classic()+
  scale_color_manual(values = c("#EB4136", "#009044", "#1C75B8", "#7F3F94"))+
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = expansion(mult = c(0.15, 0.15)))
dev.off()

for (i in 1:length(unique(df$Threshold))) {
  print(compare_means(Balanced_accuracy ~ Group, data = df[df$Threshold==unique(df$Threshold)[i],], method = "t.test"))
}





##### E5a and table S2 #####
OSN_10_rm = readRDS("data/OR_OSN_10_rm.rds")
AG_genes = grep("Eph|Nrp|Kirrel|Efn|Pcdh|Sema|Plxn|Robo|Slit|Ntn|Wnt", rownames(OSN_10_rm), value = T)
AG_genes = grep("Pcdha|Pcdhb|Pcdhg", AG_genes, value = T, invert = T)

ORs = as.character(unique(OSN_10_rm$OR_identity))
train.index <- createDataPartition(OSN_10_rm$OR_identity, p = 0.75, list = FALSE)
test_cells = Cells(OSN_10_rm)[-train.index]
train_cells = Cells(OSN_10_rm)[train.index]

train_data = subset(OSN_10_rm, cells = train_cells)
test_data = subset(OSN_10_rm, cells = test_cells)
train_metadata = train_data@meta.data
#take min number of cells
ncells = min(table(train_data$OR_identity))

# split each OSN group of min number
chunks <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
ncells_list = list()
for (i in 1:length(ORs)) {
  sub =  train_metadata[train_metadata$OR_identity==ORs[i],]
  ncells_list[[i]] = chunks(sample(rownames(sub)), ncells)
}

# average the normalized transcriptome into avergae cells
train_metadata$group = NA
for (i in 1:ncells) {
  cell = unlist(lapply(ncells_list, `[[`, i))
  train_metadata$group[rownames(train_metadata)%in%cell]=i
}
train_metadata$ave_cell = paste0(train_metadata$OR_identity,"_",train_metadata$group)

train_data@meta.data = train_metadata
Idents(train_data) = train_data$ave_cell
train_data_ave = AverageExpression(train_data, return.seurat = T)
train_data_ave = FindVariableFeatures(train_data_ave)
train_data_ave$OR_identity = gsub("_.*","",colnames(train_data_ave))

#pca with the averaged cells
train_data_ave_norm = train_data_ave@assays$RNA@data[AG_genes,]
AG_genes = AG_genes[-as.numeric(which(apply(train_data_ave_norm, 1, var) == 0))]
train_data_ave_norm = train_data_ave@assays$RNA@data[AG_genes,]
prin_comp <- prcomp(t(train_data_ave_norm), center = TRUE,scale. = TRUE)
# project all the train data onto the PCA space generated by the averaged train cells
train_data_norm = train_data@assays$RNA@data[AG_genes,]
projected_train_data = scale(t(as.matrix(train_data_norm)), prin_comp$center, prin_comp$scale) %*% prin_comp$rotation 
all(rownames(projected_train_data)==colnames(train_data))
projected_train_data = as.data.frame(projected_train_data)
projected_train_data = projected_train_data[,1:50]
projected_train_data$OR = train_data$OR_identity
saveRDS(prin_comp, paste0("data/pca_AG.rds"))
write.csv(projected_train_data, paste0("data/train_AG.csv"))

# project all the test data onto the PCA space generated by the averaged train cells
test_data_norm = test_data@assays$RNA@data[AG_genes,]
projected_test_data = scale(t(as.matrix(test_data_norm)), prin_comp$center, prin_comp$scale) %*% prin_comp$rotation 
all(rownames(projected_test_data)==colnames(test_data))
projected_test_data = as.data.frame(projected_test_data)
projected_test_data = projected_test_data[,1:50]
projected_test_data$OR = test_data$OR_identity
write.csv(projected_test_data, paste0("data/test_AG.csv"))

#AG_OSN_SVM.ipynb was used to get AG_OSN_coef.csv
pca = readRDS("data/pca_AG.rds")
loading = pca$rotation[,c(1:8)]
SVM_coef = read.csv("data/SVM_AG_coef.csv")
loading = abs(loading)
for (i in 1:8) {
  loading[,i] = loading[,i]*SVM_coef[,2][i]
}
loading = as.data.frame(loading)
loading$contribution = rowSums(loading)
loading = loading[order(loading$contribution,decreasing = T),]
res = loading[,9,drop=F]
res$gene = rownames(res)
openxlsx::write.xlsx(res, file = "data/Supplementary table 2.xlsx", overwrite = T, )


# Get the biological process & Molecular function GO term with at least 150 genes
k = keys(org.Mm.eg.db, keytype="GO")
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("ONTOLOGY"), keytype="GO")
df = df[df$ONTOLOGY=="BP",]
BPGO = df$GO
k = unique(BPGO)
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("SYMBOL"), keytype="GO")
df = df[,c("GO","SYMBOL")]
df = df[df$SYMBOL%in%rownames(OSN),]

x = df %>%
  group_by(GO) %>%
  summarise(count = n_distinct(SYMBOL))
x = x[x$count>=150,]
ref = x$GO
saveRDS(ref, "data/Go_ID_BP.rds")

k = keys(org.Mm.eg.db, keytype="GO")
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("ONTOLOGY"), keytype="GO")
df = df[df$ONTOLOGY=="MF",]
MFGO = df$GO
k = unique(MFGO)
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("SYMBOL"), keytype="GO")
df = df[,c("GO","SYMBOL")]
df = df[df$SYMBOL%in%rownames(OSN),]

x = df %>%
  group_by(GO) %>%
  summarise(count = n_distinct(SYMBOL))
x = x[x$count>=150,]
ref = x$GO
saveRDS(ref, "data/Go_ID_MF.rds")

# get 10 OR groups
# The OSN_10.rds, Go_ID_BP.rds, Go_ID_MF.rds  was used as the input of classifier_input_generation_10_OR_GO_hpcc.R to generate the input of classifier on hpcc,
# Each GO, 100 different input was generated by ramdom sampling the genes in that GO term as the variable genes for PCA
# The SVM_10.py was used to run the prediction result of each GO. the results were downloaded to the data/SVM_10_BP SVM_10_MF
# note that the same analysis was done with the 654 group, the result are the same that Axon guidance genes shown as the top contributor

files = list.files("data/SVM_10_BP", full.names = T)
GO10 = list()
for (i in 1:length(files)) {
  tmp = read.csv(files[i])
  GO = gsub(".*/","",files[i])
  GO_name = strsplit(GO,"_")[[1]][2]
  rep = strsplit(GO,"_")[[1]][4]
  rep = gsub(".csv","",rep)
  tmp = tmp[,c("barcode","observed","predicted")]
  tmp$rep = rep
  acc = grouped_acc(tmp, "rep")
  tmp = c(GO_name,rep,acc)
  GO10[[i]] = tmp
}
GO10 = do.call(rbind.data.frame,GO10)
colnames(GO10) = c("GO","rep","acc")
goterms <- as.data.frame(Term(GOTERM))
rownames(goterms) = gsub("GO:","",rownames(goterms) )
GO10$GO = goterms[GO10$GO,]
GO10$GO = Hmisc::capitalize(as.character(GO10$GO))
GO10$acc = as.numeric(GO10$acc)

files = list.files("data/AG_SVM_10_AM/", full.names = T)
AM10 = list()
for (i in 1:length(files)) {
  tmp = read.csv(files[i])
  tmp = tmp[,c("barcode","observed","predicted")]
  tmp$rep = i
  acc = grouped_acc(tmp, "rep")
  tmp = c(as.character(i),acc)
  AM10[[i]] = tmp
}
AM10 = do.call(rbind.data.frame,AM10)
colnames(AM10) = c("rep","acc")
AM10$GO = "Axon guidance and cell adhesion (curated)"
tmp = rbind.data.frame(GO10,AM10)
tmp$acc = as.numeric(tmp$acc)
pdf("plots/FigureE5a.pdf", width = 20, height = 8)
ggplot(tmp, aes(x = forcats::fct_reorder(GO, acc, .desc = T), y = acc))+geom_boxplot()+
  xlab(NULL)+ylab("Balanced Accuracy")+ylim(c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))+
  theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), legend.position = "none")
dev.off()

bootdata=list()
for (i in 1:length(unique(tmp$GO))) {
  sub_tmp = tmp[tmp$GO==unique(tmp$GO)[i],]
  bootdata[[i]]=sub_tmp[sample(nrow(sub_tmp), 1000, replace=TRUE), ]
}
bootdata = do.call(rbind.data.frame,bootdata)
pdf("plots/FigureE5a.pdf", width = 20, height = 8)
ggplot(bootdata, aes(x = forcats::fct_reorder(GO, acc, .desc = T), y = acc))+geom_boxplot()+
  xlab(NULL)+ylab("Balanced Accuracy")+ylim(c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))+
  theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), legend.position = "none")
dev.off()


##### E5b #####
# Get the biological process GO term with at least 150 genes
#same procedure as E5a but use the genes in MF
files = list.files("data/SVM_10_MF", full.names = T)
GO10 = list()
for (i in 1:length(files)) {
  tmp = read.csv(files[i])
  GO = gsub(".*/","",files[i])
  GO_name = strsplit(GO,"_")[[1]][2]
  rep = strsplit(GO,"_")[[1]][4]
  rep = gsub(".csv","",rep)
  tmp = tmp[,c("barcode","observed","predicted")]
  tmp$rep = rep
  acc = grouped_acc(tmp, "rep")
  tmp = c(GO_name,rep,acc)
  GO10[[i]] = tmp
}
GO10 = do.call(rbind.data.frame,GO10)
colnames(GO10) = c("GO","rep","acc")
goterms <- as.data.frame(Term(GOTERM))
rownames(goterms) = gsub("GO:","",rownames(goterms) )
GO10$GO = goterms[GO10$GO,]
GO10$GO = Hmisc::capitalize(as.character(GO10$GO))
GO10$acc = as.numeric(GO10$acc)

bootdata=list()
for (i in 1:length(unique(GO10$GO))) {
  sub_tmp = GO10[GO10$GO==unique(GO10$GO)[i],]
  bootdata[[i]]=sub_tmp[sample(nrow(sub_tmp), 1000, replace=TRUE), ]
}
bootdata = do.call(rbind.data.frame,bootdata)
pdf("plots/FigureE5b.pdf", width = 20, height = 8)
ggplot(bootdata, aes(x = forcats::fct_reorder(GO, acc, .desc = T), y = acc))+geom_boxplot()+
  xlab(NULL)+ylab("Balanced Accuracy")+ylim(c(0,1))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))+
  theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), legend.position = "none")
dev.off()

##### 2e and E5c #####
OSN = readRDS("data/OR_OSN.rds")
OR_info = read.table("data/OR_cluster_info.txt", header = T)
Simi = readRDS("data/OR_similarity_score_default.rds")
merged = read.table("data/F2d_654_raw_result.txt")
OR_info = OR_info[!duplicated(OR_info$gene_name),]
rownames(OR_info) = OR_info$gene_name
OSN$OR_cluster = OR_info[OSN$OR_identity,"OR_cluster"]
OSN$chr = OR_info[OSN$OR_identity,"seqnames"]
OSN$class = OR_info[OSN$OR_identity,"Class"]
metadf = OSN@meta.data
metadf = metadf[metadf$OR_identity%in%unique(merged$observed),]
Simi = Simi[unique(merged$observed),unique(merged$observed)]
wrong = merged[!(merged$observed==merged$predicted),]

Simi_long = melt(Simi)
rownames(Simi_long) = paste(Simi_long$Var1,Simi_long$Var2, sep = "_")
wrong$test_chr = OR_info[wrong$observed,"seqnames"]
wrong$pred_chr = OR_info[wrong$predicted,"seqnames"]
wrong$test_cls = OR_info[wrong$observed,"OR_cluster"]
wrong$pred_cls = OR_info[wrong$predicted,"OR_cluster"]
wrong$pairs = paste(wrong$observed,wrong$predicted, sep = "_")
wrong$similarity = Simi_long[wrong$pairs,"value"]
wrong$chr = as.numeric(wrong$test_chr==wrong$pred_chr)
wrong$OR_cluster = as.numeric(wrong$test_cls==wrong$pred_cls)
wrong = wrong[!(is.na(wrong$OR_cluster)),]
wrong$cuts = cut(wrong$similarity, seq(0,1,0.1))
wrong$above_0.75 = wrong$similarity>=0.75

write.table(wrong, "data/prediction_error_info.txt", sep = "\t", col.names = NA)
wrong = read.table("data/prediction_error_info.txt", header = T)
wrong$above_0.75 = as.numeric(wrong$similarity>=0.75)

summary_df = wrong %>% group_by(group,rep) %>% summarise(simi = mean(above_0.75), chr2 = mean(chr), cls2 = mean(OR_cluster))

summary_df_s = summary_df[summary_df$group=="sh",]
summary_df_t = summary_df[summary_df$group=="OR",]

fcdf = cbind.data.frame("value" = c(summary_df_t$simi/summary_df_s$simi,summary_df_t$chr2/summary_df_s$chr2,summary_df_t$cls2/summary_df_s$cls2),
                        "group" = c(rep("simi",nrow(summary_df_t)),rep("chr",nrow(summary_df_t)),rep("cls",nrow(summary_df_t))))

pairwise_t_test(data = fcdf, value ~ group, p.adjust.method = "BH")
fcdf$fold_change = log2(as.numeric(fcdf$value))
fcdf$group = factor(fcdf$group,levels = c("chr","cls","simi"))
pdf("plots/Figure2e.pdf", width = 5, height = 7)
ggboxplot(as.data.frame(fcdf), x = "group", y = "fold_change", color = "group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "none")+
  ylab("fold change")+scale_color_manual(values = c("#EB4136", "#009044", "#1C75B8"))
dev.off()

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
summary_df = wrong %>% group_by(group, rep) %>% summarise(Bin0_1 = mean(Bin0_1),
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


fcdf = as.data.frame(cbind(fold_change = c(unlist(bin_list)),
                           group = c(namev)))
fcdf$fold_change = log2(as.numeric(fcdf$fold_change))

for (i in 1:length(unique(fcdf$group))) {
  inputs = unique(fcdf$group)[-i]
  a = fcdf[fcdf$group==inputs[1],"fold_change"]
  b = fcdf[fcdf$group==inputs[2],"fold_change"]
  print(inputs)
  x = t.test(a,b, alternative = "two.sided")
  print(x)
  print(x$p.value)
}

pdf("plots/FigureE5c.pdf", width = 5, height = 7)
ggboxplot(as.data.frame(fcdf), x = "group", y = "fold_change", color = "group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "none")+
  ylab("fold change")
dev.off()



##### E5d #####
#reload Simi
Simi = as.data.frame(as.matrix(Simi))
OSN654 = readRDS("data/OR_OSN_654_rm.rds")
df = as.data.frame(OSN654@reductions$harmony@cell.embeddings)
df$OR = OSN654$OR_identity
write.table(df, "data/OR_OSN_654_50PC_3000g.txt", col.names = NA)
df_mean = aggregate(df[,1:50], list(df$OR), mean)
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

Simi2$cuts = cut(Simi2$protein_sequence_similarity,seq(0,1,0.1))

pdf("plots/FigureE5d.pdf", width = 4, height = 5)
ggplot(Simi2, aes(x=cuts, y=correlation, color = cuts))+
  geom_errorbar(stat="summary", fun.data="mean_se",width=.2)+theme(legend.position = "none")+ylab("transcriptome correlation")+xlab("protein sequence similarity score")+theme_classic()+theme(legend.position = "none")
dev.off()

##### 2f  #####
#reload OSN, Simi, OR_info
OR_info = OR_info[!is.na(OR_info$OR_cluster),]
OR_set1 = unique(OSN$OR_identity)
OR_set2 = unique(OR_info$gene_name)
OR_set3 = unique(rownames(as.matrix(Simi)))
common_OR = intersect(OR_set1,intersect(OR_set2,OR_set3))

OR_info = OR_info[OR_info$gene_name%in%common_OR,]
OR_info = OR_info[!duplicated(OR_info$gene_name),]

Simi2 = Simi[common_OR,common_OR]

Idents(OSN) = OSN$OR_identity
OSN = subset(OSN, idents = common_OR)
metadf = OSN@meta.data
EM = OSN@assays$RNA@counts
EM = EM[-grep("Olfr", rownames(EM)),]
OSN = CreateSeuratObject(EM)
OSN@meta.data = metadf
OSN = NormalizeData(OSN)
OSN = FindVariableFeatures(OSN, nfeatures = 3000)
OSN = ScaleData(OSN)
OSN = RunPCA(OSN)
OSN@meta.data = metadf
OSN = RunHarmony(OSN, group.by.vars = "orig.ident",verbose = F)
df = as.data.frame(OSN@reductions$harmony@cell.embeddings)
df$OR_identity = metadf$OR_identity
write.csv(df, "data/FE5_common_OSN_50SPCs_3000G.csv")
df_mean = aggregate(df[,1:50], list(df$OR_identity), mean)
rownames(df_mean) = df_mean$Group.1
df_mean = df_mean[,-1]
cor_df = cor(t(df_mean))

tss = c()
for (i in 1:nrow(OR_info)) {
  if (OR_info[i,"strand"]=="+") {
    tss[i] = OR_info[i,"start"] 
  } else {tss[i] = OR_info[i,"end"]}
}
OR_info$tss = tss
frequency_list = as.data.frame(table(OR_info$OR_cluster))
target_cluster = as.character(frequency_list$Var1[frequency_list$Freq>=10])

result_list = list()
for (i in 1:length(unique(target_cluster))) {
  ORs = unique(OR_info$gene_name[OR_info$OR_cluster==target_cluster[i]])
  subdf = OR_info[OR_info$gene_name%in%ORs,c("gene_name","tss")]
  mat = as.matrix(dist(subdf$tss))
  colnames(mat) = subdf$gene.name
  rownames(mat) = subdf$gene.name
  longdf1 = melt(mat)
  longOR_info = melt(Simi2[ORs,ORs])
  longdf3 = melt(cor_df[ORs,ORs])
  merged = cbind(longdf1,longOR_info$value,longdf3$value)
  colnames(merged) = c("OR1","OR2","distance","similarity","correlation")
  merged = merged[!(merged$OR1==merged$OR2),]
  merged$cluster = target_cluster[i]
  result_list[[i]] = merged
}
all_result = do.call(rbind.data.frame,result_list)

write.table(all_result,"data/TTS_results.txt", quote = F, sep = "\t")

cor.test(all_result$distance, all_result$similarity)$p.value
cor.test(all_result$distance, all_result$correlation)$p.value
p1 = ggscatter(all_result, x = "distance", y = "similarity", size = 0, alpha = 0)+
  xlab("TSS distance")+ylab("Protein sequence similarity")+
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_viridis()+
  stat_cor(method = "pearson", label.y = 1.1)+ theme(legend.position = "right")+labs(fill = "Counts")

p2 = ggscatter(all_result, x = "distance", y = "correlation", size = 0, alpha = 0)+
  xlab("TSS distance")+ylab("Trancriptome correlation")+
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_viridis()+
  stat_cor(method = "pearson", label.y = 1.1)+ theme(legend.position = "right")+labs(fill = "Counts")

pdf("plots/Figure2f.pdf", width = 10, height = 5)
p1|p2
dev.off()

##### E5f #####
#reload OSN, Simi, OR_info
OR_set1 = unique(OSN_rm$OR_identity)
OR_set2 = unique(OR_info$gene_name)
OR_set3 = unique(rownames(as.matrix(Simi)))
common_OR = intersect(OR_set1,intersect(OR_set2,OR_set3))
df = as.data.frame(OSN_rm@reductions$harmony@cell.embeddings)
df$OR = OSN_rm$OR_identity
df_mean = aggregate(df[,1:50], list(df$OR), mean)
rownames(df_mean) = df_mean$Group.1
df_mean = df_mean[,-1]
cor_df = cor(t(df_mean))
cor_df = as.data.frame(cor_df)

Simi = as.data.frame(as.matrix(Simi))
Simi = Simi[rownames(cor_df),rownames(cor_df)]
all(colnames(Simi) == colnames(cor_df))
all(rownames(Simi) == rownames(cor_df))

OR_info = OR_info[!is.na(OR_info$OR_cluster),]
OR_info = OR_info[OR_info$gene_name%in%common_OR,]
OR_info = OR_info[!duplicated(OR_info$gene_name),]

breaksList = seq((-1), 1, by = 0.01)

OR_info = OR_info[order(OR_info$seqnames,OR_info$start, decreasing = F),]
cor_df2 = cor_df[OR_info$gene_name,OR_info$gene_name]
Simi2 = Simi[OR_info$gene_name,OR_info$gene_name]
cor_df2[lower.tri(cor_df2)] <- Simi2[lower.tri(Simi2)]
all(rownames(cor_df)==df$gene_name)
anndf = OR_info[,c("gene_name","seqnames","OR_cluster")]
rownames(anndf) = anndf$gene_name
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

split_name = OR_info$seqnames
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
pdf("plots/FigureE5f_chr10.pdf", height = 10, width = 12)
Heatmap(as.matrix(OR_ch10), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
dev.off()
jpeg("plots/FigureE5f_chr10.jpeg", height = 10, width = 12, res = 300, units = "in")
Heatmap(as.matrix(OR_ch10), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
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
pdf("plots/FigureE5f_chr2.pdf", height = 10, width = 12)
Heatmap(as.matrix(OR_ch2), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
dev.off()
jpeg("plots/FigureE5f_chr2.jpeg", height = 10, width = 12, res = 300, units = "in")
Heatmap(as.matrix(OR_ch2), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
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
pdf("plots/FigureE5f_chr7.pdf", height = 10, width = 12)
Heatmap(as.matrix(OR_ch7), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
dev.off()
jpeg("plots/FigureE5f_chr7.jpeg", height = 10, width = 12, res = 300, units = "in")
Heatmap(as.matrix(OR_ch7), cluster_columns = F, cluster_rows = F, name = "score",col = inferno(50,direction = 1), show_row_names = F, show_column_names = F, top_annotation = ha,column_title_rot = 90, row_title_rot = 0, column_title_side = "bottom")
dev.off()

##### 3b #####
dt = readRDS("data/OB2_slideseq/Slide11.rds")
genes = c("Kctd12","Calb2", "Cdhr1","Pcp4","Sox11")
ldf = cbind.data.frame("x"=c(1,2,3,4),"y"=c(5,6,7,8))
plots = list()
for (i in 1:length(genes)) {
  plots[[i]] = FeaturePlot(dt, genes[i], order = T, reduction = "Spatial", max.cutoff = "q95", pt.size = 0.5)+
    NoLegend()&scale_color_viridis(option = "C", direction = 1)
}
wrap_plots(plots, nrow = 1)
ggsave("plots/Figure3b.jpeg", width = 27, height = 5, dpi = 300, units = "in")

##### 3c #####
set = c(1,1,1,2,2,2)
slide = c(6,9,10,11,12,15)
ORs = c("Olfr791","Olfr16","Olfr1018","Olfr479","Olfr517","Olfr1250")
plots = list()
for (i in 1:6) {
  dt = readRDS(paste0("data/OB",set[i],"_slideseq/Slide",slide[i],".rds"))
  EM = GetAssayData(dt, "counts")
  OR = ORs[i]
  EM = EM[OR,,drop = F]
  cell_list = colnames(EM)[colSums(as.matrix(EM))>0]
  p1 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+NoLegend()+NoAxes()& 
    scale_color_viridis(option = "D", direction = -1)
  p2 = p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =2)
  plots[[i]] = p2
}
wrap_plots(plots, nrow = 1)
ggsave("plots/Figure3c.jpeg", width = 33, height = 5, dpi = 300, units = "in")

##### 3d #####
df = read.table("data/map_obs.txt", header = T)
pdf("plots/Figure3d.pdf", width = 10, height = 10)
ggplot(df, aes(x = x, y = y, label = OR))+background_image(img)+geom_point()+
  geom_text_repel()+xlab("A---P")+ylab("V---D")+xlim(c(0,2450))+ylim(c(-3130,300))+theme_classic()
dev.off()

##### 3e #####
# reload OR_OSN_rm.rds
OSN_rm_ave = AverageExpression(OSN_rm, return.seurat = T)
OSN_rm_ave = FindVariableFeatures(OSN_rm_ave, nfeatures = 3000)
v3k = VariableFeatures(OSN_rm_ave)
VariableFeatures(OSN_rm) = v3k
OSN_rm = ScaleData(OSN_rm)
OSN_rm = RunPCA(OSN_rm, npcs = 50)
OSN_rm = RunHarmony(OSN_rm, group.by.vars = "orig.ident")
df = as.data.frame(OSN_rm@reductions$harmony@cell.embeddings)
df$observed = OSN_rm$OR_identity
write.table(df, "data/LR_OSN_all_50PCs_3000G.csv", sep = ",")

#LR_OSN_all_50PCs_3000G.csv and LR_LOO.ipynb was used to generate leave one out prediction

dfx = read.csv("data/LR_LOO_x.csv", header = T)
dfy = read.csv("data/LR_LOO_y.csv", header = T)
df = data.frame(cbind(rep(dfx$OR,2),
                      c(dfx$observed,dfx$predicted),
                      c(dfy$observed,dfy$predicted),
                      c(rep("obs",nrow(dfx)),rep("pre",nrow(dfx)))))
colnames(df) = c("OR", "x", "y", "group")
df$x = as.numeric(df$x)
df$y = as.numeric(df$y)




OR_candidates = c("Olfr1217", "Olfr1364" ,"Olfr1030"  , "Olfr47" ,  "Olfr1250"  ,"Olfr598" )
plots = list()
for (i in 1:length(OR_candidates)) {
  plots[[i]] = ggplot(df[df$OR==OR_candidates[i],], aes(x = x, y = y, color = group))+
    background_image(img)+geom_point(size =10)+xlab("A---P")+ylab("V---D")+xlim(c(0,3000))+ylim(c(-3000,400))+theme_classic()+NoLegend()+NoAxes()
}
pdf("plots/Figure3e.pdf", width = 20, height = 30)
wrap_plots(plots, ncol = 2)
dev.off()

##### 3f #####
p1 = ggplot(dfx,aes(predicted, observed)) + geom_point(color = "black") + 
  geom_smooth(method=lm) + ggtitle(paste("A-P axis","R^2 =",round(R2_Score(dfx$predicted,dfx$observed),2),"MAE =",round(MAE(dfx$predicted,dfx$observed),2), sep = " ")) +
  xlab("Predicted") + ylab("Observed")+theme_classic()



p2 = ggplot(dfy,aes(predicted, observed)) + geom_point(color = "black")+
  geom_smooth(method=lm) + ggtitle(paste("D-V axis","R^2 =",round(R2_Score(dfy$predicted,dfy$observed),2),"MAE =",round(MAE(dfy$predicted,dfy$observed),2), sep = " ")) +
  xlab("Predicted") + ylab("Observed")+theme_classic()

pdf("plots/Figure3f.pdf", width = 4, height = 8)
p1 /
  p2
dev.off()

dfp1 = cbind.data.frame(dfx[,2],dfy[,2])
dfp2 = cbind.data.frame(dfx[,3],dfy[,3])
colnames(dfp1) = c("x","y")
colnames(dfp2) = c("x","y")
pdist = c()
for (i in 1:nrow(dfp1)) {
  pdist[i] = dist(rbind.data.frame(dfp1[i,],dfp2[i,]))
}
mean(pdist) #310.8029

##### E6a #####
Sets = c()
Reps = c()
RNA_count = c()
Gene_count = c()
for (x in 1:2) {
  for (i in 1:20) {
    slide = readRDS(paste0("data/OB",x,"_slideseq/Slide",i,".rds"))
    RNA_count = c(RNA_count,slide$nCount_RNA)
    Gene_count = c(Gene_count,slide$nFeature_RNA)
    Sets = c(Sets,rep(x,ncol(slide)))
    Reps = c(Reps,rep(i,ncol(slide)))
  }
}
df = cbind.data.frame(Sets,Reps,RNA_count,Gene_count)
df$Sets = factor(df$Sets)
df$Reps = factor(df$Reps)


write.csv(df, "data/slideseq_qc.csv")
ggplot(df, aes(y = log10(RNA_count),x = Reps,fill = Sets))+geom_boxplot()+theme_classic()+xlab("slides")
ggsave("plots/FigureE6a_umi.pdf", width = 15, height = 10, units = "in")


ggplot(df, aes(y = log10(Gene_count),x = Reps,fill = Sets))+geom_boxplot()+theme_classic()+xlab("slides")
ggsave("plots/FigureE6a_gene.pdf", width = 15, height = 10, units = "in")

##### E6b #####
dt = readRDS("data/OB2_slideseq/Slide11.rds")
jpeg("plots/FigureE6b.jpeg", res = 300, units = "in", width = 5, height = 5)
FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
dev.off()

##### E6c #####
# ss_3_auto_gl_call.R and select_GL.R are used to do the automatical glomeruli calling

##### E6d #####
S1 = read.csv("data/Slideseq_1_glomeruli_position_in_um.csv")
S2 = read.csv("data/Slideseq_2_glomeruli_position_in_um.csv")

duplicated_OR = c(S1$OR,S2$OR)[duplicated(c(S1$OR,S2$OR))]
S1_com = S1[S1$OR%in%duplicated_OR,]
S1_com = S1_com[order(S1_com$OR),]
S2_com = S2[S2$OR%in%duplicated_OR,]
S2_com = S2_com[order(S2_com$OR),]

ORs = S1_com$OR 
mouse1_AP = c()
mouse1_DV = c()
mouse2_AP = c()
mouse2_DV = c()
  
for (o in 1:length(ORs)) {
  S1_com_13 = S1_com[!S1_com$OR%in%ORs[o],]
  S2_com_13 = S2_com[!S2_com$OR%in%ORs[o],]
  P = t(S1_com_13[,c("AP","DV")])
  Q = t(S2_com_13[,c("AP","DV")])
  K <- kabsch(P, Q)
  P = t(S1_com[,c("AP","DV")])
  Q = t(S2_com[,c("AP","DV")])
  new1 = t(K$U %*% P + c(K$R))
  new2 = t(Q)
  new1 = as.data.frame(new1)
  new2 = as.data.frame(new2)
  colnames(new1) = c("AP","DV")
  colnames(new2) = c("AP","DV")
  new1$mouse = 1
  new2$mouse = 2
  new1$OR = ORs
  new2$OR = ORs
  mouse1_AP[o] = new1$AP[new1$OR==ORs[o]]
  mouse2_AP[o] = new2$AP[new2$OR==ORs[o]]
  mouse1_DV[o] = new1$DV[new1$OR==ORs[o]]
  mouse2_DV[o] = new2$DV[new2$OR==ORs[o]]
}

result = cbind.data.frame(ORs,mouse1_AP,mouse1_DV,mouse2_AP,mouse2_DV)




p1 = ggscatter(result, x = "mouse1_AP", y = "mouse2_AP", 
               add = "reg.line", conf.int = TRUE, add.params = list(color = "blue", fill = "lightblue"),
               cor.coef = TRUE, cor.method = "pearson", title = "A-P")
p2 = ggscatter(result, x = "mouse1_DV", y = "mouse2_DV", 
               add = "reg.line", conf.int = TRUE, add.params = list(color = "blue", fill = "lightblue"),
               cor.coef = TRUE, cor.method = "pearson", title = "D-V")

MAE(result$mouse1_AP,mouse2_AP) # 143.1931
MAE(result$mouse1_DV,mouse2_DV) # 110.6756

pdist = c()
for (i in 1:nrow(result)) {
  x = result[i,c(2,3)]
  colnames(x) = c("AP","DV")
  y = result[i,c(4,5)]
  colnames(y) = c("AP","DV")
  pdist[i] = dist(rbind.data.frame(x,y))
}
mean(pdist) #199.3201

pdf("plots/FigureE6d.pdf", width = 10, height = 5)
p1|p2 
dev.off()

##### E6e #####
dt = readRDS("data/OB1_slideseq/Slide8.rds")

EM = GetAssayData(dt, "counts")
OR = "Olfr414"
EM = EM[OR,,drop = F]
cell_list = colnames(EM)[colSums(as.matrix(EM))>0]
p1 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
p2 = p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =2)

dt = readRDS("data/OB1_slideseq/Slide9.rds")
EM = GetAssayData(dt, "counts")
OR = "Olfr16"
EM = EM[OR,,drop = F]
cell_list = colnames(EM)[colSums(as.matrix(EM))>0]
p3 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
p4 = p3+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =2)

jpeg("plots/FigureE6e.jpeg", res = 300, units = "in", width = 10, height = 5)
p2|p4
dev.off()

##### 4a #####
#load map_pre
pdf("plots/Figure4a.pdf", width = 5, height = 5)
ggplot(map_pre, aes(x = x, y = y))+background_image(img)+ xlim(c(0,2450))+ylim(c(-3130,300))+geom_point(color = "snow3")+
  theme(legend.position = "none")+
  geom_point(data=map_pre[map_pre$OR == "Olfr160",],color="magenta",size=2)+
  geom_point(data=map_pre[map_pre$OR == "Olfr1507",],color="green",size=2)+
  geom_point(data=map_pre[map_pre$OR == "Olfr73",],color="green",size=2)+
  theme_classic()+NoAxes()
dev.off()

##### 4b #####
pdf("plots/Figure4b.pdf", width = 5, height = 5)
plot_features(map_pre,"class")+NoAxes()+NoLegend()
dev.off()

##### 4c #####
p1 = plot_features(map_pre,"Acsm4")+NoAxes()+NoLegend()
p2 = plot_features(map_pre,"Nrp2")+NoAxes()+NoLegend()
p3 = plot_features(map_pre,"Plxna1")+NoAxes()+NoLegend()
p4 = plot_features(map_pre,"Nrp1")+NoAxes()+NoLegend()

pdf("plots/Figure4c.pdf", width = 10, height = 10)
wrap_plots(list(p1,p2,p3,p4), ncol = 2)
dev.off()

##### 4d #####
rownames(map_pre) = map_pre$OR
m = -3
c1 = 600
c2 = 1800
c3 = 3000
spots = read.csv("data/MERFISH_spot_summary.csv", row.names = 1)
spots$pair= factor(spots$pair)
spots$type= factor(spots$type)
spots$group= factor(spots$group)
spots$set = factor(spots$set)
ML = spots[spots$ML=="L",]
MR = spots[spots$ML=="R",]
distv = c()
for (i in 1:7) {
  distv[i] = dist(ML[c(i,i+7),c("x","y")])
}
MAE = mean(distv)
cols = hue_pal()(3)
p1 =ggplot(ML, aes(x = x, y = y, color = set))+background_image(img)+
  theme(legend.position = "none")+coord_cartesian(ylim=c(-3130,300), xlim = c(0,2450))+
  theme_classic()+scale_shape_manual(values = c(16,17))+
  geom_abline(intercept = c1, slope = m, color="pink", 
              linetype="dotted", size=1)+
  geom_abline(intercept = c2, slope = m, color="lightgreen", 
              linetype="dotted", size=1)+
  geom_abline(intercept = c3, slope = m, color="lightblue", 
              linetype="dotted", size=1)+
  geom_point(aes(shape = type),size=5)+
  geom_line(aes(group = pair), color = "black", alpha = 0.5, linetype = "dashed")+NoAxes()+
  ggtitle(paste0("MAE = ",round(MAE,2),"um"))+scale_shape_discrete(labels = c("predicted","observed"))+
  scale_color_manual(values=c(cols))+scale_shape_manual(values = c(16,1))

distv = c()
for (i in 1:7) {
  distv[i] = dist(MR[c(i,i+7),c("x","y")])
}
MAE = mean(distv)
cols = hue_pal()(3)
p2 =ggplot(MR, aes(x = x, y = y, color = set))+background_image(img)+
  theme(legend.position = "none")+coord_cartesian(ylim=c(-4584.6427*0.65,292.5073*0.65), xlim = c(1.00000*130,21.29802*130))+
  theme_classic()+scale_shape_manual(values = c(16,17))+
  geom_abline(intercept = c2, slope = m, color="lightgreen", 
              linetype="dotted", size=1)+
  geom_abline(intercept = c3, slope = m, color="lightblue", 
              linetype="dotted", size=1)+
  geom_point(aes(shape = type),size=5)+
  geom_line(aes(group = pair), color = "black", alpha = 0.5, linetype = "dashed")+NoAxes()+
  ggtitle(paste0("MAE = ",round(MAE,2),"um"))+scale_shape_discrete(labels = c("observed"))+
  scale_color_manual(values=c(cols[c(2,3)]))+scale_shape_manual(values = c(16,1))

p1+p2
ggsave("plots/Figure4d.pdf", width = 10, height = 4)

##### 4e #####
genes = c("Sema3a","Robo2","Mycbp2","Dpysl2","Ntf3","Ccdc141")
plots = list()
for (i in 1:length(genes)) {
  plots[[i]] = plot_features(map_pre,genes[i])+NoAxes()+NoLegend()+ggtitle(genes[i])
}
pdf("plots/Figure4e.pdf", width = 15, height = 10)
wrap_plots(plots, ncol = 3)
dev.off()

##### 4f #####
Freq = as.data.frame(table(map_pre$OR_cluster))
OR_cluster = as.character(Freq$Var1[Freq$Freq>10]) # only check the cluster with more than 10 glomeruli
OR_cluster = OR_cluster[order(nchar(OR_cluster), OR_cluster)]
g2dist_df = map_pre[,c("x","y")]
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
pdf("plots/Figure4f.pdf", width = 4, height = 4)
ggplot(qdf, aes(x = quantile, y = value))+ 
  stat_summary(geom = "point", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, position = "dodge", alpha = 0.5)+geom_smooth(method='glm')+
  xlab("ranked glomeruli position by qunatile (close to distant)")+ylab("Average protein similarity score")+theme_classic()+ scale_x_continuous(breaks=seq(0,10,1))
dev.off()

##### 4g #####
g2dist_df = map_pre[,c("x","y")]
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
pdf("plots/Figure4g.pdf", width = 5, height = 5)
ggplot(long_c2df, aes(x = Var2, y = value))+ 
  stat_summary(geom = "point", fun = mean, position = "dodge")+
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, position = "dodge", alpha = 0.2)+geom_smooth(method = 'glm')+
  xlab("Ranked glomeruli position (min to max)")+ylab("Average protein similarity score")+theme_classic()
dev.off()

##### 4h ######
index = grep("Ligand",colnames(map_pre))
ketones = c(11,12,17,18,19,20)
FA = c(22,37,45,48,54,58)
targets = c(FA,ketones)
cols = hue_pal()(6)
cols2 = c(rep(cols[1],6),rep(cols[3],6))
ligand_plots = list()
for (i in 1:length(targets)) {
  ligand_plots[[i]] = ggplot(map_pre, aes(x, y))+background_image(img)+xlim(c(0,2450))+ylim(c(-3130,300))+geom_point(color = "snow3")+
    geom_point(data=map_pre[map_pre[,index[targets[i]]] == "Yes",],color=cols2[i],size=2)+theme_classic()+ggtitle(paste0(gsub("Ligand: ","",colnames(map_pre)[index[targets[i]]])))+NoAxes()
}
pdf("plots/Figure4h.pdf", width = 30, height = 10)
wrap_plots(ligand_plots, ncol = 6)
dev.off()

##### E7c #####
# MERFISH.ipynb was used to generate the plots
##### E8a #####
#reload OR_OSN
#coefficients.txt and coefficient_p_values.txt were generated with ss_6_map_reconstruction.ipynb
coefficients = read.table("data/coefficients.txt")
coefficients_p = read.table("data/coefficient_p_values.txt")
colnames(coefficients_p) = c("axis","PC","p_value")
coefficients_p$coef = c(coefficients$V1,coefficients$V2)
coefficients_p_sig = coefficients_p[coefficients_p$p_value<0.05,]

loading = Loadings(object = OSN_rm[["harmony"]])
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
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
foo <- getBM(attributes=c('entrezgene_id',
                          'mgi_symbol'),mart = ensembl)
target_genes = (rownames(x_load)[!is.na(rownames(x_load))])[1:350]
backgroundgenes = rownames(OSN)
backgroundid = as.character(foo$entrezgene_id[foo$mgi_symbol%in%backgroundgenes])
targetid = as.character(foo$entrezgene_id[foo$mgi_symbol%in%target_genes])

xBPgenes = enrichGO(gene = targetid, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff=1, pAdjustMethod="BH", universe = backgroundid, qvalueCutoff=1, minGSSize=5)
xBPgenes_df = xBPgenes@result
xBPgenes_df$GeneRatio = as.numeric(sub('\\/.*', '', xBPgenes_df$GeneRatio))/as.numeric(sub('.*\\/', '', xBPgenes_df$GeneRatio))


target_genes = (rownames(y_load)[!is.na(rownames(y_load))])[1:350]
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

pdf("plots/FigureE8a.pdf", width = 7, height = 5)
ggplot(df, 
       aes(x = group, y = Description)) + 
  geom_count(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(high = "blue" ,low="red") +
  ylab(NULL) +xlab(NULL)+
  ggtitle("GO biological pathway enrichment")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size_area()
dev.off()

##### E8b #####
genes = c("Tubb3","Dclk1","Efna5","Dlx5","Fezf1","Sema6c","Nexn","Etv1")
plots = list()
for (i in 1:length(genes)) {
  plots[[i]] = plot_features(map_pre,genes[i])+NoAxes()+NoLegend()+ggtitle(genes[i])
}
pdf("plots/FigureE8b.pdf", width = 40, height = 5)
wrap_plots(plots, ncol = 8)
dev.off()

##### E9 #####
OR_cluster = unique(map_pre$OR_cluster)
OR_cluster = OR_cluster[order(nchar(OR_cluster), OR_cluster)]
OR_cluster = OR_cluster[-1]
cols = hue_pal()(length(OR_cluster))
plots = list()
for (i in 1:length(OR_cluster)) {
  plots[[i]] = ggplot(map_pre, aes(x,y))+background_image(img)+xlim(c(0,2450))+ylim(c(-3130,300))+ geom_point(color = "snow3")+
    geom_point(data=map_pre[map_pre$OR_cluster == OR_cluster[i],],color=cols[i],size=2)+theme_classic()+NoLegend()+ggtitle(OR_cluster[i])+NoAxes()
}  
pdf("plots/FigureE9.pdf", width = 50, height = 30)
wrap_plots(plots, ncol = 10)
dev.off()

##### E10a #####
map_preII = map_pre[map_pre$class=="II",]
g2dist_df = map_preII[,c("x","y")]
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
pdf("plots/FigureE10a.pdf", width = 5, height = 5)
ggplot(long_c2df, aes(x = Var2, y = value))+ 
  stat_summary(geom = "point", fun = mean, position = "dodge")+
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, position = "dodge", alpha = 0.2)+
  geom_smooth(method="lm", formula = y~log1p(x))+
  xlab("Ranked glomeruli position (min to max)")+ylab("Average protein similarity score")+theme_classic()
dev.off()

##### E10b #####
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
  ligand_plots[[x+1]] = ggplot(map_pre, aes(x, y))+background_image(img)+xlim(c(0,2450))+ylim(c(-3130,300))+geom_point(color = "snow3")+
    geom_point(data=map_pre[map_pre[,index[i]] == "Yes",],color= "blue",size=2)+theme_classic()+ggtitle(paste0(gsub("Ligand: ","",colnames(map_pre)[index[i]])))+NoAxes()
}
pdf("plots/FigureE10b.pdf", width = 24, height = 12)
wrap_plots(ligand_plots, ncol = 8)
dev.off()

##### E10c #####
#reload OR_OSN.rds annd map_pred.rds
Freq = as.data.frame(table(OSN$OR_identity))
ORs = Freq$Var1[Freq$Freq>=7]
Idents(OSN) = OSN$OR_identity
OSN = subset(OSN, ident = ORs)
metadf = OSN@meta.data
OSN = GetAssayData(OSN,slot = "counts")
OSN = OSN[grep("Olfr",rownames(OSN), invert = T),]
OSN = CreateSeuratObject(OSN)
OSN@meta.data = metadf
OSN = NormalizeData(OSN)
OSN = FindVariableFeatures(OSN, nfeatures = 3000)
VariableFeatures(OSN) = v3k
OSN = ScaleData(OSN)
OSN = RunPCA(OSN, npcs = 50)
OSN = RunHarmony(OSN, group.by.vars = "orig.ident")
OSN = FindNeighbors(OSN, reduction = "harmony")
df = OSN@reductions$harmony@cell.embeddings
df_mean = aggregate(df[,1:50], list(OSN$OR_identity), mean)
rownames(df_mean) = df_mean$Group.1
df_mean = df_mean[,-1]
saveRDS(df_mean,"data/OSN_harmony_agg.rds")
hardist = dist(df_mean)
hc = hclust(hardist, method = "average")
den <- as.phylo(hc)
den$tip.label

rownames(map_pre) = map_pre$OR
colnames(map_pre)[4] = "px"
colnames(map_pre)[5] = "py"
index = grep("Ligand",colnames(map_pre))
ligand_name = grep("Ligand",colnames(map_pre), value = T)
ketones = c(11,12,17,18,19,20)
FA = c(22,37,45,48,54,58)
Ligands_ketones = map_pre[,ligand_name[ketones]]
map_pre$Ketone = unname(ifelse(apply(Ligands_ketones, 1, function(r) any(r %in% c("Yes"))),"Ketone",NA))
Ligands_FA = map_pre[,ligand_name[FA]]
map_pre$`Fatty acid` = unname(ifelse(apply(Ligands_FA, 1, function(r) any(r %in% c("Yes"))),"Fatty acid",NA))
map_pre = map_pre[den$tip.label,]

circ = ggtree(den, layout='fan', branch.length='none', ladderize=F, open.angle = 10, right = T)
gheatmap(circ, map_pre[,c("Nrp1","Nrp2","Plxna1","Nqo1")], offset=.2, width=.2,
         colnames_angle=85,hjust = 1, font.size = 3)+scale_fill_viridis(limits=c(min(map_pre[,c("Nrp1","Nrp2","Plxna1","Nqo1")]), 2))

p1 = gheatmap(circ, map_pre[,c("px"),drop = F], offset=0, width=.1,color = NA,
              colnames_angle=85,hjust = 2, font.size = 3)+scale_fill_viridis(name = "A-P axis", breaks = c(seq(500,2500,500)))

p2 <- p1 + new_scale_fill()
p3 = gheatmap(p2, map_pre[,c("py"),drop = F], offset=2.5, width=.1,color = NA,
         colnames_angle=85,hjust = 2, font.size = 3) +scale_fill_viridis(name = "D-V axis", breaks = c(seq(-500,-2500,-500)))

p4 <- p3 + new_scale_fill()
map_pre[map_pre=="No"]=NA
cols2 = c("Yes" = "Blue", "No" = NA)
gheatmap(p4, map_pre[,c("Ketone","Fatty acid"),drop = F], offset=7, width=.1,color = "gray",
         colnames_angle=85,hjust = 1, font.size = 3)
ggsave("plots/FigureE10c.pdf", width = 12, height = 10)

##### E10d #####
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

genes = rownames(ligand_markers)[ligand_markers$p_val_adj<10e-3&ligand_markers$avg_log2FC>.5|ligand_markers$p_val_adj<10e-3&ligand_markers$avg_log2FC<(-.5)]
df2 = df[df$SYMBOL%in%genes,]
df2 = df2[df2$GO_term%in%c("lipid metabolic process","oxidation-reduction process"),]
subgenes = c(unique(df2$SYMBOL),"Anxa5","Glyatl3","Acss2") # genes known to be invovled in lipid metabolic process but haven't updated to the GO
pdf("plots/FigureE10d.pdf", width = 8, height = 12)
EnhancedVolcano(ligand_markers,
                lab = rownames(ligand_markers),
                selectLab = subgenes,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Fatty acid vs Ketone",
                subtitle = "",
                pCutoff = 10e-3,
                FCcutoff = .5,
                pointSize = 2,
                labSize = 6,
                drawConnectors = T,
                labCol = 'black',
                boxedLabels = TRUE,
                colAlpha = 1)
dev.off()
