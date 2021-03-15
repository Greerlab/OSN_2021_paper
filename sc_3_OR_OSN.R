library(Seurat)
library(patchwork)
library(ggplot2)
library(scales)
library(harmony)
library(viridis)
library(monocle3)
cols = c("gray","yellow","orange","red")
OSN = readRDS("data/OSN.rds")

#remove OSN expressing non-canonical receptor
EM = GetAssayData(OSN, slot = "counts")
Receptor = grep("^Olfr|^Taar|^Vmn1|^Vmn2|^Fpr|^Ms4a|^Tas1r|^Tas2r", rownames(EM), value = T)
REM = EM[Receptor,]
REM[REM==1]=0
REM[REM>0]=1
CleanOSN = REM[,colSums(REM)==1] 
ncol(CleanOSN)# 24231

#Get the single OR expressing OSN population
OR = grep("^Olfr", rownames(EM), value = T)
OREM = EM[OR,colnames(CleanOSN)]
OREM[OREM==1]=0
OREM[OREM>1]=1
Clean_OR_OSN_cells = colnames(OREM[,colSums(OREM)==1]) # 24081
EM = GetAssayData(OSN, slot = "data")
tmp = EM[OR,Clean_OR_OSN_cells]
tmp = colSums(as.matrix(tmp))
tmp= as.data.frame(tmp)
ggplot(tmp, aes(x = tmp))+geom_histogram(binwidth = 0.01)+labs(title = "Normalized OR expression in all 1 OR OSN", x = "expression")


tmp = tmp[order(tmp$tmp, decreasing = T),,drop=F]
tmp = tmp[1:round(nrow(tmp)*.9,0),,drop = F]
min(tmp$tmp)# 3.320634
Clean_OR_OSN_cells = rownames(tmp)#21672

Clean_OR_OSN = subset(OSN, cells = Clean_OR_OSN_cells)
Clean_OR_OSN = GetAssayData(Clean_OR_OSN, slot = "counts")
Clean_OR_OSN = CreateSeuratObject(Clean_OR_OSN)
Clean_OR_OSN <- NormalizeData(Clean_OR_OSN)
Clean_OR_OSN <- FindVariableFeatures(Clean_OR_OSN)
Clean_OR_OSN <- ScaleData(Clean_OR_OSN)
Clean_OR_OSN <- RunPCA(Clean_OR_OSN)
Clean_OR_OSN = RunHarmony(Clean_OR_OSN, group.by.vars = "orig.ident",verbose = F)
EM = GetAssayData(Clean_OR_OSN, slot = "counts")
EM = EM[OR,]
EM[EM==1]=0
EM[EM>1]=1

OR_vac = c()
for (i in 1:ncol(EM)) {
  tmp2 = EM[,i]
  OR_vac[i] = names(which(tmp2>0))
}

Clean_OR_OSN$OR_identity = OR_vac

#Increase the UMI threshold to keep the high confidence OSN
ggplot(OSN@meta.data, aes(x = nCount_RNA))+geom_histogram(binwidth = 50)+geom_vline(xintercept = c(2600,15000))
Clean_OR_OSN = subset(Clean_OR_OSN, nCount_RNA>=2600 & nCount_RNA<=15000)
metadf = Clean_OR_OSN@meta.data
Clean_OR_OSN = CreateSeuratObject(Clean_OR_OSN@assays$RNA@counts)
Clean_OR_OSN = NormalizeData(Clean_OR_OSN)
Clean_OR_OSN = FindVariableFeatures(Clean_OR_OSN)
Clean_OR_OSN = ScaleData(Clean_OR_OSN)
Clean_OR_OSN = RunPCA(Clean_OR_OSN)
ElbowPlot(Clean_OR_OSN)
Clean_OR_OSN@meta.data = metadf
Clean_OR_OSN = RunHarmony(Clean_OR_OSN, group.by.vars = "orig.ident")
Clean_OR_OSN <- RunUMAP(Clean_OR_OSN, reduction = "harmony", dims = 1:18)
DimPlot(Clean_OR_OSN)
saveRDS(Clean_OR_OSN, "data/OR_OSN.rds")
