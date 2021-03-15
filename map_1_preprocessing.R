library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)
library(Seurat)

OSN = readRDS("data/OR_OSN.rds")
map = read.table("data/map_654.txt")
colnames(map) = c("OR", "cell", "type", "x", "y")

#ref
map_obs = map[map$type =="slide-seq",]
write.table(map_obs, "map_obs.txt", col.names = NA)
map_pre = map[map$type =="predicted",]

#add chr, cluster, class
OR_loci = read.table("data/OR_loci.txt", header = T)
OR_loci = OR_loci[OR_loci$gene.name%in%map_pre$OR,]
OR_loci = OR_loci[!duplicated(OR_loci$gene.name),]
map_pre$OR_cluster = "unknown"
map_pre$chr = "unknown"
map_pre$class = "unknown"
for (i in 1:nrow(OR_loci)) {
  OR = OR_loci$gene.name[i]
  map_pre[map_pre$OR==OR,"OR_cluster"]=OR_loci$OR_cluster[i]
  map_pre[map_pre$OR==OR,"chr"]=OR_loci$chr[i]
  map_pre[map_pre$OR==OR,"chr"]=OR_loci$class[i]
}
map_pre$OR_cluster[is.na(map_pre$OR_cluster)]="unknown"


#add sequence similarity
Simi = readRDS("data/OR_similarity_score_defult.rds")
Simi = Simi[map_pre$OR,map_pre$OR]
Simi = melt(Simi)
OR6 = map_pre$OR

sim_list = list()
for (i in 1:length(OR6)) {
  OR = OR6[i]
  tmp = Simi[Simi$Var1==OR,]
  rownames(tmp) = tmp$Var2
  tmp = tmp[map_pre$OR,]
  sim_list[[i]] = tmp$value
}
sim_list = do.call(cbind,sim_list)
colnames(sim_list) = paste("sim",OR6, sep = "_")
map_pre = cbind(map_pre,sim_list)

#add gene expression
ORs = unique(map_pre$OR)
Idents(OSN) = OSN$OR_identity
OSN = subset(OSN, idents = ORs)
metadf = OSN@meta.data
EM = OSN@assays$RNA@counts
EM = EM[-grep("Olfr", rownames(EM)),]
OSN = CreateSeuratObject(EM)
OSN@meta.data = metadf
OSN.list <- SplitObject(OSN, split.by = "orig.ident")
for (i in 1:length(OSN.list)) {
  OSN.list[[i]] <- SCTransform(OSN.list[[i]], verbose = FALSE)
}
OSN.features <- SelectIntegrationFeatures(object.list = OSN.list, nfeatures = 3000)
OSN.list <- PrepSCTIntegration(object.list = OSN.list, anchor.features = OSN.features, 
                               verbose = FALSE)
OSN.anchors <- FindIntegrationAnchors(object.list = OSN.list, normalization.method = "SCT", 
                                      anchor.features = OSN.features, verbose = FALSE)
OSN.integrated <- IntegrateData(anchorset = OSN.anchors, normalization.method = "SCT", 
                                verbose = FALSE)
OSN.integrated <- RunPCA(OSN.integrated, verbose = FALSE)
saveRDS(OSN.integrated, "data/OSN_654int.rds")
OSN.integrated = readRDS("data/OSN_654int.rds")
Idents(OSN.integrated) =  OSN.integrated$OR_identity
ave  = AverageExpression(OSN.integrated, assays = "integrated", return.seurat = T)
scaled_mtx = ave@assays$integrated@scale.data
scaled_mtx = scaled_mtx[,map_pre$OR]
map_pre = cbind(map_pre,t(scaled_mtx))

#add ligand information
ligand = read.xlsx("data/ligand_Summary.xlsx")
ligand$value = 1

df = reshape(ligand, idvar = "Or", timevar = "Ligand", direction = "wide")
df[is.na(df)] = 0
rownames(df) = df[,1]
df = df[,-1]
colnames(df) = gsub("value.","",colnames(df))
df = df[sort(rownames(df)),sort(colnames(df))]

df_tmp = df
df_tmp = df_tmp[rownames(df_tmp)%in%map_pre$OR,]
df_tmp = df_tmp[colSums(df_tmp)>0]
df_tmp = df_tmp[rowSums(df_tmp)>0] 
df_tmp[rownames(df_tmp)[which(df_tmp$`N-Heptanoic Acid`==1)],"Heptanoic Acid"] = 1
df_tmp[rownames(df_tmp)[which(df_tmp$Ethylisobutyrate==1)],"Ethyl Isobutyrate"] = 1
df_tmp[rownames(df_tmp)[which(df_tmp$Allylphenylacetate==1)],"Allyl Phenylacetate"] = 1
df_tmp[rownames(df_tmp)[which(df_tmp$Prenylacetate==1)],"Prenyl Acetate"] = 1
df_tmp[rownames(df_tmp)[which(df_tmp$`Carvone(-)`==1)],"(-)-Carvone"] = 1
df_tmp[rownames(df_tmp)[which(df_tmp$`Carvone(+)`==1)],"(+)-Carvone"] = 1
df_tmp[rownames(df_tmp)[which(df_tmp$Benzylacetate==1)],"Benzyl Acetate"] = 1
CID_df = read.xlsx("data/Ligand_CID.xlsx") # from Pubchem
df_tmp = df_tmp[,!colnames(df_tmp)%in%c("Allylphenylacetate","Ethylisobutyrate","N-Heptanoic Acid","Prenylacetate","Carvone(-)","Carvone(+)","Benzylacetate")]
colnames(df_tmp)[colnames(df_tmp)%in%CID_df$Chemicals==F]="Hydroxycitronellal Dimethyl Acetal"
write(colnames(df_tmp),"data/chemical_name.txt")


df_tmp[df_tmp==1] ="Yes"
df_tmp[df_tmp==0] ="No"

ORs_in = rownames(df_tmp)
ORs_ex = map_pre$OR[!map_pre$OR%in%ORs_in]

x = matrix(1:(length(ORs_ex)*ncol(df_tmp)), nrow = length(ORs_ex), ncol = ncol(df_tmp))
rownames(x) = ORs_ex
colnames(x) = colnames(df_tmp)
x[x>0]="No"
tmpx = rbind(df_tmp,x)
colnames(tmpx) = gsub("^","Ligand: ",colnames(tmpx))
tmpx = tmpx[map_pre$OR,]
map_pre = cbind(map_pre,tmpx)  #43 OR vs 66 ligands

saveRDS(map_pre,"data/map_pre.rds")