#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(harmony)
library(lsa)
library(cluster)

start_n = as.integer(args[1])
end_n = as.integer(args[2])

OSN = readRDS("OR_OSN_654.rds") # change to the input if need it and do it in different folder
Freq = as.data.frame(table(OSN$OR_identity))
ORs = as.character(Freq$Var1[Freq$Freq>=7])
matrix_input = OSN
Ident_input = "OR_identity"
pair_input = ORs
nGenes = 2000
ncell = NULL

Idents(matrix_input) = matrix_input@meta.data[,Ident_input]
Pairs = combn(unique(pair_input),2)
Silhouette = c()

for (i in start_n:end_n) {
  
  target_cell_A = WhichCells(matrix_input, idents = Pairs[1,i])
  target_cell_B = WhichCells(matrix_input, idents = Pairs[2,i])
  if (is.null(ncell)) {
    ncells = min(length(target_cell_A),length(target_cell_B))
  }
  Pair_mtx = subset(matrix_input, cells = c(sample(target_cell_A,ncells),sample(target_cell_B,ncells)))
  Pair_mtx = NormalizeData(Pair_mtx, verbose = FALSE)
  Pair_mtx = FindVariableFeatures(Pair_mtx, verbose = FALSE, nfeatures = nGenes)
  Pair_mtx = ScaleData(Pair_mtx, verbose = FALSE)
  npcs = 50
  if (ncol(Pair_mtx)<51) {
    npcs = ncol(Pair_mtx)-1
  }
  Pair_mtx <- RunPCA(Pair_mtx, verbose = FALSE,npcs = npcs)
  if (ncol(Pair_mtx)<20) {
    center = ncol(Pair_mtx)-1
    Pair_mtx <- RunHarmony(Pair_mtx, group.by.vars = "orig.ident", verbose = FALSE, nclust = center)
  } else {
    Pair_mtx <- RunHarmony(Pair_mtx, group.by.vars = "orig.ident", verbose = FALSE)
  }
  df = as.data.frame(Pair_mtx@reductions$harmony@cell.embeddings[,1:2])
  df2 = Pair_mtx@meta.data
  df2$cluster = NA
  df2[df2[,Ident_input]==Pairs[1,i],]$cluster=1
  df2[df2[,Ident_input]==Pairs[2,i],]$cluster=2
  Dist = cosine(t(df))
  Dist = 1-Dist
  Sil = silhouette(df2$cluster, Dist)
  res = round(summary(Sil)[[4]],2)
  Silhouette = c(Silhouette,res)
}

saveRDS(Silhouette, paste0("OR_",start_n,"_",end_n,".rds"))
