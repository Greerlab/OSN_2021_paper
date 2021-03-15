#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(Seurat)
library(harmony)
library(org.Mm.eg.db)
library(GO.db)
library(AnnotationDbi)
library(dplyr)

a = as.integer(args[1])
z = as.integer(args[2])


OSN10 = readRDS("OSN10.rds")
Go_ID = readRDS("Go_ID.rds")

dirname = gsub(":","_",Go_ID[a])
dir.create(paste(dirname),showWarnings = F)
all_genes = unique(AnnotationDbi::select(org.Mm.eg.db, keys=Go_ID[a], columns=c("SYMBOL"), keytype="GO")[,"SYMBOL"])
all_genes = all_genes[all_genes%in%rownames(OSN10)]
genes = sample(all_genes, 150)
OSN10 = ScaleData(OSN10, features = genes)
OSN10 = RunPCA(OSN10, npcs = 50)
OSN10 = RunHarmony(OSN10, group.by.vars = "orig.ident")
OSN10 = FindNeighbors(OSN10, reduction = "harmony")
OSN10 = RunSPCA(OSN10, graph = "RNA_snn")
df = as.data.frame(OSN10@reductions$spca@cell.embeddings)
df$observed = OSN10$OR_identity
write.table(df, paste0(dirname,"/SPCA_",z,".csv"), sep = ",")





