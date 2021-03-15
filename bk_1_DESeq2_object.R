library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(patchwork)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Rgraphviz)
library(gridExtra)
library(ggpubr)
library(openxlsx)

source("DEseq2_fun.R")
mtx = read.table("data/genes_expression_expected_count.tsv", sep = "\t", header = T, row.names = 1)
mtx = mtx[,-1]
df2 <- data.frame(apply(mtx, 2, function(x) round(x,0)))
coldata = as.data.frame(cbind(genotype = c(rep("M72",3),rep("MOR174_9",3),rep("MOR28",3)),
                              batch = c(1:9)))
rownames(coldata) = colnames(mtx)

##DESeq object
dds <- DESeqDataSetFromMatrix(countData = df2,
                              colData = coldata,
                              design= ~ genotype)
##filtering
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]
dds <- DESeq(dds)
#Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)

saveRDS(dds, file = "data/dds.rds")
saveRDS(vsd, file = "data/vsd.rds")
