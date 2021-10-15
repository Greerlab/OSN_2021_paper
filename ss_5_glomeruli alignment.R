library(spdep)
library(Seurat)
library(patchwork)
library(pracma)
S1 = read.csv("data/Slideseq_1_glomeruli_position_in_um.csv")
S2 = read.csv("data/Slideseq_2_glomeruli_position_in_um.csv")
S1$mouse = 1
S2$mouse = 2
duplicates = c(S1$OR,S2$OR)[duplicated(c(S1$OR,S2$OR))]
S1_com_tmp = S1[S1$OR%in%duplicates,]
S2_com_tmp = S2[S2$OR%in%duplicates,]
S1_com_tmp = S1_com_tmp[order(S1_com_tmp$OR),]
S2_com_tmp = S2_com_tmp[order(S2_com_tmp$OR),]

P = t(S1_com_tmp[,c("AP","DV")])
Q = t(S2_com_tmp[,c("AP","DV")])
K <- kabsch(P, Q)

P = t(S1[,c("AP","DV")])
Q = t(S2[,c("AP","DV")])

new1 = t(K$U %*% P + c(K$R))
new2 = t(Q)
S1$AP = new1[,1]
S1$DV = new1[,2]
S2$AP = new2[,1]
S2$DV = new2[,2]
Newdf = rbind.data.frame(S1,S2)
Newdf = Newdf[,c("OR","AP","DV","mouse")]

shared_OR = Newdf$OR[duplicated(Newdf$OR)]
Newdftmp = Newdf[!Newdf$OR%in%shared_OR,]
hold = list()
for (x in 1:length(shared_OR)) {
  trg = Newdf[Newdf$OR==shared_OR[x],c("AP","DV")]
  hold[[x]] = c(shared_OR[x],colMeans(trg),"shared")
}
hold = do.call(rbind.data.frame,hold)
colnames(hold) = colnames(Newdftmp)
transformed = rbind.data.frame(Newdftmp,hold)
write.csv(transformed, "data/Aligned_glomeruli.csv")
