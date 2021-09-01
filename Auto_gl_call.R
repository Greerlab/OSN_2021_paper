library(Seurat)
library(ggplot2)
library(viridis)
source("select_cells.R")

# select the spot within the glomeruli layer
s = 2
i = 20
dt = readRDS(paste0("data/OB",s,"_slideseq/Slide",i,".rds"))
pd = as.data.frame(dt@reductions$Spatial@cell.embeddings)
colnames(pd) = c("x","y")
pd$mito = dt$percent.mt 
pd$cell = rownames(pd)
pd = pd[order(pd$mito,decreasing = F),]
runGadget(ui, server)

#make a df to call the glomeruli
# 1. at least 2 OR spots
# 2. the number of OR spots outside the glomeruli layer =< the OR spot inside the gloermuli layer
# 3. Kmeans with 2 cluster to split the medial site glomeruli on different OB
# 4. calcualte the mean disatnce to the centronoid
# 5. record the biggest number
# 6. removed the spots if the mean disatnce to the centronoid is more than 170 um
# 7. for the same OR spot detected in the different sections pick the one with smaller mean disatnce to the centronoid
summary_df = list()
for (s in 1:2) {
for (i in 1:20) {
  dt = readRDS(paste0("data/OB",s,"_slideseq/Slide",i,".rds"))
  pd = as.data.frame(dt@reductions$Spatial@cell.embeddings)
  colnames(pd) = c("x","y")
  EM = GetAssayData(dt, "counts")
  
  selected = unique(scan(paste0("OB",s,"_slide",i,"medial_GL.txt"), what = "character"))
  sel = EM[grep("Olfr",rownames(EM)),selected,drop = F]
  sel[sel>0]=1
  sel = sel[rowSums(sel)>=2,,drop = F]
  rowSums(sel)
  
  nsel = EM[rownames(sel),!colnames(EM)%in%selected,drop = F]
  nsel = nsel[,unique(colnames(nsel)),drop = F]
  nsel[nsel>0]=1
  df = cbind.data.frame(rowSums(sel),rowSums(nsel))
  colnames(df)  = c("sel","nsel")
  candidates = rownames(df)[df$sel>=df$nsel]
  
  if (length(candidates) ==0 ) {
    next
  } else {
    OR_detected = c()
    OR_dist2_cen =c()
    for (x in 1:length(candidates)) {
      targets = colnames(sel)[sel[candidates[x],]>0]
      if (length(targets)>2) {
        group = kmeans(pd[targets,c("x","y")],2, nstart = 25)
        centro1 = group$centers[1,]
        points1 = pd[names(group$cluster)[group$cluster==1],c("x","y")]
        dist2cen1 = c()
        for (y in 1:nrow(points1)) {
          dist2cen1[y] = dist(rbind.data.frame(centro1,points1[y,]))
        }
        
        centro2 = group$centers[2,]
        points2 = pd[names(group$cluster)[group$cluster==2],c("x","y")]
        dist2cen2 = c()
        for (y in 1:nrow(points2)) {
          dist2cen2[y] = dist(rbind.data.frame(centro2,points2[y,]))
        }
        print(paste(candidates[x],round(max(c(mean(dist2cen1),mean(dist2cen2)))*0.65,2)))
        OR_detected[x] = candidates[x]
        OR_dist2_cen[x] = round(max(c(mean(dist2cen1),mean(dist2cen2)))*0.65,2)
      } else {
        centro = colMeans(pd[targets,c("x","y")])
        points = pd[targets,c("x","y")]
        dist2cen = c()
        for (y in 1:nrow(points)) {
          dist2cen[y] = dist(rbind.data.frame(centro,points[y,]))
        }
        print(paste(candidates[x],round(mean(dist2cen)*0.65,2)))
        OR_detected[x] = candidates[x]
        OR_dist2_cen[x] = round(mean(dist2cen)*0.65,2)
      }
    }
  }
  summary_df[[(s-1)*20+i]] =  cbind.data.frame(rep(s,length(OR_detected)),rep(i,length(OR_detected)),OR_detected,OR_dist2_cen)
}
  
}
summary_df = do.call(rbind.data.frame,summary_df)
colnames(summary_df) = c("set","slide","OR","dist")
summary_df = summary_df[summary_df$dist<170,]

duplicated(summary_df$OR[summary_df$set==1]) # 48 OR for set1, no duplicated
duplicated(summary_df$OR[summary_df$set==2]) # 52 OR for set2 2 duplicated
summary_df$OR[summary_df$set==2][duplicated(summary_df$OR[summary_df$set==2])]
#set 2 Olfr987 and Olfr772 are duplicated
summary_df = summary_df[order(summary_df$dist, decreasing = F),]
summary_df = summary_df[!duplicated(summary_df[,c("set","OR")]),]
summary_df = summary_df[ order(as.numeric(rownames(summary_df))), ]
write.table(summary_df, "data/Glomeruli_candidates.txt")


df = read.table("data/Glomeruli_candidates.txt")