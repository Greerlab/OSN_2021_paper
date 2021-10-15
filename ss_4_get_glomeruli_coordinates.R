library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)
df = read.table("data/Glomeruli_candidates.txt")
colnames(df) = c("mouse","slide","OR","dist")
target = df[,c("mouse","slide","OR")]
S1 = df$mouse =1
S2 = df$mouse =2
S1$OR = NA

n = 48

m = target$mouse[n]
s = target$slide[n]
OR = target$OR[n]
dt = readRDS(paste0("data/OB",m,"_slideseq/Slide",s,".rds"))
EM = dt@assays$RNA@counts[OR,]
cell_list = names(which(EM>0))
p1 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+coord_fixed()+ggtitle(paste0(OR))+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =1)
coord = dt@reductions$Spatial@cell.embeddings[cell_list,]
coord

cells = rownames(coord)[-c(1,4)]
cells = rownames(coord)[c(1)]
cells = rownames(coord)[coord[,2]<4300&coord[,2]>4000]
cells = rownames(coord)[coord[,1]>2300]
p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cells,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =1)

S1[S1$slide==s&S1$barcode%in%cells,]$OR = OR
#S1[S1$slide==s&S1$barcode%in%cell_list,]$OR = OR

#project left to right
n = 48

m = target$mouse[n]
s = target$slide[n]
OR = target$OR[n]

tmp = S1[S1$slide==s,]
p1 = ggplot(tmp, aes(x,y, color = percent.mt))+geom_point(size = 0.1)+coord_fixed()+ggtitle(paste0(OR))+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
p1 + geom_point(data = tmp[tmp$OR==OR,,drop = F], aes(x = x, y =y), col = 'magenta', size =1)




p1 + geom_point(data = tmp[tmp$OR==OR,,drop = F], aes(x = x+400, y =y), col = 'magenta', size =1)

S1[S1$slide==s&S1$OR%in%OR,]$x = S1[S1$slide==s&S1$OR%in%OR,]$x+400
S1[S1$slide==s&S1$OR%in%OR,]$y = S1[S1$slide==s&S1$OR%in%OR,]$y-200
# c(8,1000,-200)  c(10,1000,0)  c(20,700,0) c(23,300,0) c(26,700) c(27,1200,0) c(30,500,0) c(34,600,0) c(36,700,0) c(41,400,0) c(47,500,0) c(48,400,0)

S1 = S1[!is.na(S1$OR),]
ORs = unique(S1$OR)
slide = c()
Meanx = c()
Meany = c()
OR = c()
for (i in 1:length(ORs)) {
  tmp = S1[S1$OR==ORs[i],]
  Meanx[i] = mean(tmp$x)
  Meany[i] = mean(tmp$y)
  OR[i] = ORs[i]
  slide[i] = unique(tmp$slide)
}
S1coord = cbind.data.frame(OR,slide,Meanx,Meany)
colnames(S1coord) = c("OR","slide","ML","DV")
write.csv(S1coord, "Slideseq_1_glomeruli_position.csv", row.names = F)


S2$OR = NA

n = 98

m = target$mouse[n]
s = target$slide[n]
OR = target$OR[n]
dt = readRDS(paste0("data/OB",m,"_slideseq/Slide",s,".rds"))
EM = dt@assays$RNA@counts[OR,]
cell_list = names(which(EM>0))
p1 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+coord_fixed()+ggtitle(paste0(OR))+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =1)
coord = dt@reductions$Spatial@cell.embeddings[cell_list,]
coord

cells = rownames(coord)[-c(4)]
cells = rownames(coord)[c(3,4)]
cells = rownames(coord)[coord[,2]<4300&coord[,2]>1650]
cells = rownames(coord)[coord[,1]>3000&coord[,1]<4000&coord[,2]<3000&coord[,2]>2000]
p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cells,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =1)

S2[S2$slide==s&S2$barcode%in%cells,]$OR = OR
#S2[S2$slide==s&S2$barcode%in%cell_list,]$OR = OR

#project left to right
n = 98

m = target$mouse[n]
s = target$slide[n]
OR = target$OR[n]

tmp = S2[S2$slide==s,]
p1 = ggplot(tmp, aes(x,y, color = percent.mt))+geom_point(size = 0.1)+coord_fixed()+ggtitle(paste0(OR))+NoLegend()+NoAxes()& 
  scale_color_viridis(option = "D", direction = -1)
p1 + geom_point(data = tmp[tmp$OR==OR,,drop = F], aes(x = x, y =y), col = 'magenta', size =1)




p1 + geom_point(data = tmp[tmp$OR==OR,,drop = F], aes(x = x+800, y =y), col = 'magenta', size =1)

S2[S2$slide==s&S2$OR%in%OR,]$x = S2[S2$slide==s&S2$OR%in%OR,]$x+800
S2[S2$slide==s&S2$OR%in%OR,]$y = S2[S2$slide==s&S2$OR%in%OR,]$y-200
# c(56,500,0), c(58,700,0), c(67,700,0), c(69,800,0), c(72,600,0), c(82,700,0), c(84,700,0), c(88,500,0), c(92,300,0), c(98,800,0)

S2 = S2[!is.na(S2$OR),]
ORs = unique(S2$OR)
slide = c()
Meanx = c()
Meany = c()
OR = c()
for (i in 1:length(ORs)) {
  tmp = S2[S2$OR==ORs[i],]
  Meanx[i] = mean(tmp$x)
  Meany[i] = mean(tmp$y)
  OR[i] = ORs[i]
  slide[i] = unique(tmp$slide)
}
S2coord = cbind.data.frame(OR,slide,Meanx,Meany)
colnames(S2coord) = c("OR","slide","ML","DV")
write.csv(S2coord, "Slideseq_2_glomeruli_position.csv", row.names = F)
