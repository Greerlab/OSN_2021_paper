library(Seurat)
library(ggplot2)
library(spdep)
library(pracma)
source("select_cells.R")
m = 2
s = 20
dt = readRDS(paste0("data/OB",m,"_slideseq/Slide",s,".rds"))
pd = as.data.frame(dt@reductions$Spatial@cell.embeddings)
colnames(pd) = c("x","y")
pd$mito = dt$percent.mt
pd$cell = rownames(pd)
pd = pd[order(pd$mito,decreasing = F),]
runGadget(ui, server)


m = 2
s = 20
dt = readRDS(paste0("data/OB",m,"_slideseq/Slide",s,".rds"))
bc = scan(paste0("data/OB",m,"_slideseq/Slide",s,"clean.txt"), "character")
df = dt@meta.data[bc,]
df$x = dt@reductions$Spatial@cell.embeddings[bc,1]
df$y = dt@reductions$Spatial@cell.embeddings[bc,2]

ggplot(df, aes(x, y, color = percent.mt)) +
  geom_point(size = 3, shape=20, alpha=1) +
  scale_color_viridis(option = "D", direction = -1) +
  theme(legend.position = "none")+geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+ coord_fixed()

#rotate 
df2 = as.data.frame(Rotation(df[,c('x','y')],-105*pi/180))
df2$percent.mt = df$percent.mt

p = ggplot(df2, aes(V1, V2, color = percent.mt)) +
  geom_point(size = 0.5, shape=20, alpha=1) +
  scale_color_viridis(option = "D", direction = -1) +geom_vline(xintercept = 0, color = "Red")+
  geom_hline(yintercept = 0, color = "Blue")+ coord_fixed()

toWebGL(p)

#translation
xshift = 00
yshift = 7500
df3 = df2
df3$V1 = df2$V1+xshift
df3$V2 = df2$V2+yshift
df3$label = rownames(df3)
p = ggplot(df3, aes(V1, V2, color = percent.mt, label = label)) +
  geom_point(size = 0.5, shape=20, alpha=1) +
  scale_color_viridis(option = "D", direction = -1) +geom_vline(xintercept = 0, color = "Red")+
  geom_hline(yintercept = 0, color = "Blue")+ coord_fixed()
toWebGL(p)


df$x = df3$V1
df$y = df3$V2
df = df[,c("x","y","percent.mt")]
write.csv(df, paste0("data/OB",m,"_slideseq/Slide",s,"_rotated.csv"))


## anchors
source("select_anchors.R")
m = 2
s = 20
df = read.csv(paste0("data/OB",m,"_slideseq/Slide",s,"_rotated.csv"))
df$anchor = NA
shinyApp(ui, server)

anchors = c("TCCTAAATGGCCAA","CCGAACATCCAAAC","CTGGATCACTAATA")
df[df$X%in%anchors,]$anchor = "anchors"

p = ggplot(df, aes(x, y, color = percent.mt, label = X)) +
  geom_point(size = 0.5, shape=20, alpha=1) +
  scale_color_viridis(option = "D", direction = -1) +
  coord_fixed()+geom_point(data = df[!is.na(df$anchor),] ,size = 5, shape=20, alpha=1, color = "magenta")
toWebGL(p)

write.csv(df, paste0("data/OB",m,"_slideseq/Slide",s,"_rotated.csv"))


df = list()
for (s in 1:20) {
  m = 1
  df[[s]] = read.csv(paste0("data/OB",m,"_slideseq/Slide",s,"_rotated.csv"))
  df[[s]]$slide = s
}
df = do.call(rbind.data.frame, df)

df2 = list()
for (s in 1:20) {
  m = 2
  df2[[s]] = read.csv(paste0("data/OB",m,"_slideseq/Slide",s,"_rotated.csv"))
  df2[[s]]$slide = s
}
df2 = do.call(rbind.data.frame, df2)

df2$x = df2$x+4000
df$x = df$x+500
df_all = rbind.data.frame(df,df2)

plots = list()
for (i in 1:20) {
  sub = df_all[df_all$slide==i,]
  plots[[i]] = ggplot(sub, aes(x, y, color = percent.mt, label = X)) +
    geom_point(size = 0.5, shape=20, alpha=1) +
    scale_color_viridis(option = "D", direction = -1)+
    scale_y_continuous(minor_breaks = seq(1000, 9000, 100), breaks = seq(1000, 9000, 500), limits = c(1000,9000))+
    scale_x_continuous(minor_breaks = seq(0, 8000, 100), breaks = seq(0, 8000, 500), limits = c(0,8000))+
    geom_point(data = sub[!is.na(sub$anchor),] ,size = 5, shape=20, alpha=1, color = "magenta")+ coord_fixed()+ theme(legend.position = "none")
}

library(patchwork)
wrap_plots(plots, ncol = 5)
ggsave("plots/OB.jpeg", width = 40, height = 20, units = "in")

#### OB1
df$z = df$slide*130*100/65
tmp = df[!is.na(df$anchor),]
tmp_list = list()
for (i in 1:20) {
  tmp_tmp = tmp[tmp$slide==i,]
  tmp_tmp$anchor = "Mid"
  tmp_tmp[which.max(tmp_tmp$y),"anchor"]="Top"
  tmp_tmp[which.min(tmp_tmp$y),"anchor"]="But"
  tmp_list[[i]] = tmp_tmp
}
tmp = do.call(rbind.data.frame,tmp_list)

plot_ly(tmp, x = ~z, y = ~y, color = ~anchor, marker = list(size = 10))
tmp2 = tmp
slide = 20
tmp2[tmp2$slide==slide,]$y = tmp2[tmp2$slide==slide,]$y+200
plot_ly(tmp2, x = ~z, y = ~y, color = ~anchor, marker = list(size = 10))

tmp3 = tmp2
slide = 19
tmp3[tmp3$slide==slide,]$x = tmp3[tmp3$slide==slide,]$x-100
plot_ly(tmp3, x = ~x, y = ~slide, color = ~anchor, marker = list(size = 10))


for (i in 1:20) {
  xd = unique(tmp3$x[tmp3$slide==i]-tmp$x[tmp$slide==i])[1]
  yd = unique(tmp3$y[tmp3$slide==i]-tmp$y[tmp$slide==i])[1]
  df[df$slide==i,]$x = df[df$slide==i,]$x+xd
  df[df$slide==i,]$y = df[df$slide==i,]$y+yd
  
}
df$x = df$x+200
df$y = df$y-1500
df[df$slide==15,]$y = df[df$slide==15,]$y-300
plots = list()
for (i in 1:20) {
  sub = df[df$slide==i,]
  plots[[i]] = ggplot(sub, aes(x, y, color = percent.mt, label = X)) +
    geom_point(size = 0.5, shape=20, alpha=1) +
    scale_color_viridis(option = "D", direction = -1)+
    scale_y_continuous(minor_breaks = seq(0, 8000, 100), breaks = seq(0, 8000, 500), limits = c(0,8000))+
    scale_x_continuous(minor_breaks = seq(0, 4500, 100), breaks = seq(0, 4500, 500), limits = c(0,4500))+
    geom_point(data = sub[!is.na(sub$anchor),] ,size = 5, shape=20, alpha=1, color = "magenta")+ coord_fixed()+ theme(legend.position = "none")
}

plot_ly(df, x = ~x, y = ~y, z = ~z, color = ~percent.mt, marker = list(size = 0.5))

write.table(df, "data/OB1_slideseq/All_slides_3d.txt")

library(patchwork)
wrap_plots(plots, ncol = 5)
ggsave("plots/OB1.jpeg", width = 25, height = 20, units = "in")

#### OB2
df2$z = df2$slide*130*100/65
tmp = df2[!is.na(df2$anchor),]
tmp_list = list()
for (i in 1:20) {
  tmp_tmp = tmp[tmp$slide==i,]
  tmp_tmp$anchor = "Mid"
  tmp_tmp[which.max(tmp_tmp$y),"anchor"]="Top"
  tmp_tmp[which.min(tmp_tmp$y),"anchor"]="But"
  tmp_list[[i]] = tmp_tmp
}
tmp = do.call(rbind.data.frame,tmp_list)

plot_ly(tmp, x = ~z, y = ~y, color = ~anchor, marker = list(size = 10))
tmp2 = tmp
slide = 20
tmp2[tmp2$slide==slide,]$y = tmp2[tmp2$slide==slide,]$y-100
plot_ly(tmp2, x = ~z, y = ~y, color = ~anchor, marker = list(size = 10))

tmp3 = tmp2
slide = 20
tmp3[tmp3$slide==slide,]$x = tmp3[tmp3$slide==slide,]$x-300
plot_ly(tmp3, x = ~x, y = ~slide, color = ~anchor, marker = list(size = 10))

plot_ly(tmp3, x = ~z, y = ~y, color = ~anchor, marker = list(size = 10))

for (i in 1:20) {
  xd = unique(tmp3$x[tmp3$slide==i]-tmp$x[tmp$slide==i])[1]
  yd = unique(tmp3$y[tmp3$slide==i]-tmp$y[tmp$slide==i])[1]
  df2[df2$slide==i,]$x = df2[df2$slide==i,]$x+xd
  df2[df2$slide==i,]$y = df2[df2$slide==i,]$y+yd
  
}
df2$x = df2$x+500
df2$y = df2$y-5000
plots = list()
for (i in 1:20) {
  sub = df2[df2$slide==i,]
  plots[[i]] = ggplot(sub, aes(x, y, color = percent.mt, label = X)) +
    geom_point(size = 0.5, shape=20, alpha=1) +
    scale_color_viridis(option = "D", direction = -1)+
    scale_y_continuous(minor_breaks = seq(0, 6000, 100), breaks = seq(0, 6000, 500), limits = c(0,6000))+
    scale_x_continuous(minor_breaks = seq(0, 4000, 100), breaks = seq(0, 4000, 500), limits = c(0,4000))+
    geom_point(data = sub[!is.na(sub$anchor),] ,size = 5, shape=20, alpha=1, color = "magenta")+ coord_fixed()+ theme(legend.position = "none")
}

plot_ly(df2, x = ~x, y = ~y, z = ~z, color = ~percent.mt, marker = list(size = 0.5))

write.table(df, "data/OB2_slideseq/All_slides_3d.txt")

library(patchwork)
wrap_plots(plots, ncol = 5)
ggsave("plots/OB2.jpeg", width = 25, height = 20, units = "in")

