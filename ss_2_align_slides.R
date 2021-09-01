#install.packages("spdep")
library(spdep)
library(ggplot2)
library(viridis)
library(plotly)
library(Seurat)
library(patchwork)
i=20 #change i one by one for alignment

dt = readRDS(paste0("data/OB2_slideseq/Slide",i,".rds"))
df = dt@meta.data
df$x = dt@reductions$Spatial@cell.embeddings[,1]
df$y = dt@reductions$Spatial@cell.embeddings[,2]

ggplot(df, aes(x, y, color = percent.mt)) +
  geom_point(size = 0.5, shape=20, alpha=1) +
  scale_color_viridis(option = "D", direction = -1) +
  theme_void()+NoLegend()


df2 = as.data.frame(Rotation(df[,c('x','y')],-110*pi/180))
df2$percent.mt = df$percent.mt


xshift = 1000
yshift = -200
p = ggplot(df2, aes(V1, V2, color = percent.mt)) +
  geom_point(size = 0.5, shape=20, alpha=1) +
  scale_color_viridis(option = "D", direction = -1) +
  theme_void()+ geom_vline(xintercept = xshift, color = "Red")+
  geom_hline(yintercept = yshift, color = "Blue")

toWebGL(p)

df2$V1 = df2$V1-xshift
df2$V2 = df2$V2-yshift
ggplot(df2, aes(V1, V2, color = percent.mt)) +
  geom_point(size = 0.5, shape=20, alpha=1) +
  scale_color_viridis(option = "D", direction = -1)+ geom_vline(xintercept = 0, color = "Red")+
  geom_hline(yintercept = 0, color = "Blue")+
  theme_void()+NoLegend()



df$x = df2$V1
df$y = df2$V2
df = df[,c("x","y","percent.mt")]
write.table(df, file = paste0("data/OB2_slideseq/Slide",i,"_3d_coor.txt"), sep = "\t", col.names = NA)

plots = list()
for (i in 1:20) {
  df = read.table(paste0("data/OB2_slideseq/Slide",i,"_3d_coor.txt"), sep = "\t", header = T, row.names = 1)
  plots[[i]] = ggplot(df, aes(x, y, color = percent.mt)) +
    geom_point(size = 0.5, shape=20, alpha=1) +
    scale_color_viridis(option = "D", direction = -1)+ geom_vline(xintercept = 0, color = "Red")+
    geom_hline(yintercept = 0, color = "Blue")+NoLegend()+xlim(c(-3000,5000))+ylim(c(-6000,2000))
}

jpeg("data/OB2_slideseq/Aligment.jpeg", width = 20, height = 16, res = 150, units = "in")
wrap_plots(plots, ncol = 5)
dev.off()

