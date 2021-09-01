library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
library(scales)
maindir = "data/OB2_slideseq/Raw_data/"
outdir = "data/OB2_slideseq/"
Slides = 1:20
dirs <- list.files(include.dirs = T,pattern = "puck")
dirs = dirs[order(nchar(dirs), dirs)]
for (i in 1:20) {
  targetdir = paste0(maindir,"/",dirs[i])
  setwd(paste0(maindir,"/",dirs[i]))
  Matrix = grep("AllIllumina",list.files(targetdir, pattern = ".digital_expression.txt.gz"),invert = T, value = T)
  Location = read.csv("BeadLocationsForR.csv", row.names = 1)
  dt = read.table(gzfile(Matrix), header = T, row.names = 1, sep = "\t")
  dt = CreateSeuratObject(dt)
  Location = Location[colnames(dt),]
  colnames(Location) = c("Spatial_1","Spatial_2")
  dt = SCTransform(dt, ncells = 3000, verbose = FALSE) #check
  dt = RunPCA(dt)
  dt <- RunUMAP(dt, dims = 1:10)
  dt <- FindNeighbors(dt, dims = 1:10)
  dt <- FindClusters(dt, resolution = 0.5, verbose = FALSE)
  dt[["Spatial"]] <- CreateDimReducObject(embeddings = as.matrix(Location), key = "Spatial_", assay = DefaultAssay(dt))
  DefaultAssay(dt) = "RNA"
  dt$logUMI = log10(dt$nCount_RNA)
  saveRDS(dt, file = paste0(outdir,"/Slide",i,".rds"))
}

for (i in 1:20) {
  dt = readRDS(paste0(outdir,"/Slide",i,".rds"))
  dt[["percent.mt"]] <- PercentageFeatureSet(dt, pattern = "^mt-")
  saveRDS(dt, file = paste0(outdir,"/Slide",i,".rds"))
  EM = GetAssayData(dt, "counts")
  OR = grep("Olfr", rownames(EM), value = T)
  EM = EM[OR,]
  EM[EM>0]=1
  EM = as.matrix(EM)
  EM = EM[rowSums(EM)>1,]
  n = nrow(EM)
  plot_list = list()
  for (y in 1:n) {
    sub_dt = EM[y,]
    cell_list = names(sub_dt)[sub_dt>0]
    p1 =  FeaturePlot(dt, features = "percent.mt", reduction = "Spatial", order = T, pt.size = 0.1)+ggtitle(paste0(rownames(EM)[y]))+NoLegend()+NoAxes()& 
      scale_color_viridis(option = "D", direction = -1)
    plot_list[[y]] = p1+geom_point(data = as.data.frame(dt@reductions$Spatial@cell.embeddings)[cell_list,], aes(x = Spatial_1, y =Spatial_2), col = 'magenta', size =1)
  }
  w = 32
  h = (ceiling(n/8)*4)+1
  z = 8
  if (n<9) {
    h = 5
    w = 4*n
    z = n
  }
  jpeg(paste0("data/OB2_slideseq/OR_plots/Silde",i,"_OR_distribution_mito.jpeg"), width = w, height = h, res = 300, units = "in")
  print(wrap_plots(plot_list, ncol = z)+ plot_layout(guides = "collect") & theme(legend.position = 'bottom'))
  dev.off()
}

Spot_count =c()
Slide_number = c()
OR_name = c()
position = 0
for (i in 1:20) {
  dt = readRDS(paste0(outdir,"Slide",i,".rds"))
  EM = GetAssayData(dt, "counts")
  OR = grep("Olfr", rownames(EM), value = T)
  EM = EM[OR,]
  EM[EM>0]=1
  EM = as.matrix(EM)
  EM = EM[rowSums(EM)>1,]
  n = nrow(EM)
  for (y in 1:n) {
    position = position+1
    Spot_count[position] = sum(EM[y,]>0)
    Slide_number[position] = i
    OR_name[position] = rownames(EM)[y]
  }
}
df = cbind.data.frame(Spot_count,Slide_number,OR_name)

jpeg(paste0("data/OB2_slideseq/OR_plots/OR_distribution.jpeg"), width = 25, height = 13, res = 150, units = "in")
ggplot(df, aes(x = Slide_number, y = Spot_count))+geom_line()+
  geom_point()+facet_wrap(~ OR_name, scales = "free_y", ncol = 15)
dev.off()