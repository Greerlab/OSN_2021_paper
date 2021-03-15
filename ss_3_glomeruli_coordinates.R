library(openxlsx)
library(ggplot2)
library(ggrepel)
score = read.xlsx("data/OB2_slideseq/OB2_score.xlsx") #The manual scoring result 
tmp = score[score$AVE>=4,] 
write.xlsx(tmp, "data/OB2_slideseq/OB2_score_4up.xlsx")

score = read.xlsx("data/OB2_slideseq/OB2_score_4up.xlsx")
OBM = score[score$MD=="M",c(1,2,5)] #take OR, slide, y column
OBM$index = 1:nrow(OBM)
to_remove = c(28,53,54,55,11) # remove duplicates with lower score
OBM = OBM[!OBM$index%in%to_remove,]

df_list = list()
for (i in 1:20) {
  tmp = read.table(paste0("data/OB2_slideseq/Slide",i,"_3d_coor.txt"), header = T, row.names = 1)
  tmp$slide = paste0(i)
  df_list[[i]] = tmp
}
df = do.call(rbind,df_list)
df$slide = as.numeric(df$slide)
df$z = df$slide*4000 # the estiminate distance between slides
df$OR = "none"

for (i in 1:20) {
  print(paste0("Slide",i))
  dt = readRDS(paste0("Slide",i,".rds"))
  EM = dt@assays$RNA@counts
  EM = EM[grep("Olfr",rownames(EM)),]
  EM[EM>0]=1
  EM = as.matrix(EM)
  EM = EM[rowSums(EM)>1,] # only show OR with 2 more spot per slide
  print(paste0("max OR expressed per spot: ", max(colSums(EM))))
  EM = EM[,colSums(EM)>0]
  for (x in 1:ncol(EM)) {
    tmp = EM[,x]
    OR_name = names(tmp[tmp>0])[1]
    df[colnames(EM)[x],"OR"] = OR_name
  }
  #fact spot for more than 2 different OR
  if (max(colSums(EM))>1) {
    EM2 = EM[,colSums(EM)>1,drop=F]
    EM2 = EM2[rowSums(EM2)>0,,drop=F]
    print(EM2)
    for (y in 1:ncol(EM2)) {
      tmp = EM2[,y]
      OR_name = names(tmp[tmp>0])[2]#checked, 2 is the max for all slides
      barcode = colnames(EM2)[y]
      tmp_data = df[barcode, ,drop=F]
      tmp_data$OR = OR_name
      rownames(tmp_data) = gsub("$","x",barcode)
      df = rbind(df, tmp_data)
    }
  }
}
df$index = 1:nrow(df)
write.table(df, "data/OB2_slideseq/All_coord.txt", sep = "\t", col.names = NA)

OBM$OR = gsub("^","Olfr",OBM$OR)
df$target = "NO"
for (n in 1:nrow(OBM)) {
  df[df$slide==OBM$Slide[n]&df$OR==OBM$OR[n],]$target= "Yes"
}

#remove the noise spots manually before averaging the glomeruli signal position
n = 50
i = OBM$Slide[n]
OR = OBM$OR[n]
df_tmp = df[df$slide==i,]
df_tmp_OR = df_tmp[df_tmp$OR==OR,]
p = ggplot(df_tmp, aes(x, y, color = percent.mt)) +
  geom_point(size = 0.5, shape=20, alpha=1) +
  scale_color_viridis(option = "D", direction = -1) +
  theme_void()
p +geom_point(data = df_tmp_OR, aes(x = x, y =y), col = 'magenta', size =1)+geom_vline(xintercept = -1523)
#+geom_point(data = df_tmp_OR, aes(x = x, y =y-300), col = 'black', size =1)
df_tmp_OR
rows = c(186482,199991,210958,313090,317284,302394,363388,369036,348232,365029,381292,419681,470222,461620,507020,508272,
         504428,509007,511677,517194,542507,544902,557403,550823,588790,549391,593843,615918,642640,648131,651882,648564,
         642484,701319,714189,721552,740474,699027,773477,744363,744437,744810,747229,748567,752349,753645,771322,776440,
         754535,752529,761466,769390,756142)
df[df$index%in%rows,]$OR = "none"

write.table(df, "data/OB2_slideseq/All_coord_OBM_labeled.txt", sep = "\t", col.names = NA)

mean_x = c()
mean_y = c()
mean_z = c()
for (n in 1:nrow(OBM)) {
  i = OBM$Slide[n]
  OR = OBM$OR[n]
  df_tmp = df[df$slide==i,]
  df_tmp_OR = df_tmp[df_tmp$OR==OR,]
  mean_x[n] = mean(df_tmp_OR$x)
  mean_y[n] = mean(df_tmp_OR$y)
  mean_z[n] = mean(df_tmp_OR$z)
}

OBM$ML = mean_x
OBM$DV = mean_y
OBM$AP = mean_z
OBM = OBM[,c("OR","Slide","AVE","ML","DV","AP")]
write.table(OBM, "data/OB2_slideseq/Mean_OBM2_coord.txt", sep = "\t")

OSN = readRDS("data/OR_OSN.rds")
df = OSN@reductions$harmony@cell.embeddings
df$OR_identity = OSN$OR_identity
write.table(df, "data/All_OSN_PC.txt", sep = "\t")

