pair_analysis_harmony <- function(
  print_plots = F ,
  matrix_input,
  pair_input,
  Ident_input,
  nGenes = 2000,
  Prefix = "pair_analysis",
  use_title = F, ncell = NULL
) 
{ require(ggpubr)
  require(gridExtra)
  require(cluster)
  Idents(matrix_input) = matrix_input@meta.data[,Ident_input]
  Pairs = combn(unique(pair_input),2)
  Silhouette = list()
  ave_Silhouette = c()
  for (i in 1:ncol(Pairs)) {
    target_cell_A = WhichCells(matrix_input, idents = Pairs[1,i])
    target_cell_B = WhichCells(matrix_input, idents = Pairs[2,i])
    if (is.null(ncell)) {
      ncell = min(length(target_cell_A),length(target_cell_B))
    }
    Pair_mtx = subset(matrix_input, cells = c(sample(target_cell_A,ncell),sample(target_cell_B,ncell)))
    Pair_mtx = NormalizeData(Pair_mtx, verbose = FALSE)
    Pair_mtx = FindVariableFeatures(Pair_mtx, verbose = FALSE, nfeatures = nGenes)
    Pair_mtx = ScaleData(Pair_mtx, verbose = FALSE)
    npcs = 50
    if (ncol(Pair_mtx)<50) {
      npcs = ncol(Pair_mtx)-1
    }
    Pair_mtx <- RunPCA(Pair_mtx, verbose = FALSE,npcs = 50)
    Pair_mtx <- RunHarmony(Pair_mtx, group.by.vars = "orig.ident",npcs = 50)
    df = as.data.frame(Pair_mtx@reductions$harmony@cell.embeddings[,1:2])
    ElbowPlot(Pair_mtx, reduction = "harmony")
    df2 = Pair_mtx@meta.data
    df2$cluster = NA
    df2[df2[,Ident_input]==Pairs[1,i],]$cluster=1
    df2[df2[,Ident_input]==Pairs[2,i],]$cluster=2
    Dist = cosine(t(df))
    Dist = 1-Dist
    Sil = silhouette(df2$cluster, Dist)
    round(summary(Sil)[[4]],2)
    ave_Silhouette[i] = round(summary(Sil)[[4]],2)
    
  }
  
  return(ave_Silhouette)
}

