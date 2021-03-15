#DESseq2 function
get_sig <- function(res, bm, thresh=0.05){
  resOrdered <- res[order(res$pvalue),]
  resOrdered <- resOrdered[!is.na(resOrdered$padj), ]
  res_sig <- resOrdered[resOrdered$padj <= thresh, ]
  return(res_sig)
}

plot_gene <- function(dds, gene, intgroup, log = F){
  require(ggplot2)
  d <- plotCounts(dds, gene, intgroup, returnData = TRUE, normalized=TRUE, size =10)
  ylabs = "Normalized count"
  if(log==T){
    d[,"count"] = log2(d[,"count"])
    ylabs = "Normalized count (log2)"}
  g <- suppressWarnings(
    ggplot(d, aes_string(x=intgroup, y="count", fill=intgroup))+
      geom_dotplot(binaxis='y', stackdir='center', size =10)+
      ylab(paste0(ylabs))+theme_classic()+ theme(legend.position = "none"))
  return(g)
}

plot_vol <- function(res_sig, cond, lab = rownames(res_sig)) {
  p1 <- EnhancedVolcano(
    res_sig,
    lab,
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 1e-3,
    title = cond,
    subtitle = "",
  )
  return(p1 + labs(subtitle = NULL) + ylab(expression('-Log'['10']*' P adjusted'))+ theme(aspect.ratio=1, legend.position = "none"))
}
