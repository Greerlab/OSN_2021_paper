plot_features <- function(df, targets, options = "viridis", directions = -1) {
  require(ggplot2)
  require(viridis)
  require(jpeg)
  require(ggpubr)
  img = readJPEG("data/projection3.jpg")
  df = df[order(df[,paste0(targets)], decreasing = F),]
  if (grepl("^[[:digit:]]+", targets)) {
    colnames(df)[which(colnames(df)==paste0(targets))] = paste0("gene_",targets)
    targets = paste0("gene_",targets)
  }
  if (grepl("-", targets)) {
    new_name  = gsub("-","_", targets)
    colnames(df)[which(colnames(df)==paste0(targets))] = new_name
    targets = new_name
  }
  if (is.numeric(df[,paste0(targets)])) {
    ggplot(df, aes_string("x", "y", colour = paste0(targets)))+background_image(img)+geom_point()+
      xlim(c(0,2450))+ylim(c(-3130,300))+scale_color_viridis(direction = directions, option = options)+theme_classic()
  } else {
    ggplot(df, aes_string("x", "y", colour = paste0(targets)))+background_image(img)+geom_point()+
      xlim(c(0,2450))+ylim(c(-3130,300))+theme_classic()}
}