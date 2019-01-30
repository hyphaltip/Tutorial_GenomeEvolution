install.packages("ggplot2","gridExtra","dplyr","RColorBrewer","pheatmap")
 if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi", version = "3.8")
BiocManager::install("tximport", version = "3.8")

