install.packages("ggplot2","gridExtra","dplyr","RColorBrewer")

 if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") i
BiocManager::install("AnnotationDbi", version = "3.8") 
BiocManager::install("tximport", version = "3.8") 

