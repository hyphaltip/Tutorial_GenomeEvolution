# Tutorial_GenomeEvolution
Genome Evolution Lecture for Phylogenomics Workshop

Framework for learning some basic Evolutionary Genomics and Comprative steps

Requirements:
 - R
 - R packages “ggplot2”,”gridExtra”,”dplyr”,”RColorBrewer”, "pheatmap"
 - R/Bioconductor packages AnnotationDbi and tximport
 - install in R console with this code
 install.packages("ggplot2","gridExtra","dplyr","RColorBrewer","pheatmap")
 if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi", version = "3.8")
BiocManager::install("tximport", version = "3.8")
 - or run the lines in the script scripts/install_pkg.R

Data
 - there are data in the data folder and you should run  download.sh 
   for an example of how to download data from FungiDB
   data/dna_species.dat lists the names of the species to download. These are the prefixes for the files
   in this folder at FungiDB http://fungidb.org/common/downloads/Current_Release/ - note you want the name
   of folders which list the strain eg Umaydis521 not Umaydis which is where the genome data are located.

 - pre-run OrthoFinder results are in analysis/ortho_set1/Results
    see the Orthogroups.csv and Orthogroups.GeneCount.csv for orthologs contain
    See the Orthogroups_UnassignedGenes.csv for genes which are not in a cluster

scripts 
 - plot_chroms_1.R is an R script to generate some summary graphics from GFF files
 - plot_heatmap_family.R is an R script to generate a heat map for gene family sizes
 - extract_orthologs_single-copy_Afum.py is a python script to 


Run the R script
Rscript scripts/plot_chroms_1.R
 - see the output in plots
Rscript scripts/plot_heatmap_family.R
 - see output in plots - also walk through each line one at a time in Rstudio or elsewhere

run the pythonscript to see how to extract gene families which are single-copy for the two Afumigatus strains
./scripts/extract_orthologs_single-copy_Afum.py 
