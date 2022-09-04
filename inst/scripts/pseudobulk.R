# To construct the transcriptomic reference for the human placenta, an 
# scRNA-seq dataset was downloaded from ArrayExpress, with experiment code 
# E-MTAB-6701 (droplet-based data) 
#
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6701/samples/ 
#
# and then transferred to Seurat to create a Seurat object with global-scaling 
# normalization. Cells with less than 500 detected genes and genes detected in 
# less than 5 cells were removed.
#
# The function scRef in the package scDeconv can make the transcriptomic 
# reference from this scRNA-seq dataset, and this process contains a step to 
# synthesize pseudo-bulk RNA-seq data from the scRNA-seq data. It can be 
# performed with scRef directly, but using another function prepseudobulk in 
# the package is more suggested to synthesize the pseudo-bulk data first and 
# then transfer the pseudo-bulk data to scRef. In this case, scRef can 
# directly generate the RNA reference from the pseudo-bulk data. The pseudo-
# bulk dataset can also be saved to be used repetitively for different tasks, 
# such as the data here.
