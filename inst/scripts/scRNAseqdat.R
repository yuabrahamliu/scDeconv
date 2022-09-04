# To construct the transcriptomic reference for the human placenta, an 
# scRNA-seq dataset was downloaded from ArrayExpress, with experiment code 
# E-MTAB-6701 (droplet-based data) 
#
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6701/samples/ 
#
# and then transferred to Seurat to create a Seurat object with global-scaling 
# normalization. Cells with less than 500 detected genes and genes detected in 
# less than 5 cells were removed.


