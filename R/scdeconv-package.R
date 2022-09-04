#' scDeconv
#'
#'Package to deconvolve bulk DNA methylation data with scRNA-seq data in a 
#'multi-omics manner
#'
#'@docType package
#'@author Yu Liu <yuabrahamliu@gmail.com>
#'@import enrichR
#'@import foreach
#'@importFrom foreach `%dopar%` foreach
#'@import Seurat
#'@importFrom Seurat GetAssayData Idents FindMarkers
#'@import plyr
#'@importFrom plyr ddply
#'@import parallel
#'@importFrom parallel detectCores makeCluster stopCluster
#'@import doParallel
#'@importFrom doParallel registerDoParallel
#'@import AnnotationDbi
#'@importFrom AnnotationDbi select
#'@import org.Hs.eg.db
#'@importFrom org.Hs.eg.db org.Hs.eg.db
#'@import scater
#'@importFrom scater calculateFPKM calculateTPM
#'@import sva
#'@importFrom sva ComBat
#'@import ggplot2
#'@importFrom ggplot2 ggplot aes geom_point geom_line xlab ylab ggtitle theme_bw 
#'geom_boxplot theme annotate element_text geom_abline facet_wrap vars geom_hline 
#'geom_vline element_blank unit guides guide_legend
#'@import scales
#'@importFrom scales hue_pal
#'@import pheatmap
#'@importFrom pheatmap pheatmap
#'@import glmnet
#'@importFrom glmnet cv.glmnet glmnet
#'@import limma
#'@importFrom limma lmFit eBayes topTable
#'@import IlluminaHumanMethylation450kanno.ilmn12.hg19
#'@importFrom IlluminaHumanMethylation450kanno.ilmn12.hg19 Other Islands.UCSC Locations
#'@import IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#'@importFrom IlluminaHumanMethylationEPICanno.ilm10b4.hg19 Other Islands.UCSC Locations
#'@useDynLib scDeconv
#'@name scDeconv
#'@keyword cell deconvolution
#'@keyword scRNA-seq
#'@keyword DNA methylation
#'@keyword co-training
#'@keyword ensemble
#'@keyword cell-type-specific inter-group differential features
NULL
