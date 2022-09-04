
#Interaction model####

covarpcs <- function(covardat = pcdat){

  pcs <- prcomp(covardat, scale. = TRUE)

  cumvar <- cumsum(pcs$sdev^2/sum(pcs$sdev^2))
  pcnum <- which(x = cumvar > 0.8)[1]
  pcnum <- min(pcnum, 5)

  toppcs <- pcs$x[, c(1:pcnum), drop = FALSE]

  #Note, The signs of the columns of the rotation matrix are
  #arbitrary, and so may differ between different programs for
  #PCA, and even between different builds of R.
  #Hence, need to check the correlation between pc1 and the
  #row means of expsub, if it is less than 0, it means the
  #sign of pc1 is reversed, and need to reverse it back
  pc1cor <- cor(rowMeans(apply(covardat, 2, scale)), toppcs[,1])
  if(pc1cor < 0){
    toppcs <- -toppcs
  }

  return(toppcs)

}


#'Find cell type specific differential features from bulk data
#'
#'Find cell type specific inter-group differential features from bulk data and
#'cell content information using a regression model containing an interaction
#'term between the sample group and the specific cell type content.
#'
#'@param cellconts A matrix recording the cell contents of the samples. Each
#'  row is a sample and each column is a cell type. The row names are sample
#'  IDs and the column names are cell type names. The deconvolution result
#'  from \code{refDeconv} and \code{methylpredict} can be directly used here.
#'@param vardat A data frame recording the sample group information, and must
#'  include 2 columns. One is named as "sampleid", recording the sample IDs
#'  same as the row names of \code{cellconts}, the other is "Samplegroup",
#'  recording the sample group to which each sample belongs. If there are any
#'  additional columns in this data frame, they will be deemed as confounding
#'  factors when selecting the cell type specific inter-group differential
#'  features using a linear regression model.
#'@param responsedat A matrix reconding the feature values of the samples, and
#'  each row is a feature and each column is a sample. The column names are
#'  sample IDs and the row names are feature names. The cell type specific
#'  differential features will be selected from the features in this matrix.
#'@param padjcutoff The adjusted p-value cutoff to select the significantly
#'  differential features. Default is NULL, and in this case, the original
#'  p-value will be used to select the differential features instead and its
#'  cutoff is defined by the parameter \code{pcutoff}, but if \code{pcutoff}
#'  is also NULL, the criterion of adjusted p-value < 0.05 will be used to
#'  find the differential features.
#'@param pcutoff If \code{padjcutoff} is set as NULL, this parameter will be
#'  used to set a cutoff on the original p-value to select the significantly
#'  differential features. It default value is also NULL, and in the case that
#'  both \code{padjcutoff} and \code{pcutoff} are NULL, \code{padjcutoff} will
#'  be set as 0.05 automatically and used to find the differential features.
#'@param gradientcutoff The cutoff on the gradient (i.e. the coefficient of
#'  the interaction term in the regression model between a specific cell type
#'  content and the sample group, which reflects the group difference on the
#'  partial derivative of the feature to the cell content) to find the cell
#'  type specific differential features. The default is NULL, and it will be
#'  deemed as 0. The cell type specific differential features between sample
#'  groups are selected using a regression model containing an interaction
#'  term between the sample group and the specific cell type content.
#'@param threads Number of threads need to be used to do the computation. Its
#'  default value is 1.
#'@param int A logical value indicating whether the cell contents need to be
#'  transformed to fit an inverse normal distribution before the differential
#'  feature regression model construction, so that the collinearity caused by
#'  the interaction term in the model can be reduced. Default is TRUE.
#'@return A list with various slots recording the inter-group differential
#'  feature selection results for each cell type. Each slot contains a matrix
#'  and the feature names, p-values, adjusted p-values, etc, are included. In
#'  addition, volcano plots for the differential features for each cell type
#'  will also be generated.
#'@examples 
#'scRNA <- system.file('extdata', 'scRNAseqdat.rds', package = 'scDeconv')
#'scRNA <- readRDS(scRNA)
#'
#'pRNA <- system.file('extdata', 'pairedRNAdat.rds', package = 'scDeconv')
#'pRNA <- readRDS(pRNA)
#'
#'pDNAm <- system.file('extdata', 'pairedDNAmdat.rds', package = 'scDeconv')
#'pDNAm <- readRDS(pDNAm)
#'
#'externalDNAm <- system.file('extdata', 'externalDNAmdat.rds', package = 'scDeconv')
#'externalDNAm <- readRDS(externalDNAm)
#'
#'DNAmpd <- system.file('extdata', 'DNAmpd.rds', package = 'scDeconv')
#'DNAmpd <- readRDS(DNAmpd)
#'
#'pseudobulk <- system.file('extdata', 'pseudobulk.rds', package = 'scDeconv')
#'pseudobulk <- readRDS(pseudobulk)
#'
#'refres <- scRef(Seuratobj = scRNA,  
#'                targetcelltypes = c('EVT', 'FB', 'HB', 'VCT'),  
#'                celltypecolname = 'annotation',  
#'                pseudobulkdat = pseudobulk, 
#'                targetdat = pRNA, 
#'                targetlogged = TRUE)
#'
#'dnamres <- epDeconv(rnaref = refres$ref, 
#'                    rnamat = refres$targetnolog, 
#'                    rnamatlogged = FALSE, 
#'                    methylmat = pDNAm, 
#'                    learnernum = 10, 
#'                    resscale = TRUE, 
#'                    targetmethyldat = externalDNAm)
#'                    
#'totalcellconts <- rbind(dnamres$methylcellconts, 
#'                        dnamres$methyltargetcellcounts)
#'
#'totalvardat <- DNAmpd[c("sampleid", "Samplegroup", "Gestwk")]
#'
#'totalresponse <- cbind(pDNAm, externalDNAm)
#'
#'celldiffprobes <- celldiff(cellconts = totalcellconts,
#'                           vardat = totalvardat, 
#'                           responsedat = totalresponse, 
#'                           pcutoff = 0.01, 
#'                           gradientcutoff = 0)
#'@export
celldiff <- function(cellconts,
                     vardat,
                     responsedat,
                     padjcutoff = NULL,
                     pcutoff = NULL,
                     gradientcutoff = NULL,
                     threads = 1,
                     int = TRUE){

  useM <- FALSE

  converttoM <- function(beta = beta){
    M <- log2(beta/(1 - beta))
    M <- M[complete.cases(M),]
    M <- M[is.finite(rowSums(M)),]
    return(M)
  }

  if(useM == TRUE){
    responsedat <- converttoM(beta = responsedat)
  }

  if(int == TRUE){

    int <- function(x){

      intval <- qnorm((rank(x, na.last = "keep") - 0.5)/sum(!is.na(x)))
      return(intval)

    }

    cellconts <- apply(X = cellconts, MARGIN = 2, int)

  }

  celltypes <- colnames(cellconts)
  varname <- names(vardat)[2]
  covarnames <- names(vardat)[3:ncol(vardat)]

  pddat <- cbind(vardat, cellconts[vardat$sampleid,])
  row.names(pddat) <- pddat$sampleid
  pddat <- pddat[-1]

  responsedat <- responsedat[,row.names(pddat)]
  designdat <- pddat

  for(i in 1:ncol(designdat)){

    if(!is.numeric(designdat[,i])){

      if(!is.factor(designdat[,i])){

        designdat[,i] <- factor(designdat[,i])

      }

      designdat[,i] <- as.numeric(designdat[,i])

      if(length(unique(designdat[,i])) > 1){
        designdat[,i] <- designdat[,i] - 1
      }
    }

  }

  #designdat <- apply(X = designdat, MARGIN = 2, FUN = scale, scale = TRUE)
  row.names(designdat) <- row.names(pddat)
  designdat <- as.data.frame(designdat, stringsAsFactors = FALSE)

  designdatslist <- list()
  i <- 1
  for(i in 1:length(celltypes)){
    celltype <- celltypes[i]

    interterm <- paste0(celltype, ':', varname)
    designdat[interterm] <- designdat[varname]*designdat[celltype]

    pcvars <- setdiff(colnames(designdat), c(interterm, varname, celltype))
    pcdat <- designdat[, pcvars, drop = FALSE]

    if(ncol(pcdat) > 0){
      pcdat <- covarpcs(covardat = pcdat)
      pcnames <- colnames(pcdat)
      designvars <- c(interterm, varname, celltype, pcnames)

      designdatslist[[i]] <- cbind(designdat, pcdat)[designvars]
    }else{
      designvars <- c(interterm, varname, celltype)
      designdatslist[[i]] <- designdat[designvars]
    }

  }

  finddiff <- function(designdats,
                       responsedat,
                       padjcut = NULL,
                       pcut = NULL,
                       gradientcutoff = NULL,
                       annotextsize = 14,
                       textsize = 14,
                       titlesize = 15,
                       face = 'bold'){

    design <- designdats

    if(is.list(design)){
      design <- unlist(design)
    }

    design <- as.numeric(design)
    design <- c(rep(1, nrow(pddat)), design)
    design <- matrix(design, ncol = ncol(designdats) + 1, byrow = FALSE)

    colnames(design) <- c('(Intercept)', colnames(designdats))
    row.names(design) <- row.names(designdats)

    intermidx <- c(1:ncol(design))[grepl(pattern = ':', x = colnames(design),
                                         fixed = TRUE)]
    celltype <- colnames(design)[intermidx]
    celltype <- gsub(pattern = ':.*$', replacement = '', x = celltype)


    fit1 <- limma::lmFit(object = responsedat, design = design)
    fit2 <- limma::eBayes(fit1)

    allg.limma <- limma::topTable(fit2, coef = intermidx, n = dim(fit1)[1])

    if(is.null(gradientcutoff)){
      gradientcutoff <- 0
    }

    if(gradientcutoff < 0){
      gradientcutoff <- 0
    }

    if(!is.null(padjcut)){
      sigg.limma <- subset(allg.limma,
                           adj.P.Val < padjcut & abs(logFC) > gradientcutoff)
      nonsigg.limma <- subset(allg.limma,
                              adj.P.Val >= padjcut | abs(logFC) <= gradientcutoff)
    }else if(!is.null(pcut)){
      sigg.limma <- subset(allg.limma,
                           P.Value < pcut & abs(logFC) > gradientcutoff)
      nonsigg.limma <- subset(allg.limma,
                              P.Value >= pcut | abs(logFC) <= gradientcutoff)
    }else{
      sigg.limma <- subset(allg.limma,
                           adj.P.Val < 0.05 & abs(logFC) > gradientcutoff)
      nonsigg.limma <- subset(allg.limma,
                              adj.P.Val >= 0.05 | abs(logFC) <= gradientcutoff)
    }


    #Draw probe valcano
    if(nrow(sigg.limma) > 0){
      both_pos <- sigg.limma[c('logFC', 'P.Value')]
      both_pos <- as.data.frame(both_pos, stringsAsFactors = FALSE)
      both_pos$Type <- 'NonSig'
    }else{
      both_pos <- data.frame(logFC = numeric(), P.Value = numeric(),
                             Type = character(),
                             stringsAsFactors = FALSE)
    }

    both_pos$Type[both_pos$logFC > 0] <- 'UP'
    both_pos$Type[both_pos$logFC < 0] <- 'DN'

    both_pos_nonsig <- subset(both_pos, Type == 'NonSig')
    both_pos <- subset(both_pos, Type != 'NonSig')

    #library(ggplot2)
    #library(RColorBrewer)

    if(nrow(both_pos) > 50){
      myColor <- densCols(both_pos$logFC, -log10(both_pos$P.Value),
                          colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
    }else{
      myColor <- rep('blue', nrow(both_pos))
    }

    both_pos$denscolor <- myColor
    rgbmat <- t(col2rgb(myColor))
    rgbmat <- as.data.frame(rgbmat, stringsAsFactors = FALSE)
    both_pos <- cbind(both_pos, rgbmat)
    both_pos <- both_pos[order(-both_pos$blue, both_pos$red, both_pos$green),]
    both_pos1 <- subset(both_pos, blue >= red)
    both_pos2 <- subset(both_pos, blue < red)
    both_pos2 <- both_pos2[order(-both_pos2$blue, both_pos2$red, -both_pos2$green),]
    both_pos <- rbind(both_pos1, both_pos2)

    both_nonsig <- nonsigg.limma[c('logFC', 'P.Value')]
    both_nonsig$Type <- 'NonSig'
    both_nonsig <- rbind(both_nonsig, both_pos_nonsig)
    both_nonsig$denscolor <- '#C0C0C0'
    both_nonsig$red <- 192
    both_nonsig$green <- 192
    both_nonsig$blue <- 192
    nonsignum <- nrow(both_nonsig)

    both_pos <- rbind(both_nonsig, both_pos)

    upnum <- nrow(subset(both_pos, Type == 'UP'))
    dnnum <- nrow(subset(both_pos, Type == 'DN'))
    nonsignum <- sum(both_pos$Type == 'NonSig')



    title <- paste0(celltype, ' Differential Features (', celltype, ')')
    if(nrow(sigg.limma) > 0){

      ycut <- -log10(max(sigg.limma$P.Value))

    }else{
      ycut <- NULL
    }


    p <- ggplot2::ggplot(both_pos, ggplot2::aes(x = logFC, y = -log10(P.Value)))

    p <- p + ggplot2::geom_point(color = both_pos$denscolor, position = 'jitter') +
      ggplot2::xlab('Gradient') + ggplot2::ylab('-log10(P.Value)') +
      ggplot2::ggtitle(title,
                       subtitle = paste0('(UP = ', upnum, ', DN = ', dnnum,
                                         ', NonSig = ', nonsignum, ')')) +
      ggplot2::geom_hline(yintercept = ycut, color = 'red') +
      ggplot2::geom_vline(xintercept = gradientcutoff, color = 'blue') +
      ggplot2::geom_vline(xintercept = -gradientcutoff, color = 'blue') +
      ggplot2::geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
      ggplot2::geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = titlesize, face = face),
                     plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                     axis.title = ggplot2::element_text(size = titlesize, face = face),
                     axis.text = ggplot2::element_text(size = textsize, face = face),
                     plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"))

    if(nrow(sigg.limma) > 0){

      sigg.limma$feature <- row.names(sigg.limma)
      row.names(sigg.limma) <- 1:nrow(sigg.limma)
      sigg.limma$celltype <- celltype

    }else{

      sigg.limma <- data.frame(logFC = numeric(),
                               AveExpr = numeric(),
                               t = numeric(),
                               P.Value = numeric(),
                               adj.P.Val = numeric(),
                               B = numeric(),
                               feature = character(),
                               celltype = character(),
                               stringsAsFactors = FALSE)

    }

    allg.limma$feature <- row.names(allg.limma)
    row.names(allg.limma) <- 1:nrow(allg.limma)
    allg.limma$celltype <- celltype

    res <- list(allg.limma = allg.limma,
                sigg.limma = sigg.limma,
                p = p)

    return(res)

  }



  if(threads == 1){
    diffs <- list()
    j <- 1
    for(j in 1:length(designdatslist)){

      diff <- finddiff(designdats = designdatslist[[j]],
                       responsedat = responsedat,
                       padjcut = padjcutoff,
                       pcut = pcutoff,
                       gradientcutoff = gradientcutoff)

      diffs[[j]] <- diff

    }

  }else{

    #library(doParallel)

    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))

    doParallel::registerDoParallel(cl)

    #date()
    `%dopar%` <- foreach::`%dopar%`
    diffs <- foreach::foreach(designdats = designdatslist,
                              .export = NULL) %dopar% {
                                finddiff(designdats,
                                         responsedat = responsedat,
                                         padjcut = padjcutoff,
                                         pcut = pcutoff,
                                         gradientcutoff = gradientcutoff)
                              }

    parallel::stopCluster(cl)

    unregister_dopar()

  }

  reslist <- list()
  i <- 1
  for(i in 1:length(diffs)){

    reslist[[i]] <- diffs[[i]]$allg.limma
    celltype <- celltypes[i]
    celltype <- gsub(pattern = '~ ', replacement = '', x = celltype)
    names(reslist)[i] <- celltype

    reslist[[length(diffs) + i]] <- diffs[[i]]$sigg.limma
    names(reslist)[length(diffs) + i] <- paste0(celltype, '.sig')

    print(diffs[[i]]$p)

  }

  return(reslist)

}


calline <- function(dat = plotdat[plotdat$Samplegroup == 'Control',],
                    xname = 'cellcont',
                    yname = 'featureval'){

  comptab <- data.frame(Set1 = dat[,xname],
                        Set2 = dat[,yname],
                        stringsAsFactors = FALSE)

  lmfit <- lm(formula = Set2~Set1, data = comptab)
  inter <- as.vector(coefficients(lmfit))[1]
  slope <- as.vector(coefficients(lmfit))[2]

  if(inter > 0){
    equalsign <- ' + '
    absinter <- signif(inter, 3)
  }else if(inter < 0){
    equalsign <- ' - '
    absinter <- signif(abs(inter), 3)
  }else{
    equalsign <- ''
    absinter <- ''
  }

  form <- paste0('y = ', signif(slope, 3), 'x', equalsign, absinter)

  x <- min(comptab$Set1) +
    (max(comptab$Set1) - min(comptab$Set1))*0.7
  y <- max(lmfit$fitted.values)*0.9

  anno <- data.frame(x = x,
                     y = y,
                     slope = slope,
                     inter = inter,
                     form = form,
                     stringsAsFactors = FALSE)

  return(anno)

}


anovaplot <- function(featurenames,
                      cellnames,
                      cellconts,
                      vardat,
                      responsedat,
                      plot = TRUE,
                      dotsize = 2,
                      textsize = 13,
                      titlesize = 15,
                      face = 'bold',
                      annotextsize = 6){

  varname <- names(vardat)[2]
  plotdat <- vardat
  row.names(plotdat) <- 1:nrow(plotdat)
  plotdat <- plotdat[,c('sampleid', varname)]

  if(length(featurenames) != length(cellnames)){

    if(length(featurenames) == 1){
      featurenames <- rep(featurenames, length(cellnames))
    }else if(length(cellnames) == 1){
      cellnames <- rep(cellnames, length(featurenames))
    }else{
      cat('The vectors `featurenames` and `cellnames` should have the same length, or
          one of them should have length of 1')
      return(NULL)
    }
  }

  plotdatslist <- list()
  for(i in 1:length(featurenames)){

    featurename <- featurenames[i]
    cellname <- cellnames[i]

    plotdat$cellcont <- cellconts[plotdat$sampleid, cellname]
    plotdat$featureval <- responsedat[featurename, plotdat$sampleid]
    plotdat$cellname <- cellname
    plotdat$featurename <- featurename

    plotdatslist[[i]] <- plotdat

  }

  plotdats <- do.call(rbind, plotdatslist)

  if(!is.factor(plotdats[,varname])){

    plotdats[,varname] <- factor(plotdats[,varname])

  }

  plotdats$cellname <- factor(plotdats$cellname,
                              levels = unique(cellnames), ordered = TRUE)
  plotdats$featurename <- factor(plotdats$featurename,
                                 levels = unique(featurenames), ordered = TRUE)

  subtitle <- c()
  varnamelevels <- levels(plotdats[,varname])
  for(i in 1:length(varnamelevels)){
    varnamelevel <- varnamelevels[i]
    varnamelevelcount <- sum(vardat[,varname] == varnamelevel)
    subtitle <- c(subtitle, paste0(varnamelevel, ' = ', varnamelevelcount))

  }
  subtitle <- paste(subtitle, collapse = '; ')
  subtitle <- paste0('(', subtitle, ')')


  if(plot == TRUE){

    uniqvarname <- unique(plotdats[,varname])
    title <- 'Cell specific feature difference'
    annoslist <- list()

    for(i in 1:length(plotdatslist)){

      plotdat <- plotdatslist[[i]]
      annos <- list()
      for(j in 1:length(uniqvarname)){

        varval <- uniqvarname[j]
        subdat <- plotdat[plotdat[,varname] == varval,]
        anno <- calline(dat = subdat, xname = 'cellcont', yname = 'featureval')
        anno[varname] <- varval
        anno$cellname <- unique(plotdat$cellname)
        anno$featurename <- unique(plotdat$featurename)
        annos[[j]] <- anno

      }

      annos <- do.call(rbind, annos)

      annoslist[[i]] <- annos

    }

    annos <- do.call(rbind, annoslist)
    annos$cellname <- factor(x = annos$cellname,
                             levels = levels(plotdats$cellname),
                             ordered = TRUE)
    annos$featurename <- factor(x = annos$featurename,
                                levels = levels(plotdats$featurename),
                                ordered = TRUE)

    facetcols <- ceiling(sqrt(length(plotdatslist)))


    p <- ggplot2::ggplot(plotdats, ggplot2::aes(x = cellcont,
                                                y = featureval,
                                                color = plotdats[,varname]))

    p <- p + ggplot2::geom_point(position = 'jitter',
                                 size = dotsize) +
      ggplot2::xlab('Cell content') + ggplot2::ylab('Feature value') +

      ggplot2::geom_abline(data = annos,
                           mapping = ggplot2::aes(slope = slope,
                                                  intercept = inter,
                                                  color = annos[,varname]),
                           size = 1) +
      ggplot2::geom_text(data = annos,
                         ggplot2::aes(x= x,
                                      y = y,
                                      label = paste0('bolditalic("', form, '")'),
                                      color = annos[,varname]),
                         parse = TRUE,
                         size = annotextsize,
                         show.legend = FALSE) +

      ggplot2::guides(color = ggplot2::guide_legend(title = varname)) +

      ggplot2::ggtitle(title, subtitle = subtitle) +
      ggplot2::facet_wrap(ggplot2::vars(cellname, featurename), ncol = facetcols,
                          scales = 'free') +
      ggplot2::theme_bw() +

      ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face),
                     plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                     axis.title = ggplot2::element_text(size = titlesize, face = face),
                     axis.text = ggplot2::element_text(size = textsize, face = face),
                     strip.text = ggplot2::element_text(size = textsize, face = face),
                     legend.title = ggplot2::element_text(size = textsize, face = face),
                     legend.text = ggplot2::element_text(size = textsize, face = face))




    print(p)

  }

  return(plotdats)

}



probeannotation <- function(platform = 450,
                            finalprobes){

  if(platform == 27){
    annotation <- 'ilmn12.hg19'
    array <- 'IlluminaHumanMethylation27k'
  }else if(platform == 450){
    annotation <- 'ilmn12.hg19'
    array <- 'IlluminaHumanMethylation450k'
  }else if(platform == 850){
    annotation <- 'ilm10b4.hg19'
    array <- 'IlluminaHumanMethylationEPIC'
  }else{
    cat('The parameter `platform` should be provided a value from 27, 450, and 850\n')
    return(NULL)
  }

  annopackage <- paste0(array, 'anno.', annotation)

  if(!(annopackage %in% installed.packages()[,'Package'])){
    cat(paste0('Package ', annopackage, ' is needed to run this function\n'))
    return(NULL)
  }

  if(!('AnnotationDbi' %in% installed.packages()[,'Package'])){
    cat('Package AnnotationDbi is needed to run this function\n')
    return(NULL)
  }

  if(!('org.Hs.eg.db' %in% installed.packages()[,'Package'])){
    cat(paste0('Package org.Hs.eg.db is needed to run this function\n'))
    return(NULL)
  }


  if(platform == 27){

    probeinfo <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Other
    islandsinfo <- probeinfo[c("CPG_ISLAND_LOCATIONS", "CPG_ISLAND")]
    locinfo <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Locations

    selectedcol <- c("Symbol", "Distance_to_TSS")
    genecol <- 'Symbol'
    featurecol <- 'Distance_to_TSS'
    islandcol <- 'CPG_ISLAND'
  }else if(platform == 450){

    probeinfo <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
    islandsinfo <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC
    locinfo <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations

    selectedcol <- c("UCSC_RefGene_Name", "UCSC_RefGene_Group")
    genecol <- 'UCSC_RefGene_Name'
    featurecol <- 'UCSC_RefGene_Group'
    islandcol <- 'Relation_to_Island'
  }else if(platform == 850){

    probeinfo <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other
    islandsinfo <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC
    locinfo <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations

    selectedcol <- c("UCSC_RefGene_Name", "UCSC_RefGene_Group")
    genecol <- 'UCSC_RefGene_Name'
    featurecol <- 'UCSC_RefGene_Group'
    islandcol <- 'Relation_to_Island'
  }

  finalprobes <- finalprobes[finalprobes %in% row.names(probeinfo)]
  if(length(finalprobes) == 0){
    return(NULL)
  }


  probeinfo <- probeinfo[finalprobes,]
  probeinfodata <- as.data.frame(probeinfo)

  probeinfodata <- probeinfodata[selectedcol]

  islandsinfo <- islandsinfo[finalprobes,]
  islandsinfodata <- as.data.frame(islandsinfo)

  locinfo <- locinfo[finalprobes,]
  locinfo <- as.data.frame(locinfo)

  probeinfodata <- cbind(locinfo, probeinfodata, islandsinfodata)
  probeinfodata$Probe <- row.names(probeinfodata)
  row.names(probeinfodata) <- 1:nrow(probeinfodata)

  if(platform == 27){
    probeinfodata$ENTREZID <- probeinfo$Gene_ID
  }


  orgnizeprobeinfo <- function(colvec){
    parseelement <- function(element){
      elementlist <- unlist(strsplit(x = element, split = ';', fixed = TRUE))

      if(length(elementlist) == 0){
        elementlist <- ''
      }

      return(elementlist)
    }
    collist <- lapply(colvec, parseelement)
    return(collist)
  }

  genenamelist <- orgnizeprobeinfo(colvec = probeinfodata[,genecol])
  featurelist <- orgnizeprobeinfo(colvec = probeinfodata[,featurecol])
  poslist <- orgnizeprobeinfo(colvec = probeinfodata[,islandcol])

  genenamelistlens <- unlist(lapply(X = genenamelist, FUN = length))
  featurelistlens <- unlist(lapply(X = featurelist, FUN = length))
  poslistlens <- unlist(lapply(X = poslist, FUN = length))

  if(sum(genenamelistlens != 1) == 0 & platform == 27){
    probeinfodata$ENTREZID <- gsub(pattern = 'GeneID:', replacement = '',
                                   x = probeinfodata$ENTREZID, fixed = TRUE)
  }

  if(sum(genenamelistlens == 0) == 0){

    datnames <- colnames(probeinfodata)

    chrvec <- rep(probeinfodata[,1], times = genenamelistlens)
    locvec <- rep(probeinfodata[,2], times = genenamelistlens)
    strandvec <- rep(probeinfodata[,3], times = genenamelistlens)

    genenamevec <- unlist(genenamelist)
    featurevec <- unlist(featurelist)

    islandvec <- rep(probeinfodata[,6], times = genenamelistlens)
    posvec <- rep(probeinfodata[,7], times = genenamelistlens)
    probevec <- rep(probeinfodata[,8], times = genenamelistlens)

    if(platform == 27){

      geneidlist <- orgnizeprobeinfo(colvec = probeinfodata$ENTREZID)
      geneidvec <- unlist(geneidlist)
      geneidvec <- gsub(pattern = 'GeneID:', replacement = '',
                        x = geneidvec, fixed = TRUE)

      probeinfodata <- tryCatch({
        data.frame(chrvec, locvec, strandvec,
                   genenamevec,
                   featurevec,
                   islandvec, posvec, probevec,
                   geneidvec,
                   stringsAsFactors = FALSE)
      }, error = function(err){
        probeinfodata
      })

      names(probeinfodata) <- datnames[1:ncol(probeinfodata)]

      probeinfodata <- unique(probeinfodata)

    }else{

      unigeneidvec <- AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db,
                                            keys = unique(genenamevec),
                                            columns = 'ENTREZID',
                                            keytype = 'SYMBOL')
      names(unigeneidvec)[1] <- selectedcol[1]

      probeinfodata <- tryCatch({
        data.frame(chrvec, locvec, strandvec,
                   genenamevec,
                   featurevec,
                   islandvec, posvec, probevec,
                   stringsAsFactors = FALSE)
      }, error = function(err){
        probeinfodata
      })

      names(probeinfodata) <- datnames[1:ncol(probeinfodata)]
      probeinfodata <- unique(probeinfodata)

      if(sum(!(unique(probeinfodata[,selectedcol[1]]) %in%
               unigeneidvec[,selectedcol[1]])) == 0){

        probeinfodata <- merge(probeinfodata, unigeneidvec,
                               by = selectedcol[1],
                               sort = FALSE)

      }

      probeinfodata <- unique(probeinfodata)

    }

    if('Distance_to_TSS' %in% datnames){

      probeinfodata$Distance_to_TSS <- as.integer(probeinfodata$Distance_to_TSS)

    }

    if('CPG_ISLAND' %in% datnames){

      probeinfodata$CPG_ISLAND <- as.logical(probeinfodata$CPG_ISLAND)

    }

  }

  if(platform == 27){

    resnames <- c('Probe', 'chr', 'pos', 'strand',
                  'CPG_ISLAND_LOCATIONS', 'CPG_ISLAND',
                  'Symbol', 'ENTREZID', 'Distance_to_TSS')

  }else{

    resnames <- c('Probe', 'chr', 'pos', 'strand',
                  'Islands_Name', 'Relation_to_Island',
                  'UCSC_RefGene_Name', 'ENTREZID', 'UCSC_RefGene_Group')

  }

  probeinfodata <- probeinfodata[,resnames[unique(c(1:7,
                                                    (ncol(probeinfodata) - 1), 9))]]

  row.names(probeinfodata) <- 1:nrow(probeinfodata)


  return(probeinfodata)

}


geneannotation <- function(genesymbols = NULL,
                           geneentrezs = NULL){

  #summaryannotation <-
  #  readRDS('C:/Users/liuy47/Desktop/Transfer/codestransfer/deconv/files/summaryannotation.rds')

  summaryannotation <- get('summaryannotation')


  if(!is.null(genesymbols)){

    genesymbols <- genesymbols[!is.na(genesymbols)]
    genesymbols <- toupper(genesymbols)

    totalprobegeneannotation <- subset(summaryannotation,
                                       SYMBOL %in% genesymbols |
                                         preferred_name %in% genesymbols)

    part1 <- subset(totalprobegeneannotation,
                    SYMBOL %in% genesymbols)
    part1$input <- part1$SYMBOL

    part2 <- subset(totalprobegeneannotation,
                    !(SYMBOL %in% genesymbols))
    part2$input <- part2$preferred_name

    totalprobegeneannotation <- rbind(part1, part2)
    totalprobegeneannotation <- unique(totalprobegeneannotation)


    others <- genesymbols[!(genesymbols %in%
                              c(totalprobegeneannotation$SYMBOL,
                                totalprobegeneannotation$preferred_name))]


  }else if(!is.null(geneentrezs)){

    geneentrezs <- as.character(geneentrezs)
    geneentrezs <- geneentrezs[!is.na(geneentrezs)]
    totalprobegeneannotation <- subset(summaryannotation,
                                       ENTREZID %in% geneentrezs)
    totalprobegeneannotation$input <- totalprobegeneannotation$ENTREZID
    totalprobegeneannotation <- unique(totalprobegeneannotation)

    others <- geneentrezs[!(geneentrezs %in% totalprobegeneannotation$ENTREZID)]

  }else{

    return(NULL)
  }

  if(length(others) > 0){
    others <- data.frame(ENTREZID = NA,
                         SYMBOL = NA,
                         chr = NA,
                         start = NA,
                         end = NA,
                         strand = NA,
                         preferred_name = NA,
                         UniProt = NA,
                         UniProt.name = NA,
                         protein.size = NA,
                         annotation = NA,
                         input = others,
                         stringsAsFactors = FALSE)
    totalprobegeneannotation <- rbind(totalprobegeneannotation, others)

  }



  totalprobegeneannotation <- unique(totalprobegeneannotation)
  totalprobegeneannotation <- totalprobegeneannotation[
    c('input', 'ENTREZID', 'SYMBOL', 'chr', 'start', 'end', 'strand',
      'preferred_name', 'UniProt', 'UniProt.name', 'protein.size', 'annotation')
  ]

  row.names(totalprobegeneannotation) <- 1:nrow(totalprobegeneannotation)


  return(totalprobegeneannotation)


}


getcorgenes <- function(probes,
                        pairedRNA,
                        pairedmethyl,
                        cutoff = 0.5,
                        generegions = c('TSS200'),
                        platform = 450){

  probes <- intersect(probes, row.names(pairedmethyl))
  pairedprobes <- pairedmethyl[probes, , drop = FALSE]

  pairedprobes <- probetogene(betadat = pairedprobes,
                              platform = platform,
                              group450k850k = generegions,
                              includemultimatch = FALSE)

  pairedprobes <- t(pairedprobes)
  pairedgenes <- t(pairedRNA)

  posgenelist <- c()
  neggenelist <- c()

  genecors <- t(cor(pairedprobes, pairedgenes))

  posgenelist <- c(posgenelist,
                   row.names(which(x = genecors > cutoff, arr.ind = TRUE)))
  neggenelist <- c(neggenelist,
                   row.names(which(x = genecors < -cutoff, arr.ind = TRUE)))
  posgenelist <- unique(posgenelist)
  neggenelist <- unique(neggenelist)

  sharedgeneslist <- intersect(posgenelist, neggenelist)
  posgenelist <- setdiff(posgenelist, sharedgeneslist)
  neggenelist <- setdiff(neggenelist, sharedgeneslist)

  res <- list(posgenes = posgenelist,
              neggenes = neggenelist)

  return(res)


}

getRNAgenes <- function(sigprobes,
                        pairedRNA,
                        pairedmethyl,
                        abscut = 0.5,
                        generegions = c('TSS200'),
                        platform = 450){

  if(nrow(sigprobes) == 0){
    return(list(upgenes = NULL, dngenes = NULL))
  }

  ups <- subset(sigprobes, logFC > 0)
  dns <- subset(sigprobes, logFC <= 0)

  upprobes <- ups$feature
  dnprobes <- dns$feature

  upprobescorgenes <- getcorgenes(probes = upprobes,
                                  pairedRNA = pairedRNA,
                                  pairedmethyl = pairedmethyl,
                                  cutoff = abscut,
                                  generegions = generegions,
                                  platform = platform)
  dnprobescorgenes <- getcorgenes(probes = dnprobes,
                                  pairedRNA = pairedRNA,
                                  pairedmethyl = pairedmethyl,
                                  cutoff = abscut,
                                  generegions = generegions,
                                  platform = platform)

  upprobesposgenes <- upprobescorgenes$posgenes
  upprobesneggenes <- upprobescorgenes$neggenes

  dnprobesposgenes <- dnprobescorgenes$posgenes
  dnprobesneggenes <- dnprobescorgenes$neggenes

  enhancegenes <- unique(c(upprobesposgenes, dnprobesneggenes))
  inhibitgenes <- unique(c(upprobesneggenes, dnprobesposgenes))

  sharedgenes <- intersect(enhancegenes, inhibitgenes)
  enhancegenes <- setdiff(enhancegenes, sharedgenes)
  inhibitgenes <- setdiff(inhibitgenes, sharedgenes)

  res <- list(upprobescorgenes = upprobescorgenes,
              dnprobescorgenes = dnprobescorgenes,
              enhancegenes = enhancegenes,
              inhibitgenes = inhibitgenes)

  return(res)

}


getgenes <- function(sigprobes,
                     platform = 450,
                     generegions = c('TSS200')){

  if(nrow(sigprobes) == 0){
    return(list(upgenes = NULL, dngenes = NULL))
  }

  ups <- subset(sigprobes, logFC > 0)
  dns <- subset(sigprobes, logFC <= 0)

  upgenes <- probeannotation(platform = platform,
                             finalprobes = ups$feature)
  dngenes <- probeannotation(platform = platform,
                             finalprobes = dns$feature)

  if(!is.null(upgenes)){
    upgenes <- subset(upgenes, UCSC_RefGene_Group %in% generegions)

    upgenes <- unique(upgenes$UCSC_RefGene_Name)
  }else{
    upgenes <- NULL
  }

  if(!is.null(dngenes) > 0){
    dngenes <- subset(dngenes, UCSC_RefGene_Group %in% generegions)

    dngenes <- unique(dngenes$UCSC_RefGene_Name)
  }else{
    dngenes <- NULL
  }

  res <- list(upgenes = upgenes, dngenes = dngenes)

  return(res)

}


finduniqgenes <- function(sigprobes,
                          pairedRNA,
                          pairedmethyl,
                          abscut = 0.6,
                          cellnames,
                          platform = 450,
                          generegions = c('TSS200')){

  sigcellnames <- names(sigprobes)[grep(pattern = '.sig', x = names(sigprobes),
                                        fixed = TRUE)]

  enhancegenelist <- list()
  inhibitgenelist <- list()
  bothgenelist <- list()
  i <- 1
  for(i in 1:length(sigcellnames)){

    sigcellname <- sigcellnames[i]

    sigdat <- sigprobes[[sigcellname]]

    if(is.null(pairedRNA) | is.null(pairedmethyl)){

      siggenes <- getgenes(sigdat, platform = platform, generegions = generegions)

      enhancegenes <- siggenes$dngenes
      inhibitgenes <- siggenes$upgenes


    }else{

      siggenes <- getRNAgenes(sigprobes = sigdat,
                              pairedRNA = pairedRNA,
                              pairedmethyl = pairedmethyl,
                              abscut = abscut,
                              generegions = generegions,
                              platform = platform)
      enhancegenes <- siggenes$enhancegenes
      inhibitgenes <- siggenes$inhibitgenes

    }

    sharedgenes <- intersect(enhancegenes, inhibitgenes)
    enhancegenes <- setdiff(enhancegenes, sharedgenes)
    inhibitgenes <- setdiff(inhibitgenes, sharedgenes)

    enhancegenelist[[sigcellname]] <- enhancegenes
    inhibitgenelist[[sigcellname]] <- inhibitgenes
    bothgenelist[[sigcellname]] <- unique(c(enhancegenes, inhibitgenes))

  }

  uniqgenelist <- list()
  i <- 1
  for(i in 1:length(sigcellnames)){
    tmplist <- bothgenelist
    sigcellname <- sigcellnames[i]
    uniqgenelist[[sigcellname]] <- list()

    tmplist[[sigcellname]] <- NULL
    othergenes <- unique(unlist(tmplist))
    uniqgenelist[[sigcellname]][['enhancegenes']] <- setdiff(enhancegenelist[[sigcellname]],
                                                             othergenes)
    uniqgenelist[[sigcellname]][['inhibitgenes']] <- setdiff(inhibitgenelist[[sigcellname]],
                                                             othergenes)

  }

  return(uniqgenelist)

}


enrichrres <- function(genes,
                       write = FALSE,
                       prefix,
                       dbs = c('GO_Biological_Process_2018',
                               'GO_Molecular_Function_2018',
                               'BioPlanet_2019', 'WikiPathways_2019_Human',
                               'KEGG_2019_Human', 'BioCarta_2016',
                               'Reactome_2016', 'NCI-Nature_2016',
                               'Panther_2016')){
  genes <- genes[!is.na(genes)]
  genes <- unique(genes)

  if(length(genes) == 0){
    return(NULL)
  }

  if(sum(genes != '') == 0){
    return(NULL)
  }

  library(enrichR)

  enriched <- enrichr(genes = genes, databases = dbs)

  enrichedsum <- do.call(rbind, enriched)

  if(nrow(enrichedsum) == 0){
    return(NULL)
  }

  enrichedsum$Database <- row.names(enrichedsum)
  enrichedsum$Database <- gsub(pattern = '\\..*$', replacement = '',
                               x = enrichedsum$Database)
  enrichedsum$Term <- gsub(pattern = '_.*$', replacement = '',
                           x = enrichedsum$Term)

  enrichedsum <- enrichedsum[c('Term', 'Database', 'Adjusted.P.value', 'P.value',
                               'Odds.Ratio', 'Combined.Score',
                               'Overlap', 'Genes')]

  if(sum(enrichedsum$Adjusted.P.value < 0.05) >= 10){
    enrichedsum <- subset(enrichedsum, Adjusted.P.value < 0.05)
  }else{
    enrichedsum <- subset(enrichedsum, P.value < 0.01)
    if(nrow(enrichedsum) == 0){
      return(NULL)
    }
  }

  enrichedsum <- enrichedsum[order(enrichedsum$Adjusted.P.value,
                                   enrichedsum$P.value,
                                   -enrichedsum$Combined.Score,
                                   -enrichedsum$Odds.Ratio),]
  row.names(enrichedsum) <- 1:nrow(enrichedsum)
  enrichedsum$Annotation <- prefix


  if(write == TRUE){

    stamp <- Sys.time()
    stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
    stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
    stamp <- paste0('.', stamp)

    write.table(enrichedsum,
                paste0(prefix, '_enrichr', stamp, '.txt'),
                sep = '\t',
                row.names = FALSE,
                quote = FALSE)

  }

  return(enrichedsum)

}


enrich <- function(enhancegenes,
                   inhibitgenes,
                   pairedRNA = NULL,
                   pairedmethyl = NULL,
                   abscut = 0.5,
                   celltype,
                   dbs = c('GO_Biological_Process_2018',
                           'BioPlanet_2019',
                           'Reactome_2016'),
                   write = FALSE){


  #library(eClock)

  if(length(inhibitgenes) > 0){
    upanno <- geneannotation(genesymbols = inhibitgenes)

    if(write == TRUE & nrow(upanno) > 0){

      stamp <- Sys.time()
      stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
      stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
      stamp <- paste0('.', stamp)

      if(is.null(pairedRNA) | is.null(pairedmethyl)){
        write.table(upanno,
                    paste0(celltype, '.up_geneanno', stamp, '.txt'),
                    sep = '\t',
                    row.names = FALSE,
                    quote = FALSE)
      }else{
        write.table(upanno,
                    paste0(celltype, '.inhibit_geneanno', stamp, '.txt'),
                    sep = '\t',
                    row.names = FALSE,
                    quote = FALSE)
      }

    }


  }else{
    upanno <- NULL
  }


  if(length(enhancegenes) > 0){
    dnanno <- geneannotation(genesymbols = enhancegenes)

    if(write == TRUE & nrow(dnanno) > 0){

      stamp <- Sys.time()
      stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
      stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
      stamp <- paste0('.', stamp)

      if(is.null(pairedRNA) | is.null(pairedmethyl)){
        write.table(dnanno,
                    paste0(celltype, '.dn_geneanno', stamp, '.txt'),
                    sep = '\t',
                    row.names = FALSE,
                    quote = FALSE)
      }else{
        write.table(dnanno,
                    paste0(celltype, '.enhance_geneanno', stamp, '.txt'),
                    sep = '\t',
                    row.names = FALSE,
                    quote = FALSE)
      }
    }

  }else{
    dnanno <- NULL
  }



  if(is.null(pairedRNA) | is.null(pairedmethyl)){

    upenrich <- enrichrres(genes = inhibitgenes,
                           write = write,
                           prefix = paste0(celltype, '.up'),
                           dbs = dbs)

    dnenrich <- enrichrres(genes = enhancegenes,
                           write = write,
                           prefix = paste0(celltype, '.dn'),
                           dbs = dbs)


    res <- list(upgeneanno = upanno,
                dngeneanno = dnanno,
                upgeneenrich = upenrich,
                dngeneenrich = dnenrich)

  }else{

    upenrich <- enrichrres(genes = inhibitgenes,
                           write = write,
                           prefix = paste0(celltype, '.inhibit'),
                           dbs = dbs)

    dnenrich <- enrichrres(genes = enhancegenes,
                           write = write,
                           prefix = paste0(celltype, '.enhance'),
                           dbs = dbs)


    res <- list(enhancegeneanno = dnanno,
                inhibitgeneanno = upanno,
                enhancegeneenrich = dnenrich,
                inhibitgeneenrich = upenrich)

  }


  return(res)

}


#'Annotate the function of the cell-type-specific differential CpG sites
#'
#'Annotate the function of the cell-type-specific differential DNA methylation
#'probes called from \code{celldiff}.
#'
#'@param sigprobes The list result returned by \code{celldiff} with the inter-
#'  group differential probe selection results for each cell type.
#'@param celltype The cell type whose differential methylation probes needs to
#'  be annotated. Should be a string of the cell type name.
#'@param pairedRNA The RNA part of the paired RNA-DNA methylation dataset. If
#'  both this parameter and \code{pairedmethyl} are not NULL, a correlation-
#'  based method will be used to find the genes whose RNA expression level in
#'  the RNA data significantly correlated with any of the cell-type-specific
#'  differential probes in the paired methylation data, and then the enriched
#'  function of these selected genes will be analyzed and each of them will
#'  also has its individual function annotated. Should be a matrix with each
#'  column representing a sample and each row for one gene. Row names are gene
#'  symbols and column names are sample IDs. The default value is NULL, and
#'  in this case, the correlation-base method will not be used, and instead,
#'  the function annotation and enrichment analysis will be directly conducted
#'  on the cell-type-specific hyper and hypomethylated DNA methylation probes
#'  in \code{sigprobes}, with the probes mapped to corresponding genes first.
#'@param pairedmethyl The methylation part of the paired RNA-DNA methylation
#'  dataset. If both this parameter and \code{pairedRNA} are not NULL, the
#'  correlation-based method will be used to find the genes with an expression
#'  level significantly correlated to the differential methylation probes and
#'  perform function analysis on them, while if it is NULL, the analysis will
#'  be conducted directly on the genes corresponding to the cell-type-specific
#'  differential methylation probes in \code{sigprobes}. It should be a beta
#'  value matrix with each column representing one sample and each row for a
#'  probe. Row names are probe names. Column names are sample IDs. The sample
#'  IDs should be the same as the ones in \code{rnamat}, because they are data
#'  for paired samples. Default is NULL.
#'@param abscut When the correlation method is used to select the genes with a
#'  significant correlation to the cell-type-specific differential probes, the
#'  ones with an absolute Pearson correlation coefficient value greater than
#'  this parameter value will be selected and used for function analysis. The
#'  default value is 0.6.
#'@param generegions When the function analysis is directly performed on the
#'  genes mapped from the differential methylated probes, if a probe located
#'  in the regions defined by this parameter, its corresponding gene will be
#'  used for the function analysis. Default is "TSS200", so that only probes
#'  within TSS200 regions will be used for gene mapping, and it can also be a
#'  vector, such as \code{c("TSS200", "TSS1500", "1stExon")}, so that probes
#'  within any of these 3 regions will be used for gene mapping. On the other
#'  hand, when the function analysis is performed on the genes significantly
#'  correlated with the changed methylation probes, actually the methylation
#'  probes are merge into gene methylation values first, and then the gene
#'  methylation values will be used to calculated the correlation with the
#'  gene RNA values in the paired RNA data. The gene methylation values are
#'  summarized according to the regions in this parameter \code{generegions}.
#'@param platform A string indicating the platform of the differential probes
#'  recorded in \code{sigprobes}. Default is "450K", can also be "EPIC".
#'@param uniquegenes A logical value and if this parameter is set as TRUE,
#'  after getting the genes for function analysis from the cell type indicated
#'  by \code{celltype}, these genes will be used to compare with the genes
#'  obtained similarly from other cell types recorded in \code{sigprobes}, and
#'  only the ones uniquely contained in the target cell type will be used for
#'  the downstream analysis. If this parameter is FALSE, this filter step will
#'  not be performed. Default is TRUE.
#'@param dbs The databases for gene function enrichment analysis. Default is
#'  \code{c("GO_Biological_Process_2018", "BioPlanet_2019", "Reactome_2016")}.
#'@param write A logical value indicating whether the gene function results
#'  need to be written into txt files in the working directory. The default
#'  value is FALSE.
#'@return A list with four slots recording the gene annotation and function
#'  enrichment results for the enhanced and inhibited genes (when correlation-
#'  based method is used) or for the genes mapped from the hypomethylated and
#'  hypermethylated probes (when correlation-based method is not used), and if
#'  \code{write} is TRUE, txt files will be made to save these results.
#'@examples 
#'scRNA <- system.file('extdata', 'scRNAseqdat.rds', package = 'scDeconv')
#'scRNA <- readRDS(scRNA)
#'
#'pRNA <- system.file('extdata', 'pairedRNAdat.rds', package = 'scDeconv')
#'pRNA <- readRDS(pRNA)
#'
#'pDNAm <- system.file('extdata', 'pairedDNAmdat.rds', package = 'scDeconv')
#'pDNAm <- readRDS(pDNAm)
#'
#'externalDNAm <- system.file('extdata', 'externalDNAmdat.rds', package = 'scDeconv')
#'externalDNAm <- readRDS(externalDNAm)
#'
#'DNAmpd <- system.file('extdata', 'DNAmpd.rds', package = 'scDeconv')
#'DNAmpd <- readRDS(DNAmpd)
#'
#'pseudobulk <- system.file('extdata', 'pseudobulk.rds', package = 'scDeconv')
#'pseudobulk <- readRDS(pseudobulk)
#'
#'refres <- scRef(Seuratobj = scRNA,  
#'                targetcelltypes = c('EVT', 'FB', 'HB', 'VCT'),  
#'                celltypecolname = 'annotation',  
#'                pseudobulkdat = pseudobulk, 
#'                targetdat = pRNA, 
#'                targetlogged = TRUE)
#'
#'dnamres <- epDeconv(rnaref = refres$ref, 
#'                    rnamat = refres$targetnolog, 
#'                    rnamatlogged = FALSE, 
#'                    methylmat = pDNAm, 
#'                    learnernum = 10, 
#'                    resscale = TRUE, 
#'                    targetmethyldat = externalDNAm)
#'                    
#'totalcellconts <- rbind(dnamres$methylcellconts, 
#'                        dnamres$methyltargetcellcounts)
#'
#'totalvardat <- DNAmpd[c("sampleid", "Samplegroup", "Gestwk")]
#'
#'totalresponse <- cbind(pDNAm, externalDNAm)
#'
#'celldiffprobes <- celldiff(cellconts = totalcellconts,
#'                           vardat = totalvardat, 
#'                           responsedat = totalresponse, 
#'                           pcutoff = 0.01, 
#'                           gradientcutoff = 0)
#'
#'hbenrichres <- enrichwrapper(sigprobes =  celldiffprobes, 
#'                             celltype = 'HB', 
#'                             pairedRNA = pRNA, 
#'                             pairedmethyl = pDNAm, 
#'                             dbs = c('Reactome_2016'))
#'@export
enrichwrapper <- function(sigprobes,
                          celltype,
                          pairedRNA = NULL,
                          pairedmethyl = NULL,
                          abscut = 0.6,
                          generegions = c("TSS200"),
                          platform = "450K",
                          uniquegenes = TRUE,
                          dbs = c('GO_Biological_Process_2018',
                                  'BioPlanet_2019',
                                  'Reactome_2016'),
                          write = FALSE){

  if(platform == '450K'){
    platform <- 450
  }else if(platform == 'EPIC'){
    platform <- 850
  }


  sigcellname <- paste0(celltype, '.sig')

  if(uniquegenes == TRUE){

    cellnames <- grepl(pattern = '.sig', x = names(sigprobes), fixed = TRUE)
    cellnames <- names(sigprobes)[!cellnames]

    uniqgenelist <- finduniqgenes(sigprobes = sigprobes,
                                  pairedRNA = pairedRNA,
                                  pairedmethyl = pairedmethyl,
                                  abscut = abscut,
                                  cellnames = cellnames,
                                  platform = platform,
                                  generegions = generegions)

    enhancegenes <- uniqgenelist[[sigcellname]][['enhancegenes']]
    inhibitgenes <- uniqgenelist[[sigcellname]][['inhibitgenes']]

  }else{


    sigdat <- sigprobes[[sigcellname]]

    if(is.null(pairedRNA) | is.null(pairedmethyl)){

      siggenes <- getgenes(sigdat, platform = platform, generegions = generegions)

      enhancegenes <- siggenes$dngenes
      inhibitgenes <- siggenes$upgenes


    }else{

      siggenes <- getRNAgenes(sigprobes = sigdat,
                              pairedRNA = pairedRNA,
                              pairedmethyl = pairedmethyl,
                              abscut = abscut)
      enhancegenes <- siggenes$enhancegenes
      inhibitgenes <- siggenes$inhibitgenes

    }

    sharedgenes <- intersect(enhancegenes, inhibitgenes)
    enhancegenes <- setdiff(enhancegenes, sharedgenes)
    inhibitgenes <- setdiff(inhibitgenes, sharedgenes)

  }

  res <- enrich(enhancegenes = enhancegenes,
                inhibitgenes = inhibitgenes,
                pairedRNA = pairedRNA,
                pairedmethyl = pairedmethyl,
                abscut = abscut,
                celltype = celltype,
                dbs = dbs,
                write = write)

  return(res)

}




summaryfeature <- function(dat, featurecolidx){

  if(!('plyr' %in% installed.packages()[,'Package'])){
    cat('Package plyr is needed to run this function\n')
    return(NULL)
  }

  calmean <- function(block){
    gene <- unique(block$gene)
    subblock <- block[-1]
    submean <- colMeans(subblock)
    submatrix <- data.frame(submean)
    submatrix <- t(submatrix)
    row.names(submatrix) <- gene
    submatrix <- as.data.frame(submatrix)
    return(submatrix)

  }

  features <- dat[,featurecolidx]
  featurefreqs <- table(features)
  unifeatures <- names(featurefreqs[featurefreqs == 1])
  mulfeatures <- names(featurefreqs[featurefreqs > 1])

  unipart <- dat[dat[,featurecolidx] %in% unifeatures,]
  mulpart <- dat[dat[,featurecolidx] %in% mulfeatures,]

  mulpart <- plyr::ddply(.data = mulpart,
                         .variables = c(names(dat)[featurecolidx]),
                         .fun = calmean)

  row.names(mulpart) <- mulpart[,featurecolidx]
  mulpart <- mulpart[-featurecolidx]

  row.names(unipart) <- unipart[,featurecolidx]
  unipart <- unipart[-featurecolidx]

  finaldat <- rbind(unipart, mulpart)
  finaldat <- as.matrix(finaldat)

  return(finaldat)

}



#'Summarize the methylation beta values of probes to genes
#'
#'Summarize the methylation beta values of probes to genes by averaging the
#'  probes located closely to the TSS of a gene.
#'
#'@param betadat A matrix recording the beta values of methylation probes for
#'  samples. Each column represents one sample and each row represents one
#'  probe. The row names are the probe names while the column names should be
#'  sample IDs.
#'@param platform The platform of the probes. Can be set as "27K", "450K", or
#'  "EPIC". Default is "450K".
#'@param range27k A positive number or a vector with two positive numbers. If
#'  the data is from 27K platform, this parameter is needed to define which
#'  probes should be considered as related to a specific gene, and only the
#'  ones with a distance to the TSS of a gene less than the maximum value and
#'  greater than the minimum value of \code{range27k} will be considered as
#'  related to the gene, and the beta values of these probes will be averaged
#'  to get the gene beta value. If it is a single number, the probes with a
#'  distance less than this number and greater than 0 will be attributed to a
#'  gene. Default is 200.
#'@param group450k850k A vector or single string. If the data is based on 450K
#'  or EPIC platform, this parameter is needed to define which probes could be
#'  considered as related to a specific gene. Only the ones located in the
#'  gene regions included in this parameter will be considered as belong to
#'  the gene. The value of this parameter need to be selected from "TSS200",
#'  "TSS1500", "1stExon", "5'UTR", '3'UTR", and "Body". If it is a vector,
#'  such as \code{c("TSS200", "TSS1500", "1stExon")}, the probes within these
#'  3 regions of a gene will be attributed to the gene and their beta values
#'  will be averaged to get the gene beta value. Default value is the string
#'  "TSS200".
#'@param includemultimatch Some probes can be attributed to multiple genes.
#'  If this parameter is TRUE, these probes will be involved into the beta
#'  value calculation for all their related genes. Otherwise, these probes
#'  will be discarded, so that the beta values of all the genes are averaged
#'  only from their uniquely related probes. Default is FALSE.
#'@return A matrix recording the summarized gene beta values for samples.
#'@examples 
#'pDNAm <- system.file('extdata', 'pairedDNAmdat.rds', package = 'scDeconv')
#'pDNAm <- readRDS(pDNAm)
#'
#'pDNAmgenes <- probetogene(betadat = pDNAm, 
#'                          platform = "450K", 
#'                          group450k850k = "TSS200")
#'@export
probetogene <- function(betadat,
                        platform = "450K",
                        range27k = 200,
                        group450k850k = "TSS200",
                        includemultimatch = FALSE){

  if(platform == '450K'){
    platform <- 450
  }else if(platform == 'EPIC'){
    platform <- 850
  }else if(platform == '27K'){
    platform <- 27
  }

  if(min(betadat) < 0 | max(betadat) > 1){

    cat('Need methylation beta value to run this function and values < 0 or > 1 cannot be used\n')
    return(NULL)

  }

  beforeprobes <- row.names(betadat)

  tssinfo <- probeannotation(platform = platform, finalprobes = beforeprobes)

  if(is.null(tssinfo)){
    return(NULL)
  }

  if(platform == 27){
    if(length(range27k) == 1){
      range27k <- c(0, range27k)
    }else{
      range27k <- c(min(range27k), max(range27k))
    }

    tssinfor <- tssinfo[!is.na(tssinfo$Distance_to_TSS),]
    tssinfor <- tssinfor[(tssinfor$Distance_to_TSS >= as.numeric(min(range27k)) &
                            tssinfor$Distance_to_TSS < as.numeric(max(range27k))),]
    tssinfor <- tssinfor[,c('Probe',
                            'Symbol', 'ENTREZID')]
  }else{
    tssinfor <- tssinfo[tssinfo$UCSC_RefGene_Group %in% group450k850k,]
    tssinfor <- tssinfor[,c('Probe',
                            'UCSC_RefGene_Name', 'ENTREZID')]
  }

  tssinfor <- unique(tssinfor)

  probes <- tssinfor$Probe

  if(includemultimatch == FALSE){

    probes <- probes[probes %in% names(table(probes)[table(probes) == 1])]

  }

  tssinfor <- subset(tssinfor, Probe %in% probes)

  tssinfor <- tssinfor[order(tssinfor[,2], tssinfor[,3]),]
  row.names(tssinfor) <- 1:nrow(tssinfor)


  tssinfor$genename <- paste0(tssinfor[,2], '::', tssinfor[,3])
  tssinfor$genename <- gsub(pattern = '^::', replacement = '',
                            x = tssinfor$genename)
  tssinfor$genename <- gsub(pattern = '::$', replacement = '',
                            x = tssinfor$genename)
  tssinfor$genename <- gsub(pattern = '^NA::', replacement = '',
                            x = tssinfor$genename)
  tssinfor$genename <- gsub(pattern = '::NA$', replacement = '',
                            x = tssinfor$genename)


  betadat <- betadat[tssinfor$Probe, , drop = FALSE]
  betadat <- as.data.frame(betadat, stringsAsFactors = FALSE)
  betadat$genename <- tssinfor$genename
  betadat <- betadat[,c(ncol(betadat), 1:(ncol(betadat) - 1))]

  betadat <- summaryfeature(dat = betadat, featurecolidx = 1)

  if(is.null(betadat)){
    return(NULL)
  }

  genenames <- row.names(betadat)
  geneorders <- order(genenames)
  betadat <- betadat[geneorders, , drop = FALSE]

  return(betadat)

}
