
#scDeconv####

compheatmap <- function(compmat, 
                        pcut = 0.05, 
                        disnum = FALSE, 
                        title = NULL, 
                        textsize = 12){
  
  if(ncol(compmat) == 1){
    clustercols <- FALSE
  }else{
    clustercols <- TRUE
  }
  
  if(nrow(compmat) == 1){
    clusterrows <- FALSE
  }else{
    clusterrows <- TRUE
  }
  
  #library(pheatmap)
  
  devforend <- 0.0000001
  
  
  heatmappercentbreaks <- unique(as.vector(compmat))
  heatmappercentbreaks <- heatmappercentbreaks[order(heatmappercentbreaks)]
  heatmappercentbreaks[1] <- heatmappercentbreaks[1] - devforend
  heatmappercentbreaks[length(heatmappercentbreaks)] <- heatmappercentbreaks[length(heatmappercentbreaks)] + devforend
  
  tmp <- heatmappercentbreaks
  
  heatmappercentbreaks <- tmp
  
  if(sum(heatmappercentbreaks > 0) > 0){
    heatmappercentbreaks1 <- as.vector(quantile(heatmappercentbreaks[heatmappercentbreaks > 0], 
                                                probs = seq(0, 1, by = (1 - 0)/49)))
  }else{
    heatmappercentbreaks1 <- 0 + devforend
  }
  
  if(sum(heatmappercentbreaks < 0) > 0){
    heatmappercentbreaks2 <- as.vector(quantile(heatmappercentbreaks[heatmappercentbreaks < 0], 
                                                probs = seq(0, 1, by = (1 - 0)/49)))
  }else{
    heatmappercentbreaks2 <- 0 - devforend
  }
  
  
  mediatepoint1 <- sum(min(heatmappercentbreaks1), 0)/2
  mediatepoint2 <- sum(max(heatmappercentbreaks2), 0)/2
  
  if(mediatepoint1 < abs(mediatepoint2)){
    mediatepoint2 <- -mediatepoint1
  }else{
    mediatepoint1 <- -mediatepoint2
  }
  
  heatmappercentbreaks1 <- c(mediatepoint1, heatmappercentbreaks1)
  heatmappercentbreaks2 <- c(heatmappercentbreaks2, mediatepoint2)
  heatmappercentbreaks <- c(heatmappercentbreaks2, 0, heatmappercentbreaks1)
  
  heatmappercentbreaks1 <- unique(heatmappercentbreaks1)
  heatmappercentbreaks2 <- unique(heatmappercentbreaks2)
  heatmappercentbreaks <- unique(heatmappercentbreaks)
  
  if(sum(is.na(heatmappercentbreaks2)) > 0){
    colors <- c(colorRampPalette(c('blue', 'cyan'))(length(unique(heatmappercentbreaks2))-1), 'lightgray')
  }else if(sum(is.na(heatmappercentbreaks1)) > 0){
    colors <- c('lightgray', colorRampPalette(c('orange', 'red'))(length(unique(heatmappercentbreaks1))-1))
  }else{
    colors <- c(colorRampPalette(c('blue', 'cyan'))(length(unique(heatmappercentbreaks2))-1), rep('lightgray', 2), 
                colorRampPalette(c('orange', 'red'))(length(unique(heatmappercentbreaks1))-1))
    
  }
  
  heatmappercentbreaks <- heatmappercentbreaks[!is.na(heatmappercentbreaks)]
  
  if(sd(compmat) != 0){
    if(is.null(title)){
      title <- paste0('Comparison SS (significant level = ', pcut, ')')
    }else{
      title <- paste(toupper(substr(title, 1, 1)), 
                     substr(title, 2, nchar(title)), sep = '')
      title <- paste0(title, ' Comparison SS (significant level = ', pcut, ')')
    }
    
    #pdf(paste0('RPC (', tag, ').pdf'), width = 9, height = 7)
    print(
      pheatmap::pheatmap(compmat, color = colors, 
                         breaks = heatmappercentbreaks, 
                         main = title, 
                         display_numbers = disnum, 
                         cluster_rows = clusterrows, 
                         cluster_cols = clustercols, 
                         angle_col = 45, 
                         fontsize = textsize)
    )
    
    #dev.off()
  }
  
  
}


#'Draw box plot and heatmaps for cell deconvolution result
#'
#'Draw box plot and heatmaps to compare cell deconvolution results between/ 
#'among different sample groups, or a box plot for only one sample group
#'
#'@param cellcontres A data frame or matrix recording the cell contents for 
#'  the samples. Each column is one cell type and each row is one sample. The 
#'  column names are the cell type names and the row names are the sample IDs. 
#'@param pddat A data frame recording the sample group information, and must 
#'  include 2 columns. One is named as "sampleid", recording the sample IDs 
#'  same as the row names of \code{cellcontres}, the other is "Samplegroup", 
#'  recording the sample group to which each sample belongs. It can also be 
#'  set as NULL, meaning all the samples are from the same group, and in this 
#'  case, only a box plot showing the cell content for each cell type will be 
#'  made, but if a data frame is provided to this parameter showing samples 
#'  belong to different groups, for each cell type in the box plot, it will be 
#'  further divided into different groups to compare their cell compositions, 
#'  and also heatmaps will be generated to show the comparison results.
#'@param title A string to define the prefixes of the plot titles, can also be 
#'  set as NULL.
#'@param titlesize A number to define the font size of the plot titles, and 
#'  the default value is 15.
#'@param textsize A number to define the base font size for the plots, except 
#'  the title size and the annotation label size of the box plot, which are 
#'  defined by the parameters \code{titlesize} and \code{annotextsize}. The 
#'  default value is 13.
#'@param annotextsize A number used to define the font size of the annotation 
#'  labels in the box plot. Default is 6.
#'@param face A string to define whether the text in the box plot should be 
#'  bold or plain. Default is "bold".
#'@param annotextheight A number to define the y-coordinate of the annotation 
#'  text in the box plot. Default is 0.95, meaning the y-coordinate of these 
#'  labels will be 95% of the box plot height.
#'@return Generate a box plot to show the sample cell contents, and if these 
#'  samples belong to different groups, also heatmaps will be made to show 
#'  the cell content comparison results among/between different groups.
#'@export
celldeconvplots <- function(cellcontres, 
                            pddat = NULL, 
                            title = NULL, 
                            titlesize = 15, 
                            textsize = 13, 
                            annotextsize = 6, 
                            face = 'bold', 
                            annotextheight = 0.95){
  
  if(is.null(pddat)){
    pddat <- data.frame(sampleid = row.names(cellcontres), 
                        samplegroup = 'Samples', 
                        stringsAsFactors = FALSE)
    
  }
  
  sharedsamples <- intersect(pddat$sampleid, row.names(cellcontres))
  pddat <- subset(pddat, sampleid %in% sharedsamples)
  cellcontres <- cellcontres[pddat$sampleid,]
  
  if(!('Samplegroup' %in% names(pddat))){
    pddat$Samplegroup <- 'Samples'
  }
  
  if(!is.factor(pddat$Samplegroup)){
    
    unisamplegroup <- unique(pddat$Samplegroup)
    unisamplegroup <- unisamplegroup[order(unisamplegroup)]
    
    pddat$Samplegroup <- factor(pddat$Samplegroup, 
                                levels = unisamplegroup, 
                                ordered = TRUE)
    
  }
  
  samplegroups <- levels(pddat$Samplegroup)
  pddat$Samplegroup <- as.character(pddat$Samplegroup)
  
  groupsamples <- list()
  groupdeconv <- list()
  groupmeans <- list()
  
  i <- 1
  for(i in 1:length(samplegroups)){
    
    samplegroup <- samplegroups[i]
    subpd <- subset(pddat, Samplegroup == samplegroup)
    
    if(nrow(subpd) == 0){
      next()
    }
    
    idx <- length(groupsamples) + 1
    
    groupsamples[[idx]] <- subpd$sampleid
    names(groupsamples)[idx] <- samplegroup
    
    groupdeconv[[idx]] <- cellcontres[groupsamples[[idx]],]
    names(groupdeconv)[idx] <- samplegroup
    
    groupmeans[[idx]] <- colMeans(groupdeconv[[idx]])
    groupmeans[[idx]] <- groupmeans[[idx]][order(-groupmeans[[idx]])]
    names(groupmeans)[idx] <- samplegroup
    
    
  }
  
  samplegroups <- names(groupsamples)
  
  celltypes <- colnames(cellcontres)
  
  
  if(length(samplegroups) >= 2){
    
    i <- 1
    for(i in 1:length(celltypes)){
      celltype <- celltypes[i]
      
      j <- 1
      for(j in 1:(length(samplegroups)-1)){
        samplegroup1 <- samplegroups[j]
        samplegroup2 <- samplegroups[j + 1]
        grouppair <- paste0(samplegroup2, '_', samplegroup1)
        
        group1celltypeconts <- groupdeconv[[samplegroup1]][,celltype]
        group2celltypeconts <- groupdeconv[[samplegroup2]][,celltype]
        
        pval <- kruskal.test(list(group2celltypeconts, 
                                  group1celltypeconts))$p.val
        
        lss <- c(1, -1)[c(mean(group2celltypeconts) > mean(group1celltypeconts), 
                          mean(group2celltypeconts) <= mean(group1celltypeconts))]
        lss <- -lss*log10(pval)
        
        if(j == 1){
          line <- pval
          liness <- lss
          grouppairs <- grouppair
        }else{
          line <- c(line, pval)
          liness <- c(liness, lss)
          grouppairs <- c(grouppairs, grouppair)
        }
        
      }
      
      line <- matrix(line, nrow = 1)
      liness <- matrix(liness, nrow = 1)
      
      row.names(line) <- row.names(liness) <- celltype
      
      colnames(line) <- colnames(liness) <- grouppairs
      
      if(i == 1){
        compmatrix <- line
        compmatrixss <- liness
      }else{
        compmatrix <- rbind(compmatrix, line)
        compmatrixss <- rbind(compmatrixss, liness)
      }
    }
    
    
    compmatrix <- compmatrix[complete.cases(compmatrix), , drop = FALSE]
    compmatrixss <- compmatrixss[complete.cases(compmatrixss), , drop = FALSE]
    
    compmatrixss2 <- compmatrixss
    compmatrixss[abs(compmatrixss) < -log10(0.05)] <- 0
    compmatrixss2[abs(compmatrixss2) < -log10(0.01)] <- 0
    
    
    compheatmap(compmat = compmatrixss, pcut = 0.05, disnum = TRUE, title = title, 
                textsize = textsize)
    compheatmap(compmat = compmatrixss2, pcut = 0.01, disnum = TRUE, title = title, 
                textsize = textsize)
    
  }
  
  for(i in 1:ncol(cellcontres)){
    cell <- celltypes[i]
    
    for(j in 1:length(groupdeconv)){
      cluster <- names(groupdeconv)[j]
      sub <- data.frame(content = groupdeconv[[cluster]][,cell], 
                        Group = cluster, 
                        cell = cell, 
                        stringsAsFactors = FALSE)
      if(j == 1){
        subs <- sub
        
      }else{
        subs <- rbind(subs, sub)
      }
      
    }
    
    if(i == 1){
      subses <- subs
    }else{
      subses <- rbind(subses, subs)
    }
  }
  
  subses$cell <- factor(subses$cell, levels = celltypes, 
                        ordered = TRUE)
  #subses$Group <- factor(subses$Group, levels = unique(subses$Group), ordered = TRUE)
  
  
  #library(ggplot2)
  if(is.null(title)){
    title <- 'Deconvolved Cell Contents'
  }else{
    title <- paste(toupper(substr(title, 1, 1)), 
                   substr(title, 2, nchar(title)), sep = '')
    title <- paste0(title, ' Deconvolved Cell Contents')
  }
  
  if(length(samplegroups) == 1){
    subtitle <- paste0('(', nrow(pddat), ' Samples)')
  }else{
    
    groupsamplenum <- table(pddat$Samplegroup)
    subtitle <- paste0(names(groupsamplenum), ' = ', as.vector(groupsamplenum))
    subtitle <- paste(subtitle, collapse = ', ')
    subtitle <- paste0('(', subtitle, ')')
    
  }
  
  p <- ggplot2::ggplot(subses, ggplot2::aes(x = cell, y = content, fill = Group))
  
  p <- p + ggplot2::geom_boxplot() + 
    ggplot2::xlab('Cell Type') + ggplot2::ylab('Cell Content') + 
    ggplot2::ggtitle(title, subtitle = subtitle) + 
    ggplot2::theme_bw()
  
  if(length(samplegroups) == 1){
    
    p <- p + ggplot2::theme(legend.position = 'none')
    
  }else{
    
    labelheight <- as.vector(quantile(unlist(cellcontres), annotextheight))
    
    i <- 1
    for(i in 1:length(celltypes)){
      
      celltype <- celltypes[i]
      cellcontres <- cellcontres[pddat$sampleid,]
      celltypecont <- data.frame(cellcont = cellcontres[,celltype], 
                                 Samplegroup = pddat$Samplegroup, 
                                 stringsAsFactors = FALSE)
      #pval <- wilcox.test(formula = cellcont~Samplegroup, data = celltypecont)$p.val
      pval <- kruskal.test(formula = cellcont~Samplegroup, data = celltypecont)$p.val
      pval <- signif(pval, 3)
      
      p <- p + ggplot2::annotate('text', 
                                 label = pval, 
                                 x = i, 
                                 y = labelheight, 
                                 size = annotextsize, 
                                 color = c('red', 'blue')[c(pval < 0.05, 
                                                            pval >= 0.05)], 
                                 parse = FALSE, 
                                 fontface = 4)
      
    }
    
  }
  
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face), 
                          plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                          axis.title = ggplot2::element_text(size = textsize, face = face), 
                          axis.text = ggplot2::element_text(size = textsize, face = face), 
                          legend.title = ggplot2::element_text(size = textsize, face = face), 
                          legend.text = ggplot2::element_text(size = textsize, face = face))
  
  
  print(p)
  
}


lm_eqn <- function(pcc){
  eq <- substitute(~~bolditalic(R^2~"="~r2), 
                   list(r2 = format(pcc^2, digits = 3)))
  as.character(as.expression(eq));
}

lm_eqn2 <- function(pcc){
  eq <- substitute(~~bolditalic(PCC~"="~r2), 
                   list(r2 = format(pcc, digits = 3)))
  as.character(as.expression(eq));
}


getplotdat <- function(celltypes, 
                       dat1, 
                       dat2, 
                       colorful = TRUE, 
                       usePCC = FALSE){
  
  i <- 1
  for(i in 1:length(celltypes)){
    celltype <- celltypes[i]
    comptab <- data.frame(Set1 = dat1[,celltype], 
                          Set2 = dat2[,celltype], 
                          celltype = celltype, 
                          samples = row.names(dat1), 
                          stringsAsFactors = FALSE)
    
    pcc <- cor(comptab$Set1, comptab$Set2)
    pcc <- signif(pcc, 3)
    
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
    
    xrange <- range(comptab$Set1)
    yrange <- range(comptab$Set2)
    xlen <- abs(xrange[2] - xrange[1])
    ylen <- abs(yrange[2] - yrange[1])
    x1 <- xrange[1] + xlen*1/4
    y1 <- yrange[1] + ylen*4/5
    x2 <- xrange[1] + xlen*2/3
    y2 <- yrange[1] + ylen*1/4
    
    myColor <- rep('blue', nrow(comptab))
    
    if(colorful == TRUE & nrow(comptab) > 50){
      
      myColor <- densCols(x = comptab$Set1, y = comptab$Set2, 
                          colramp = colorRampPalette(rev(rainbow(10, end=4/6))))
      
    }
    
    comptab$denscolor <- myColor
    rgbmat <- t(col2rgb(myColor))
    rgbmat <- as.data.frame(rgbmat, stringsAsFactors = FALSE)
    comptab <- cbind(comptab, rgbmat)
    comptab <- comptab[order(-comptab$blue, comptab$red, comptab$green),]
    comptab1 <- subset(comptab, blue >= red)
    comptab2 <- subset(comptab, blue < red)
    comptab2 <- comptab2[order(-comptab2$blue, comptab2$red, -comptab2$green),]
    comptab <- rbind(comptab1, comptab2)
    
    if(i == 1){
      comptabs <- comptab
      x1s <- x1
      x2s <- x2
      y1s <- y1
      y2s <- y2
      slopes <- slope
      inters <- inter
      forms <- form
      pccs <- lm_eqn(pcc)
      if(usePCC == TRUE){
        pccs <- lm_eqn2(pcc)
      }
      
      
    }else{
      comptabs <- rbind(comptabs, comptab)
      x1s <- c(x1s, x1)
      x2s <- c(x2s, x2)
      y1s <- c(y1s, y1)
      y2s <- c(y2s, y2)
      slopes <- c(slopes, slope)
      inters <- c(inters, inter)
      forms <- c(forms, form)
      if(usePCC == FALSE){
        pccs <- c(pccs, lm_eqn(pcc))
      }else{
        pccs <- c(pccs, lm_eqn2(pcc))
      }
      
    }
    
  }
  
  annos <- data.frame(x1s = x1s, 
                      x2s = x2s, 
                      y1s = y1s, 
                      y2s = y2s, 
                      slopes = slopes, 
                      inters = inters, 
                      forms = forms, 
                      pccs = pccs, 
                      celltype = celltypes, 
                      stringsAsFactors = FALSE)
  
  plotdat <- list(comptabs = comptabs, 
                  annos = annos)
  return(plotdat)
  
}

#'Draw scatter plots to compare two sets of cell content values
#'
#'Draw scatter plots to compare two sets of cell content values and compare 
#'two sets of cell-cell correlation values. For each cell type, its scatter 
#'plot accounts for one facet of all the plots.
#'
#'@param dat1 A data frame or matrix recording the sample cell contents in the 
#'  first dataset. Each column is a cell type and each row is a sample. The 
#'  column names are the cell type names and the row names are the sample IDs. 
#'@param dat2 A data frame or matrix recording the cell contents for the same 
#'  samples in the second dataset. Each column is a cell type and each row is 
#'  a sample. The column names are the cell type names and the row names are 
#'  the sample IDs. 
#'@param title A string to define the prefixes of the plot titles, can also be 
#'  set as NULL.
#'@param colorful A logical value indicating whether the dots in the scatter 
#'  plots should be rainbow colored to indicate their density or not. Default 
#'  is TRUE. This parameter only works when the dot number in the plot is more 
#'  than 50, otherwise, all the dots will be colored as blue. 
#'@param corcmp A logical value indicating whether a scatter plot comparing 
#'  the cell-cell correlation in the 2 datasets should be drawn or not. The 
#'  default value is TRUE.
#'@param usePCC If this value is TRUE, the Pearson correlation coefficient 
#'  between the cell contents will be labeled in the plot. If FALSE, the R 
#'  square value will be labeled instead. Default is FALSE.
#'@param xtitle The title of the x-axes of the scatter plots. Default is "RNA 
#'  cell contents".
#'@param ytitle Title of the y-axes. Default is "Methylation cell contents". 
#'@param annotextsize A number used to define the font size of the annotation 
#'  labels in the scatter plot. Default is 6.
#'@param textsize A number to define the base font size for the plots, except 
#'  the plot title size, axis title size, and the annotation label size, which 
#'  are defined by code{titlesize} and \code{annotextsize}. The default value 
#'  is 13.
#'@param titlesize A number to define the size for the plot and axis titles, 
#'  and the default value is 15.
#'@param face A string to define whether the text in the plot should be bold 
#'  or plain. Default is "bold".
#'@param dotsize A number to define the dot size in the scatter plots. Default 
#'  is 1.
#'@return Generate a scatter plot to compare the cell contents for the same 
#'  samples recorded by two datasets. For each cell type, a facet will be made 
#'  for its plot. If the parameter \code{corcmp} is set as TRUE, an additional 
#'  scatter plot comparing the cell-cell correlation between the two datasets 
#'  will also be drawn.
#'@export
cellscatterplots <- function(dat1, 
                             dat2, 
                             title = NULL, 
                             colorful = TRUE, 
                             corcmp = TRUE, 
                             usePCC = FALSE, 
                             xtitle = 'RNA cell contents', 
                             ytitle = 'Methylation cell contents', 
                             annotextsize = 6, 
                             textsize = 13, 
                             titlesize = 15, 
                             face = 'bold', 
                             dotsize = 1){
  
  dat2 <- dat2[,colnames(dat1), drop = FALSE]
  dat2 <- dat2[row.names(dat1), , drop = FALSE]
  
  if(is.null(title)){
    title1 <- 'Cell Contents Comparison'
  }else{
    title1 <- paste(toupper(substr(title, 1, 1)), 
                    substr(title, 2, nchar(title)), sep = '')
    title1 <- paste0(title1, ' Cell Contents Comparison')
  }
  
  subtitle1 <- paste0('(', ncol(dat1), ' cell types on ', 
                      nrow(dat1), ' samples)')
  
  celltypes <- colnames(dat1)
  
  
  plotdats <- getplotdat(celltypes = celltypes, 
                         dat1 = dat1, 
                         dat2 = dat2, 
                         colorful = colorful, 
                         usePCC = usePCC)
  
  comptabs <- plotdats$comptabs
  annos <- plotdats$annos
  
  
  p <- ggplot2::ggplot(comptabs, ggplot2::aes(x = Set1, y = Set2))
  
  print(
    p + ggplot2::geom_point(color = comptabs$denscolor, position = 'jitter', 
                            size = dotsize) + 
      ggplot2::xlab(xtitle) + 
      ggplot2::ylab(ytitle) + 
      ggplot2::ggtitle(title1, subtitle = subtitle1) + 
      ggplot2::geom_abline(data = annos, 
                           mapping = ggplot2::aes(slope = slopes, 
                                                  intercept = inters), 
                           color = 'red', size = 1) + 
      #ggplot2::geom_abline(slope = slope, intercept = inter, color = 'red', size = 1) + 
      
      ggplot2::geom_text(data = annos, 
                         ggplot2::aes(x= x1s, 
                                      y = y1s, 
                                      label = pccs), 
                         parse = TRUE, 
                         color = 'red', 
                         size = annotextsize) + 
      ggplot2::geom_text(data = annos, 
                         ggplot2::aes(x = x2s, 
                                      y = y2s, 
                                      label = paste0('bolditalic("', forms, '")')), 
                         color = 'red', parse = TRUE, 
                         size = annotextsize) + 
      ggplot2::facet_wrap(ggplot2::vars(celltype), scales = 'free') + 
      ggplot2::theme_bw() + 
      
      ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face), 
                     plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                     axis.title = ggplot2::element_text(size = titlesize, face = face), 
                     axis.text = ggplot2::element_text(size = textsize, face = face), 
                     strip.text = ggplot2::element_text(size = textsize, face = face))
  )
  
  
  if(corcmp == TRUE){
    
    cellcormat1 <- cor(dat1)
    cellcormat2 <- cor(dat2)
    cellcormat <- cor(cellcormat1, cellcormat2)
    #The row names are the col names of the first matrix
    #The col names are the col names of the second matrix
    
    corplotdats <- getplotdat(celltypes = celltypes, 
                              dat1 = cellcormat1, 
                              dat2 = cellcormat2, 
                              colorful = colorful)
    corcomptabs <- corplotdats$comptabs
    corannos <- corplotdats$annos
    corcomptabs$cellpair <- paste0(corcomptabs$celltype, '-', 
                                   corcomptabs$samples)
    
    if(is.null(title)){
      title2 <- 'Cell Correlation Comparison'
    }else{
      title2 <- paste(toupper(substr(title, 1, 1)), 
                      substr(title, 2, nchar(title)), sep = '')
      title2 <- paste0(title2, ' Cell Correlation Comparison')
    }
    
    p <- ggplot2::ggplot(corcomptabs, ggplot2::aes(x = Set1, y = Set2))
    
    print(
      p + 
        ggplot2::geom_point(color = 'blue', position = 'jitter') + 
        ggplot2::xlab(paste0('PCCs with ', xtitle)) + 
        ggplot2::ylab(paste0('PCCs with ', ytitle)) + 
        ggplot2::ggtitle(title2) + 
        ggplot2::geom_abline(data = corannos, 
                             mapping = ggplot2::aes(slope = slopes, 
                                                    intercept = inters), 
                             color = 'red', size = 1) + 
        #ggplot2::geom_abline(slope = slope, intercept = inter, color = 'red', size = 1) + 
        
        ggplot2::geom_text(data = corannos, 
                           ggplot2::aes(x= x1s, 
                                        y = y1s, 
                                        label = pccs), 
                           parse = TRUE, 
                           color = 'red', 
                           size = annotextsize) + 
        ggplot2::geom_text(data = corannos, 
                           ggplot2::aes(x = x2s, 
                                        y = y2s, 
                                        label = paste0('bolditalic("', forms, '")')), 
                           color = 'red', parse = TRUE, 
                           size = annotextsize) + 
        ggplot2::geom_text(data = corcomptabs, 
                           ggplot2::aes(x = Set1, 
                                        y = Set2, 
                                        label = paste0('bolditalic("', samples, '")')), 
                           color = 'blue', parse = TRUE, 
                           check_overlap = TRUE, 
                           size = annotextsize) + 
        
        ggplot2::facet_wrap(ggplot2::vars(celltype), scales = 'free') + 
        ggplot2::theme_bw() + 
        
        ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face), 
                       plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                       axis.title = ggplot2::element_text(size = titlesize, face = face), 
                       axis.text = ggplot2::element_text(size = textsize, face = face), 
                       strip.text = ggplot2::element_text(size = textsize, face = face))
    )
    
  }
  
}


rowSD <- function(mat){
  
  for(i in 1:nrow(mat)){
    line <- mat[i,]
    linesd <- sd(line)
    if(i == 1){
      linesds <- linesd
    }else{
      linesds <- c(linesds, linesd)
    }
  }
  
  return(linesds)
}

#targetdat can be log transformed or not, and this status should be indicated 
#with the parameter targetlogged
#refdat must be non-log transformed
unifydats <- function(targetdat = simuexp, 
                      targetlogged = FALSE, 
                      refdat = scref, 
                      adjustminus = TRUE){
  
  pseudocount <- 1
  
  if(targetlogged == TRUE){
    
    if(sum(2^targetdat - 1 < 0) > 0){
      
      pseudocount <- 0
      
    }
    
    targetdat <- 2^targetdat - pseudocount
    
  }
  
  if(adjustminus == TRUE){
    
    targetdat[targetdat < 0] <- 0
    refdat[refdat < 0] <- 0
    
  }
  
  
  if(ncol(targetdat) > 1){
    
    targetdatsds <- rowSD(targetdat)
    targetdat <- targetdat[targetdatsds != 0,]
    
  }
  
  if(ncol(refdat) > 1){
    
    refdatsds <- rowSD(refdat)
    refdat <- refdat[refdatsds != 0,]
    
    
  }
  
  
  sharedgenes <- intersect(row.names(targetdat), row.names(refdat))
  
  targetdat <- targetdat[sharedgenes,]
  refdat <- refdat[sharedgenes,]
  refdat<- t(refdat)
  refdat <- t(refdat)
  
  res <- list(targetdat = targetdat, 
              refdat = refdat)
  
  return(res)
  
}



singlelasso <- function(bootmethylmat, 
                        bootsolutions, 
                        threads = 1, 
                        seednum, 
                        errortype = 'min', 
                        a = 1, 
                        lowerlimits = -Inf, 
                        upperlimits = Inf, 
                        family = "mgaussian", 
                        intercept = TRUE){
  
  #library(glmnet)
  #library(doParallel)
  #library(foreach)
  
  if(threads > 1){
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    doParallel::registerDoParallel(cl)
    parallelval <- TRUE
  }else{
    parallelval <- FALSE
  }
  
  #foldnum <- max(min(floor(nrow(bootsolutions)/10), 10), 3)
  foldnum <- max(min(round(nrow(bootsolutions)/10), 10), 3)
  
  l <- a
  set.seed(seednum)
  
  
  cv <- glmnet::cv.glmnet(x = bootmethylmat, 
                          y = bootsolutions, 
                          family = family, 
                          nfold = foldnum, 
                          type.measure = "deviance", 
                          paralle = parallelval, 
                          alpha = l, 
                          lower.limits = lowerlimits, 
                          upper.limits = upperlimits, 
                          intercept = intercept)
  
  if(threads > 1){
    parallel::stopCluster(cl)
    
    unregister_dopar()
  }
  
  cvms.1ses <- c(cv$cvm[cv$lambda == cv$lambda.1se])
  cvms.mins <- c(cv$cvm[cv$lambda == cv$lambda.min])
  lambda.1ses <- c(cv$lambda.1se)
  lambda.mins <- c(cv$lambda.min)
  
  search <- data.frame(cvms.1ses = cvms.1ses, cvms.mins = cvms.mins, 
                       lambda.1ses = lambda.1ses, lambda.mins = lambda.mins, 
                       alphas = l, roundnum = seednum)
  
  
  if(errortype == 'min'){
    parameters <- search[search$cvms.mins == min(search$cvms.mins),]
    cat(paste0('Choose the elastic net mixing parameter (alpha) as ', 
               parameters$alphas, '\n'))
    cat(paste0('Choose the regularization constant (lambda) as ', 
               signif(parameters$lambda.mins, 3), '\n'))
    
    Alpha <- parameters$alphas
    Lambda <- parameters$lambda.mins
  }else{
    parameters <- search[search$cvms.1ses == min(search$cvms.1ses),]
    cat(paste0('Choose the elastic net mixing parameter (alpha) as ', 
               parameters$alphas, '\n'))
    cat(paste0('Choose the regularization constant (lambda) as ', 
               signif(parameters$lambda.1ses, 3), '\n'))
    Alpha <- parameters$alphas
    Lambda <- parameters$lambda.1ses
  }
  
  
  #Elastic net 1ses##########
  elasticnetmodel <- glmnet::glmnet(x = bootmethylmat, 
                                    y = bootsolutions, 
                                    family = family, 
                                    lambda = Lambda, 
                                    alpha = Alpha, 
                                    lower.limits = lowerlimits, 
                                    upper.limits = upperlimits, 
                                    intercept = intercept)
  
  modelcoefs <- coef(elasticnetmodel)
  
  
  modelcoef <- tryCatch({
    modelcoefs[[1]]
  }, error = function(err){
    as.matrix(modelcoefs)
  })
  
  modelcoef <- as.matrix(modelcoef)
  
  if(lowerlimits != 0){
    
    modelcoef <- modelcoef[modelcoef[,1] != 0,]
    
  }
  
  if(family == 'mgaussian'){
    
    modelcoef <- do.call(cbind, modelcoefs)
    modelcoef <- as.matrix(modelcoef)
    
    if(lowerlimits != 0){
      
      modelcoef <- modelcoef[modelcoef[,1] != 0,]
      
    }
    
    colnames(modelcoef) <- names(modelcoefs)
    
  }
  
  
  reslist <- list(elasticnetmodel = elasticnetmodel, 
                  modelcoef = modelcoef)
  
  return(reslist)
  
}


#Ridge version

singleconstraint <- function(refdat, 
                             singletarget, 
                             oriref, 
                             method = 'ridge', 
                             resscale = FALSE, 
                             threads = 1, 
                             seednum = 2022){
  
  if(method == 'svr'){
    
    #library(e1071)
    
    singletargetvector <- as.vector(singletarget)
    
    nus <- seq(0.1, 0.9, 0.1)
    
    i <- 1
    for(i in 1:length(nus)){
      
      nu <- nus[i]
      
      singlesvr <- e1071::svm(x = refdat, 
                              y = singletargetvector, 
                              scale = FALSE, 
                              type = 'nu-regression', 
                              nu = nu)
      mse <- mean(singlesvr$residuals^2)
      
      if(i == 1){
        mses <- mse
      }else{
        mses <- c(mses, mse)
      }
    }
    
    idx <- match(min(mses), mses)
    
    mod <- e1071::svm(x = refdat, 
                      y = singletargetvector, 
                      scale = FALSE, 
                      type = 'nu-regression', 
                      nu = nus[idx])
    
    modpre <- mod$fitted
    
    W = t(mod$coefs) %*% mod$SV
    estimates <- W/(sum(W))
    estimates <- estimates[1,]
    
  }else if(method == 'lm'){
    
    mod <- lsfit(x = refdat, 
                 y = singletarget, 
                 intercept = FALSE)
    modpre <- t(mod$coefficients %*% t(refdat))
    
    estimates <- mod$coefficients
    
    if(resscale == TRUE){
      estimates <- estimates/(sum(estimates))
    }
  }else if(method == 'quadprog'){
    
    Dmat <- t(refdat) %*% refdat
    dvec <- t(singletarget) %*% refdat
    Amat <- t(diag(ncol(refdat)))
    bvec <- rep(0, ncol(refdat))
    
    sc <- norm(Dmat, '2')
    
    modres <- tryCatch({
      quadprog::solve.QP(Dmat = Dmat, 
                         dvec = dvec, 
                         Amat = Amat, 
                         bvec = bvec, 
                         meq = 0, 
                         factorized = FALSE)
    }, error = function(err){
      quadprog::solve.QP(Dmat = Dmat/sc, 
                         dvec = dvec/sc, 
                         Amat = Amat, 
                         bvec = bvec, 
                         meq = 0, 
                         factorized = FALSE)
    })
    
    #If the elements in dvec and Dmat are huge, some kinds of overflow error 
    #will happen. To work around this, can scale Dmat (unfactorized) and dvec
    #Scaling Dmat and dvec is fine to do because the constraint minimizer of 
    #(-d^T b + 1/2 b^T D b ) is the same as the constrained minimizer of 
    #sc*(-d^T b + 1/2 b^T D b) for any constant sc
    
    estimates <- modres$solution
    estimates[estimates < 0] <- 0
    names(estimates) <- colnames(refdat)
    
    if(resscale == TRUE){
      estimates <- estimates/(sum(estimates))
    }
    
  }else{
    
    if(resscale == TRUE){
      upperlimits <- Inf
    }else{
      upperlimits <- 1
    }
    
    modres <- singlelasso(bootmethylmat = refdat, 
                          bootsolutions = singletarget, 
                          threads = threads, 
                          seednum = seednum, 
                          errortype = 'min', 
                          a = 0, 
                          lowerlimits = 0, 
                          upperlimits = upperlimits, 
                          family = "gaussian", 
                          intercept = FALSE)
    estimates <- modres$modelcoef
    estimates <- estimates[2:nrow(estimates), , drop = FALSE]
    estimates <- as.vector(estimates)
    estimates[estimates < 0] <- 0
    names(estimates) <- colnames(refdat)
    
    if(resscale == TRUE){
      estimates <- estimates/(sum(estimates))
    }
    
    
    
    
  }
  
  
  names(estimates) <- colnames(refdat)
  estimates <- estimates[order(-estimates)]
  
  if(as.vector(estimates)[length(estimates)] >= 0){
    estimates <- estimates[colnames(oriref)]
    names(estimates) <- colnames(oriref)
    estimates[is.na(estimates)] <- 0
    return(estimates)
  }else{
    removecelltype <- names(estimates)[length(estimates)]
    newrefdat <- refdat[,-match(removecelltype, colnames(refdat))]
    neworiref <- oriref
    
    singleconstraint(refdat = newrefdat, 
                     singletarget = singletarget, 
                     oriref = neworiref, 
                     method = method, 
                     resscale = resscale)
  }
  
}


ensemblecellcontents <- function(cellcontlist, 
                                 normweights){
  
  i <- 1
  for(i in 1:length(cellcontlist)){
    
    cellcont <- cellcontlist[[i]]
    normweight <- normweights[i,]
    
    for(j in 1:ncol(cellcont)){
      cellname <- colnames(cellcont)[j]
      normcellcont <- cellcont[,j]*normweight[j]
      normcellcont <- data.frame(cellcont = normcellcont, stringsAsFactors = FALSE)
      names(normcellcont) <- cellname
      
      if(j == 1){
        normcellconts <- normcellcont
      }else{
        normcellconts <- cbind(normcellconts, normcellcont)
      }
      
    }
    
    if(i == 1){
      finalcellconts <- normcellconts
    }else{
      finalcellconts <- finalcellconts + normcellconts
    }
    
  }
  
  return(finalcellconts)
  
}


singlerefDeconv <- function(idx = NULL, 
                            ref, 
                            targetdat, 
                            targetlogged = FALSE, 
                            modmethod = 'ridge', 
                            resscale = FALSE, 
                            adjustminus = TRUE, 
                            
                            plot = FALSE, 
                            pddat = NULL, 
                            threads = 1){
  
  if(!is.null(idx)){
    
    set.seed(idx)
    bootgenes <- sample(x = 1:nrow(ref), size = nrow(ref), replace = TRUE)
    targetdat <- targetdat[bootgenes, , drop = FALSE]
    ref <- ref[bootgenes, , drop = FALSE]
    
  }
  
  if(is.vector(targetdat)){
    
    targetdat <- as.matrix(targetdat, ncol = 1)
    colnames(targetdat) <- 'Sample1'
    
  }
  
  unifieddats <- unifydats(targetdat = targetdat, 
                           targetlogged = targetlogged, 
                           refdat = ref, 
                           adjustminus = adjustminus)
  
  targetdat <- unifieddats$targetdat
  refdat <- unifieddats$refdat
  
  if(is.vector(targetdat)){
    
    targetdat <- as.matrix(targetdat, ncol = 1)
    colnames(targetdat) <- 'Sample1'
    
  }
  
  i <- 1
  for(i in 1:ncol(targetdat)){
    
    eachtargetdat <- targetdat[, i, drop = FALSE]
    samplename <- colnames(targetdat)[i]
    
    constraintsolution <- singleconstraint(refdat = refdat, 
                                           singletarget = eachtargetdat, 
                                           oriref = refdat, 
                                           method = modmethod, 
                                           resscale = resscale, 
                                           threads = threads, 
                                           seednum = i)
    
    cellnames <- names(constraintsolution)
    constraintsolution <- matrix(constraintsolution, nrow = 1, byrow = TRUE)
    colnames(constraintsolution) <- cellnames
    row.names(constraintsolution) <- samplename
    
    if(i == 1){
      constraintsolutions <- constraintsolution
    }else{
      constraintsolutions <- rbind(constraintsolutions, constraintsolution)
    }
    
  }
  
  
  if(plot == TRUE){
    
    celldeconvplots(cellcontres = constraintsolutions, 
                    pddat = pddat, 
                    title = NULL)
    
  }
  
  
  return(constraintsolutions)
  
}

singlersquare <- function(idx, 
                          ref, 
                          targetdat, 
                          res, 
                          targetlogged = FALSE, 
                          adjustminus = TRUE){
  
  
  unifieddats <- unifydats(targetdat = targetdat, 
                           targetlogged = targetlogged, 
                           refdat = ref, 
                           adjustminus = adjustminus)
  
  targetdat <- unifieddats$targetdat
  ref <- unifieddats$refdat
  
  if(!is.null(idx)){
    
    set.seed(idx)
    bootgenes <- sample(x = 1:nrow(ref), size = nrow(ref), replace = TRUE)
    targetdat <- targetdat[bootgenes, , drop = FALSE]
    ref <- ref[bootgenes, , drop = FALSE]
    
  }
  
  
  pretarget <- ref %*% t(res)
  pccs <- cor(targetdat, pretarget)
  rsquares <- diag(pccs)^2
  
  if(!is.null(idx)){
    
    return(mean(rsquares, na.rm = TRUE))
    
  }else{
    
    return(rsquares)
  }
  
}




#'Deconvolve cell contents using reference from the same omics
#'
#'Deconvolve cell contents using reference from the same omics.
#'
#'@param ref The reference matrix recording the signature of each cell type. 
#'  Each row represents one feature, and each column is one cell type. Each 
#'  entry should be a non-log transformed value, such as TPM value for RNA 
#'  data, and beta value for DNA methylation data. Column names are cell type 
#'  names and row names are feature names.
#'@param targetdat The target cell mixture data need to be deconvolved. Should 
#'  be a matrix with each column representing one sample and each row for one 
#'  feature. Row names are feature names and column names are sample IDs. If 
#'  the reference matrix is generated with the function \code{scRef}, and both 
#'  the scRNA-seq and the target cell mixture data were transferred to it, the 
#'  result reference matrix can be transferred to \code{ref} and the adjusted 
#'  target data returned by \code{scRef} can be transferred to this parameter.
#'@param targetlogged Whether the feature values in \code{targetdat} are log2 
#'  transformed values or not. 
#'@param resscale For each sample, whether its cell contents result should be 
#'  scaled so that the sum of different cell types is 1. Default is FALSE.
#'@param plot Whether generate a box plot and heatmaps for the cell contents 
#'  deconvolved. Default is FALSE.
#'@param pddat If set \code{plot} as TRUE, this parameter can be used to show 
#'  the sample group information, so that the box plot generated will also 
#'  compare the group difference for each cell type, and heatmaps with this 
#'  comparison will also be generated. It should be a data frame recording the 
#'  sample groups, and must include 2 columns. One is named as "sampleid", 
#'  recording the sample IDs same as the column names of \code{targetdat}, the 
#'  other is "Samplegroup", recording the sample group to which each sample 
#'  belongs. It can also be NULL, meaning all the samples are from the same 
#'  group, and in this case, only a box plot showing the cell content for each 
#'  cell type will be made when \code{plot} is set as TRUE.
#'@param threads Number of threads need to be used to do the computation. Its 
#'  default value is 1.
#'@return A matrix recording the cell composition result for the samples.
#'@export
refDeconv <- function(ref, 
                      targetdat, 
                      targetlogged = FALSE, 
                      resscale = FALSE, 
                      
                      plot = FALSE, 
                      pddat = NULL, 
                      threads = 1){
  
  learnernum <- 1
  removeoutliers <- FALSE
  adjustminus <- TRUE
  
  if(removeoutliers == TRUE){
    
    preres <- singlerefDeconv(idx = NULL, 
                              ref = ref, 
                              targetdat = targetdat, 
                              targetlogged = targetlogged, 
                              modmethod = 'lm', 
                              resscale = resscale, 
                              adjustminus = adjustminus, 
                              plot = FALSE, 
                              pddat = NULL, 
                              threads = threads)
    
    prersquares <- singlersquare(idx = NULL, 
                                 ref = ref, 
                                 targetdat = targetdat, 
                                 res = preres, 
                                 targetlogged = targetlogged, 
                                 adjustminus = adjustminus)
    
    prersquares <- prersquares[!is.na(prersquares)]
    prersquares <- prersquares[prersquares > 0.5]
    
    if(length(prersquares) == 0){
      
      return(NULL)
      
    }else{
      targetdat <- targetdat[,names(prersquares), drop = FALSE]
    }
    
  }
  
  
  reslist <- list()
  rsquaremeans <- list()
  n <- 1
  i <- 1
  
  if(learnernum == 1){
    
    constraintsolutions <- singlerefDeconv(idx = NULL, 
                                           ref = ref, 
                                           targetdat = targetdat, 
                                           targetlogged = targetlogged, 
                                           modmethod = 'lm', 
                                           resscale = resscale, 
                                           adjustminus = adjustminus, 
                                           plot = FALSE, 
                                           pddat = NULL, 
                                           threads = threads)
    
    constraintsolutions <- as.matrix(constraintsolutions)
    
  }else{
    
    while(n <= learnernum){
      
      res <- singlerefDeconv(idx = i, 
                             ref = ref, 
                             targetdat = targetdat, 
                             targetlogged = targetlogged, 
                             modmethod = 'lm', 
                             resscale = resscale, 
                             adjustminus = adjustminus, 
                             plot = FALSE, 
                             pddat = NULL, 
                             threads = threads)
      
      rsquaremean <- singlersquare(idx = i, 
                                   ref = ref, 
                                   targetdat = targetdat, 
                                   res = res, 
                                   targetlogged = targetlogged, 
                                   adjustminus = adjustminus)
      
      
      i <- i + 1
      
      if(is.na(rsquaremean)){
        cat(paste0('Base learner ', (i - 1), ' was dropped because of NA value\n'))
        
        next()
      }
      
      if(rsquaremean <= 0.5){
        
        cat(paste0('Base learner ', (i - 1), ' was dropped because of an R square <= 0.5\n'))
        
        next()
      }
      
      
      reslist[[n]] <- res
      rsquaremeans[[n]] <- rsquaremean
      
      n <- n + 1
      
      cat(paste0('Base learner ', (n - 1), ' has been finished\n'))
      
    }
    
    rsquaremeans <- unlist(rsquaremeans)
    
    weights <- 0.5*log(rsquaremeans/(1 - rsquaremeans))
    
    normweights <- weights/sum(weights)
    
    normweightsmat <- replicate(ncol(targetdat), normweights)
    
    
    
    constraintsolutions <- ensemblecellcontents(cellcontlist = reslist, 
                                                normweights = normweightsmat)
    constraintsolutions <- as.matrix(constraintsolutions)
    
  }
  
  
  if(plot == TRUE){
    
    celldeconvplots(cellcontres = constraintsolutions, 
                    pddat = pddat, 
                    title = NULL)
    
  }
  
  
  return(constraintsolutions)
  
}



#Generate RNA reference using Seurat scRNA-seq object and then deconvolve 
#target RNA data


#'Deconvolve bulk RNA data using scRNA-seq data
#'
#'Deconvolve bulk RNA-seq or bulk RNA microarray data using scRNA-seq data.
#'
#'@param Seuratobj An object of class Seurat generated with the \code{Seurat} 
#'  R package from scRNA-seq data, should contain read count data, normalized 
#'  data, and cell meta data. The meta data should contain a column recording 
#'  the cell type name of each cell.
#'@param targetcelltypes The cell types whose content need to be deconvolved. 
#'  If NULL, all the cell types included in \code{Seuratobj} will be included. 
#'  Default is NULL.
#'@param celltypecolname In the "meta.data" slot of \code{Seuratobj}, which 
#'  column records the cell type information for each cell and the name of 
#'  this column should be transferred to this parameter. Default value is 
#'  "annotation".
#'@param samplebalance At the beginning of making the cell reference matrix, 
#'  the scRNA-seq cell counts contained in \code{Seuratobj} will be sampled 
#'  and used to generate 100 pseudo-bulk RNA-seq sample, for each cell type, 
#'  and during generating them, the number of single cells can be sampled is 
#'  always different for each cell type. If want to adjust this bias and make 
#'  the single cell numbers used to make pseudo-bulk RNA-seq data same for 
#'  different cell types, set this parameter as TRUE. Then, the cell types 
#'  with too many candidate cells will be down-sampled while the ones with 
#'  much fewer cells will be over-sampled. The down-sampling is performed with 
#'  bootstrapping, and the over-sampling is conducted with SMOTE (Synthetic 
#'  Minority Over-sampling Technique). This is a time-consuming step and the 
#'  default value of this parameter is FALSE, so that no such adjustment will 
#'  be done during sampling the pseudo-bulk samples.
#'@param geneversion To calculate the TPM value of the genes in the reference 
#'  matrix, the effective length of the genes will be needed. This parameter 
#'  is used to define from which genome version the effective gene length will 
#'  be extracted. For human genes, "hg19" or "hg38" can be used, for mouse, 
#'  "mm10" can be used. Default is "hg19". 
#'@param genekey The type of the gene IDs used in the \code{Seuratobj}, it is 
#'  "SYMBOL" in most cases, and the default value of this parameter is also 
#'  "SYMBOL", but sometimes it may be "ENTREZID", "ENSEMBL", or other types.
#'@param targetdat The target cell mixture gene expression data need to be 
#'  deconvolved. Should be a matrix with each column representing one sample 
#'  and each row representing one gene. The gene ID type here should be the 
#'  same as that transferred to the parameter \code{genekey}. Row names are 
#'  gene IDs and column names are sample IDs.
#'@param targetlogged Whether the gene expression values in \code{targetdat} 
#'  are log2 transformed values or not. 
#'@param manualmarkerlist During making the reference matrix from scRNA-seq 
#'  data, for each cell type, the genes specially expressed in it with a high 
#'  level will be deemed  as markers and used to generate the reference, but 
#'  it cannot be ensured that some known classical markers can be selected, 
#'  and so if want to make sure these markers can be used for the reference, a 
#'  list can be used as an input to this parameter, with its element names as 
#'  the cell type names and the elements as vectors with the gene IDs of these 
#'  classical markers. It should be noted that before the final reference is 
#'  determined, all the marker genes need to go through several filter steps, 
#'  such as extremely highly expressed genes and collinearity contributing 
#'  genes removal, to improve the reference quality, so that the classical 
#'  genes provided via this parameter will be definitely used for reference 
#'  generation, but may also be filtered out before the final one is made. The 
#'  default value of this parameter is NULL.
#'@param markerremovecutoff During the reference matrix generation from the 
#'  scRNA-seq data, the gene expression values in the \code{targetdat} matrix 
#'  are used to calculate the correlation with the scRNA-seq selected markers 
#'  in this \code{targetdat} matrix and the ones with a high correlation to 
#'  the first principle component of these marker genes will also be used to 
#'  make the reference. The cutoff of the correlation coefficient is set by 
#'  this parameter and the default value is 0.6.
#'@param minrefgenenum Because the genes to generate the reference matrix need 
#'  to go through several filter steps and in some cases, only a small number 
#'  of them can fulfill all the filter conditions, which makes the gene number 
#'  in the reference is very small and then influences the next deconvolution. 
#'  To avoid this extreme case, a cutoff for the reference gene number need to 
#'  be defined here, so that once the gene number in the reference has been 
#'  filtered to this level, the filter process will be ended to guarantee the 
#'  gene number of the reference. This parameter is used to set this cutoff, 
#'  and its default value is 500.
#'@param saveref Whether need to save the finally generated reference matrix, 
#'  and the adjusted cell mixture matrix to be deconvolved as rds files in the 
#'  working directory automatically. Default is FALSE.
#'@param refcutoff To improve the robustness of the deconvolution result, some 
#'  extremely highly expressed genes in the reference need to be filtered out 
#'  due to their large variance. This cutoff is used to set the percent of 
#'  genes can be kept in the reference while the other genes with a higher 
#'  expression level will be filtered. The default value is 0.95, meaning the 
#'  top 5% most highly expressed genes will be removed from the reference.
#'@param refadjustcutoff For some similar cell types, their gene expressions 
#'  in the reference matrix are highly correlated, which makes the downstream 
#'  deconvolution difficult. To relive this problem, for each similar cell 
#'  pair, some genes largely contributing to their correlation will be found 
#'  and removed, so that their correlation in the reference can be reduced. 
#'  This parameter is used to set the cutoff of the cell pair correlation, and 
#'  if a cell pair has a Pearson correlation coefficient greater than it, the 
#'  contributing gene filter process will be used to reduce the coefficient 
#'  until it becomes smaller than this value. The default is 0.4.
#'@param resscale For each sample, whether its cell contents result should be 
#'  scaled so that the sum of different cell types is 1. Default is FALSE.
#'@param plot Whether generate a box plot and heatmaps for the cell contents 
#'  deconvolved. Default is FALSE.
#'@param pddat If set \code{plot} as TRUE, this parameter can be used to show 
#'  the sample group information, so that the box plot generated will also 
#'  compare the group difference for each cell type, and heatmaps with this 
#'  comparison will also be generated. It should be a data frame recording the 
#'  sample groups, and must include 2 columns. One is named as "sampleid", 
#'  recording the sample IDs same as the column names of \code{targetdat}, the 
#'  other is "Samplegroup", recording the sample group to which each sample 
#'  belongs. It can also be NULL, meaning all the samples are from the same 
#'  group, and in this case, only a box plot showing the cell content for each 
#'  cell type will be made when \code{plot} is set as TRUE.
#'@param threads Number of threads need to be used to do the computation. Its 
#'  default value is 1.
#'@return A list containing the generated RNA reference, the adjusted target 
#'  data to be deconvolved, and the cell deconvolution result for the samples. 
#'  The gene values in the adjusted target data are non-log transformed ones. 
#'@export
scDeconv <- function(Seuratobj,  
                     targetcelltypes = NULL, 
                     celltypecolname = "annotation", 
                     samplebalance = FALSE, 
                     geneversion = 'hg19', 
                     genekey = "SYMBOL", 
                     targetdat = NULL, 
                     targetlogged = FALSE, 
                     manualmarkerlist = NULL, 
                     markerremovecutoff = 0.6, 
                     minrefgenenum = 500, 
                     saveref = FALSE, 
                     refcutoff = 0.95, 
                     refadjustcutoff = 0.4, 
                     resscale = FALSE, 
                     
                     plot = FALSE, 
                     pddat = NULL, 
                     
                     threads = 1){
  
  learnernum <- 1
  removeoutliers <- FALSE
  adjustminus <- TRUE
  
  if(is.null(targetdat)){
    cat('A data matrix to be deconvolved is required by the parameter `targetdat`.\n')
    return(NULL)
  }
  
  refres <- scRef(Seuratobj = Seuratobj,  
                  targetcelltypes = targetcelltypes,  
                  celltypecolname = celltypecolname,  
                  pseudobulknum = 100, 
                  samplebalance = samplebalance, 
                  pseudobulkpercent = 0.9, 
                  geneversion = geneversion, 
                  genekey = genekey, 
                  targetdat = targetdat, 
                  targetlogged = targetlogged, 
                  manualmarkerlist = manualmarkerlist, 
                  markerremovecutoff = markerremovecutoff, 
                  minrefgenenum = minrefgenenum, 
                  savefile = saveref, 
                  threads = threads, 
                  cutoff = refcutoff, 
                  adjustcutoff = refadjustcutoff)
  
  scref <- refres$ref
  targetnolog <- refres$targetnolog
  
  deconvres <- refDeconv(ref = scref, 
                         targetdat = targetnolog, 
                         targetlogged = FALSE, 
                         resscale = resscale, 
                         plot = plot, 
                         pddat = pddat, 
                         threads = threads)
  
  res <- list(scref = scref, 
              targetnolog = targetnolog, 
              deconvres = deconvres)
  
  return(res)
  
}



