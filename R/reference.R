
#Make reference#####

#Some parallel computing going on in the background that is not getting cleaned up
#fully between runs can cause Error in summary.connection(connection) : invalid
#connection. The following function is needed to be called to fix this error.
unregister_dopar <- function(){

  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)

}


#Output count table and cell meta data from Seurat object, as well as
#cell markers via comparing one cell type with others
processSeuratobj <- function(Seuratobj,
                             targetcells = NULL,
                             celltypecol = 'annotation'){

  #library(Seurat)

  #Extract count table
  counts <- Seurat::GetAssayData(Seuratobj, slot = 'counts')

  #Extract and organize cell meta data
  meta <- Seuratobj[[]]

  meta <- meta[c('orig.ident', 'nCount_RNA', 'nFeature_RNA',
                 celltypecol)]
  names(meta)[ncol(meta)] <- 'celltype'

  if(!is.null(targetcells)){
    targetcells <- targetcells[targetcells %in% meta$celltype]
    if(length(targetcells) > 0){
      submeta <- subset(meta, celltype %in% targetcells)
      #submeta$celltype <- factor(submeta$celltype,
      #                           levels = targetcells, ordered = TRUE)
      #submeta <- submeta[order(submeta$celltype),]
      #submeta$celltype <- as.character(submeta$celltype)
    }else{
      submeta <- meta
      #submeta <- submeta[order(submeta$celltype),]
    }

  }else{
    submeta <- meta
    #submeta <- submeta[order(submeta$celltype),]
  }

  #rm(meta)

  #Organize count table according to the organized cell meta data
  #cellidents <- names(Seurat::Idents(object = Seuratobj))
  #cellidents <- cellidents[cellidents %in% row.names(submeta)]
  #submeta <- submeta[cellidents,]

  subcounts <- counts[,row.names(submeta)]
  subcounts <- t(as.matrix(subcounts))

  subcounts <- t(subcounts)

  #Subset Seurat object!###
  #Seuratobj <- Seuratobj[,Seuratobj[['annotation']][,1] %in% unique(submeta$celltype)]
  ####

  #Find differentially expressed cell markers
  Seurat::Idents(object = Seuratobj) <- submeta$celltype

  celltypes <- unique(submeta$celltype)



  cellmarkers <- list()

  i <- 1
  for(i in 1:length(celltypes)){

    celltype <- celltypes[i]

    cellmarkernum <- 0
    minpct <- 0.25
    logfcthreshold <- 1

    while(cellmarkernum < 10 &
          minpct >= 0.1 &
          logfcthreshold >= 0.25){

      cellmarker <- Seurat::FindMarkers(Seuratobj,
                                        ident.1 = celltype,
                                        min.pct = minpct,
                                        only.pos = TRUE,
                                        logfc.threshold = logfcthreshold)

      if(nrow(cellmarker) > 0){

        if(sum(cellmarker$p_val_adj < 0.05) < 10 &
           sum(cellmarker$p_val < 0.01) >= 10){

          cellmarker <- subset(cellmarker, p_val < 0.01)

        }else{
          cellmarker <- subset(cellmarker, p_val_adj < 0.05)
        }

      }


      cellmarkernum <- nrow(cellmarker)
      minpct <- minpct - 0.05
      logfcthreshold <- logfcthreshold - 0.25

    }

    if(cellmarkernum > 0){

      cellmarker$celltype <- celltype
      cellmarker$gene <- row.names(cellmarker)

      if('avg_log2FC' %in% names(cellmarker)){

        cellmarker <- cellmarker[c('p_val', 'avg_log2FC',
                                   'pct.1', 'pct.2', 'p_val_adj',
                                   'celltype', 'gene')]

      }else if('avg_logFC' %in% names(cellmarker)){

        cellmarker <- cellmarker[c('p_val', 'avg_logFC',
                                   'pct.1', 'pct.2', 'p_val_adj',
                                   'celltype', 'gene')]
        names(cellmarker)[2] <- 'avg_log2FC'


      }

      row.names(cellmarker) <- 1:nrow(cellmarker)



    }else{

      cellmarker <- data.frame(p_val = numeric(),
                               avg_log2FC = numeric(),
                               pct.1 = numeric(),
                               pct.2 = numeric(),
                               p_val_adj = numeric(),
                               celltype = character(),
                               gene = character(),
                               stringsAsFactors = FALSE)


    }

    cellmarkers[[i]] <- cellmarker
    names(cellmarkers)[i] <- celltype

    #print(i)
  }

  cellmarkerdat <- do.call(rbind, cellmarkers)
  cellmarkerdat <- unique(cellmarkerdat)
  row.names(cellmarkerdat) <- 1:nrow(cellmarkerdat)

  res <- list(counts = subcounts,
              meta = submeta,
              cellmarkers = cellmarkerdat)

  return(res)

}

#Generate preliminary reference via sampling
#nround is the round number of sampling
#In each round, for each of the cell types, samplepercent cells
#will be randomly selected and merged together as a psudo-bulk
#cell sequencing sample
makeref.cv <- function(scmatrix = Seuratobjlist$counts,
                       metainfo = Seuratobjlist$meta,
                       nround = 10,
                       samplepercent = 0.9,
                       balance = FALSE,
                       threads = 1){

  scmatrix <- t(scmatrix)
  cellid <- row.names(scmatrix)
  scmatrix <- as.data.frame(scmatrix)
  scmatrix$cellid <- cellid

  metainfo$cellid <- row.names(metainfo)
  cellidmapping <- metainfo[c('cellid', 'celltype')]


  scmatrix <- merge(cellidmapping, scmatrix, by = c('cellid'))

  cellidmapping <- scmatrix[c('cellid', 'celltype')]

  celltypes <- unique(cellidmapping$celltype)

  row.names(scmatrix) <- scmatrix$cellid



  mergeblock <- function(block){
    block <- block[-1]

    skipminus <- function(col){
      res <- sum(col)
      #Calculate sum, not mean
      return(res)
    }

    res <- apply(X = block, MARGIN = 2, FUN = skipminus)
    return(res)
  }

  singleroundsampling <- function(celltypes = celltypes,
                                  cellidmapping = cellidmapping,
                                  i = i,
                                  scmatrix = scmatrix,
                                  balance = balance,
                                  samplepercent = samplepercent){

    for(j in 1:length(celltypes)){

      CellType <- celltypes[j]
      sub <- subset(cellidmapping, celltype == CellType)

      if(balance == TRUE){

        targetnum <- ceiling(nrow(cellidmapping)/length(celltypes))
        subnum <- nrow(sub)

        if(subnum >= targetnum | subnum == 1){

          cellids <- sub$cellid

          set.seed(i)
          samples <- sample(x = cellids,
                            size = targetnum,
                            replace = TRUE)

          scmatrixsub <- scmatrix[samples,]


        }else{

          sub <- subset(scmatrix, celltype == CellType)
          subvar <- sub[-c(1, 2)]

          subvar <- as.matrix(subvar)

          dists <- dist(subvar, method = 'euclidean', diag = FALSE)
          distmin <- min(dists)
          dists <- as.matrix(dists)

          pairindces <- which(dists == distmin, arr.ind = TRUE)[1,]
          sample1 <- row.names(dists)[pairindces[1]]
          sample2 <- row.names(dists)[pairindces[2]]
          sample1var <- subvar[sample1,]
          sample2var <- subvar[sample2,]

          diff1 <- sample1var - sample2var
          diff1 <- diff1/2
          diff2 <- -diff1

          diffs <- rbind(diff1, diff2)
          row.names(diffs) <- 1:2

          set.seed(i)
          SMOTEressampleididx <- sample(x = 1:nrow(sub),
                                        size = targetnum,
                                        replace = TRUE)
          SMOTEressampleid <- sub$cellid[SMOTEressampleididx]

          set.seed(i)
          SMOTEdiffid <- sample(x = c(1, 2),
                                size = targetnum,
                                replace = TRUE)

          set.seed(i)
          zetas <- runif(n = targetnum, min = 0, max = 1)
          zetas[!duplicated(SMOTEressampleididx)] <- 0

          SMOTEsamplevar <- subvar[SMOTEressampleid,]
          SMOTEdiff <- diffs[SMOTEdiffid,]
          SMOTEdiff <- SMOTEdiff * zetas

          SMOTEsamplevar <- SMOTEsamplevar + SMOTEdiff


          SMOTEsub <- sub[SMOTEressampleididx,c(1,2)]

          suffix <- rep('', nrow(SMOTEsub))
          suffix[grepl(pattern = '\\.', x = row.names(SMOTEsub))] <-
            substring(row.names(SMOTEsub),
                      regexpr('\\.', row.names(SMOTEsub)))[grepl(pattern = '\\.',
                                                                 x = row.names(SMOTEsub))]

          SMOTEsub$cellid <- paste0(SMOTEsub$cellid, suffix)
          row.names(SMOTEsub) <- 1:nrow(SMOTEsub)
          row.names(SMOTEsamplevar) <- SMOTEsub$cellid

          scmatrixsub <- cbind(SMOTEsub, SMOTEsamplevar)

        }

        if(j == 1){
          scmatrixsubs <- scmatrixsub
        }else{
          scmatrixsubs <- rbind(scmatrixsubs, scmatrixsub)
        }

      }else{

        cellids <- sub$cellid

        samplesize <- max(1, round(length(cellids)*samplepercent))

        set.seed(i)
        samples <- sample(x = cellids,
                          size = samplesize,
                          replace = FALSE)

        scmatrixsub <- scmatrix[samples,]

        if(j == 1){
          scmatrixsubs <- scmatrixsub
        }else{
          scmatrixsubs <- rbind(scmatrixsubs, scmatrixsub)
        }



      }

    }

    scmatrixsubs <- scmatrixsubs[-1]

    #library(plyr)
    #tmp <- ddply(.data = scmatrixsub, .variables = c('cell_type'), .fun = mergeblock)
    tmp <- plyr::ddply(.data = scmatrixsubs, .variables = c('celltype'), .fun = mergeblock)

    #row.names(tmp) <- tmp$cell_type
    row.names(tmp) <- tmp$celltype
    tmp <- tmp[-1]
    tmp <- t(tmp)

    colnames(tmp) <- paste0(colnames(tmp), '_', i)

    return(tmp)

  }


  if(threads == 1){

    for(i in 1:nround){

      tmp <- singleroundsampling(celltypes = celltypes,
                                 cellidmapping = cellidmapping,
                                 i = i,
                                 scmatrix = scmatrix,
                                 balance = balance,
                                 samplepercent = samplepercent)

      if(i == 1){
        tmps <- tmp
      }else{
        tmps <- cbind(tmps, tmp)
      }
    }


  }else{

    iseqs <- 1:nround

    #library(doParallel)

    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))

    doParallel::registerDoParallel(cl)

    #date()
    `%dopar%` <- foreach::`%dopar%`
    tmplist <- foreach::foreach(i = iseqs,
                                #.export = ls(name = globalenv())) %dopar% {
                                .export = NULL) %dopar% {
                                  singleroundsampling(celltypes = celltypes,
                                                      cellidmapping = cellidmapping,
                                                      i,
                                                      scmatrix = scmatrix,
                                                      balance = balance,
                                                      samplepercent = samplepercent)
                                }

    parallel::stopCluster(cl)

    unregister_dopar()

    tmps <- do.call(cbind, tmplist)

    rm(tmplist)

  }

  tmps <- round(tmps)
  tmps[tmps < 0] <- 0

  return(tmps)

}



#'Synthesize pseudo-bulk RNA-seq data from scRNA-seq data
#'
#'Synthesize pseudo-bulk RNA-seq data for RNA reference generation.
#'
#'@param Seuratobj An object of class Seurat generated with the \code{Seurat}
#'  R package from scRNA-seq data, should contain read count data, normalized
#'  data, and cell meta data. The meta data should contain a column recording
#'  the cell type name of each cell.
#'@param targetcelltypes The cell types in \code{Seuratobj} whose content need
#'  to be deconvolved via \code{scDeconv} package. If NULL, all the cell types
#'  included in it will be included. Default is NULL.
#'@param celltypecolname In the "meta.data" slot of \code{Seuratobj}, which
#'  column records the cell type information for each cell and the name of
#'  this column should be transferred to this parameter. Default value is
#'  "annotation".
#'@param pseudobulknum The scRNA-seq cell counts contained in \code{Seuratobj}
#'  will be sampled and used to generate some pseudo-bulk RNA-seq samples, for
#'  each cell type. The parameter \code{pseudobulknum} here defines how many
#'  pseudo-bulk RNA-seq data for each cell type need to be generated. Default
#'  is 10.
#'@param samplebalance During generating the pseudo-bulk RNA-seq data, the
#'  number of single cells can be sampled is always different for each cell
#'  type. If want to adjust this bias and make the single cell numbers used to
#'  make pseudo-bulk RNA-seq data same for different cell types, set this
#'  parameter as TRUE. Then, the cell types with too many candidate cells will
#'  be down-sampled while the ones with much fewer cells will be over-sampled.
#'  The down-sampling is performed using bootstrapping, and the over-sampling
#'  is conducted with SMOTE (Synthetic Minority Over-sampling Technique). This
#'  is a time-consuming step and the default value of this parameter is FALSE.
#'@param pseudobulkpercent If the parameter \code{samplebalance} is FALSE, for
#'  the pseudo-bulk sampling for each cell type, a percent of single cells for
#'  each cell type will be randomly sampled and this parameter is used to set
#'  this percent value and should be a number between 0 and 1, but if the
#'  parameter \code{samplebalance} is set as TRUE, bootstrapping and SMOTE
#'  will be performed to do the sampling and this parameter will be omitted.
#'@param threads Number of threads need to be used. Its default value is 1.
#'@param savefile Whether need to save the generated pseudo-bulk matrix as an
#'  rds file in the working directory automatically. Default is FALSE.
#'@return A pseudo-bulk RNA-seq matrix with pseudo-bulk samples as columns and
#'  genes as features. The gene values in this matrix are pseudo-bulk RNA-seq
#'  read counts. This matrix can be transferred to the functions \code{scRef},
#'  \code{scDeconv}, or \code{epDeconv}. Their parameter \code{pseudobulkdat}
#'  can accept this matrix, so that they can skip their own pseudo-bulk data
#'  synthesis step and directly use this matrix as their pseudo-bulk data to
#'  further generate the RNA deconvolution reference. Because if the scRNA-seq
#'  dataset need to be converted to the RNA referece is large, generating the
#'  pseudo-bulk data can be time-consuming and if the scRNA-seq data need to
#'  be repeatedly used to deconvolve different datasets, to avoid repeating
#'  this pseudo-bulk data generation process, this function can be used to
#'  synthesize and save the data in advance, then the data can be repeatedly
#'  used and the synthesis step can always be skipped.
#'@export
prepseudobulk <- function(Seuratobj,
                          targetcelltypes = NULL,
                          celltypecolname = "annotation",
                          pseudobulknum = 10,
                          samplebalance = FALSE,
                          pseudobulkpercent = 0.9,
                          threads = 1,
                          savefile = FALSE){


  Seuratobjlist <- processSeuratobj(Seuratobj = Seuratobj,
                                    targetcells = targetcelltypes,
                                    celltypecol = celltypecolname)

  if(!is.null(targetcelltypes)){

    #Extract and organize cell meta data
    meta <- Seuratobj[[]]
    meta <- meta[c('orig.ident', 'nCount_RNA', 'nFeature_RNA',
                   celltypecolname)]
    names(meta)[ncol(meta)] <- 'celltype'

    targetcelltypenum <- length(unique(targetcelltypes[targetcelltypes %in% meta$celltype]))

    if(length(unique(Seuratobjlist$cellmarkers$celltype)) < ceiling(targetcelltypenum*0.5)){

      Seuratobjlist <- processSeuratobj(Seuratobj = Seuratobj,
                                        targetcells = NULL,
                                        celltypecol = celltypecolname)


    }

  }

  refcounts <- makeref.cv(scmatrix = Seuratobjlist$counts,
                          metainfo = Seuratobjlist$meta,
                          nround = pseudobulknum,
                          samplepercent = pseudobulkpercent,
                          balance = samplebalance,
                          threads = threads)

  tag <- ''

  if(savefile == TRUE){

    if(samplebalance == TRUE){
      tag <- '.balance'
    }else{
      tag <- '.nobalance'
    }

    stamp <- Sys.time()
    stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
    stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
    stamp <- paste0('.', stamp)

    if(is.null(pseudobulkdat)){
      saveRDS(refcounts, file = paste0('orirefcounts', stamp, tag, '.rds'))
    }

  }

  return(refcounts)

}




#Remove the '.' in gene names and if after removing, some original
#different gene names become the same, merge them together by
#getting their means
orggeneids <- function(oriref){

  oriref <- as.data.frame(oriref)
  geneids <- row.names(oriref)
  geneids <- gsub(pattern = '\\..*$', replacement = '',
                  x = geneids)
  oriref$geneids <- geneids

  dupgeneids <- unique(geneids[duplicated(geneids)])
  oriref <- oriref[colnames(oriref)[c(ncol(oriref),
                                      1:(ncol(oriref) - 1))]]

  uniref <- subset(oriref, !(geneids %in% dupgeneids))
  dupref <- subset(oriref, geneids %in% dupgeneids)

  mergeblock <- function(block){
    block <- block[-1]

    skipminus <- function(col){
      res <- mean(col)
      #Calculate mean, not sum
      return(res)
    }

    res <- apply(X = block, MARGIN = 2, FUN = skipminus)
    return(res)
  }

  #library(plyr)
  tmp <- plyr::ddply(.data = dupref,
                     .variables = c('geneids'),
                     .fun = mergeblock)

  row.names(tmp) <- tmp$geneids
  tmp <- tmp[-1]

  rm(dupref)

  uniref <- uniref[-1]

  oriref <- rbind(uniref, tmp)
  oriref <- as.matrix(oriref)

  return(oriref)

}

#Remove genes with unusual gene symbols
#1) without an entrez id
#2) with multiple entrez ids
#3) with a shared entrez id
#It is because the following TPM calculation step needs gene entrez id
#to get the exon length of the genes, if an entrez id is NA,
#the gene exon length cannot be gotten, and if
#multiple mapping exists between gene symbols and entrez ids,
#the gene counts in the count matrix of one gene need to be decomposed
#to different entrez ids, which cannot be done, or
#the gene counts in the count matrix of different genes need to be
#merged to one shared entrez id, which will lead different genes with
#different gene symbols in the original dataset be merged, hence,
#all the genes without entrez id or with a multiple gene symbol-
#entrez id relationship will be removed.
orggenes <- function(oriref,
                     version,
                     genekey){

  #library(org.Hs.eg.db)
  #library(AnnotationDbi)

  oriref <- orggeneids(oriref = oriref)

  if(genekey == 'ENTREZID'){
    return(oriref)
  }else{

    if(version == 'mm10'){

      #Use the database org.Mm.eg.db for mm10 annotation,
      #not org.Mmu.eg.db!
      #org.Mmu.eg.db is for monkey!
      refids <- AnnotationDbi::select(x = org.Mm.eg.db::org.Mm.eg.db,
                                      keys = row.names(oriref),
                                      columns = 'ENTREZID',
                                      keytype = genekey)

    }else{

      refids <- AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db,
                                      keys = row.names(oriref),
                                      columns = 'ENTREZID',
                                      keytype = genekey)

    }

    #Remove genes without entrez id
    regs <- refids[complete.cases(refids),]

    finals <- regs

    names(finals)[1] <- genekey

    if(length(unique(finals[,1])) < nrow(finals)){

      dupkeys <- unique(finals[,genekey][duplicated(finals[,genekey])])
      finals <- finals[!(finals[,1] %in% dupkeys),]

    }

    if(length(unique(finals$ENTREZID)) < nrow(finals)){

      dupentrezid <- unique(finals$ENTREZID[duplicated(finals$ENTREZID)])
      finals <- subset(finals, !(ENTREZID %in% dupentrezid))
    }


    newref <- oriref[finals[,1],]

    res <- list(newref = newref, finals = finals)

  }

  return(res)

}


#Get the effective gene length for TPM calculating
gettpmlength <- function(version = 'hg19'){

  if(version == 'hg38'){

    exons.db <- GenomicFeatures::exonsBy(
      TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
      by = 'gene'
    )

  }else if(version == 'mm10'){

    exons.db <- GenomicFeatures::exonsBy(
      TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
      by = 'gene'
    )

  }else{

    exons.db <- GenomicFeatures::exonsBy(
      TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
      by = 'gene'
    )

  }

  tpmlength <- sapply(names(exons.db),
                      function(eg){
                        exons <- exons.db[[eg]]
                        exons <- GenomicRanges::reduce(exons)
                        print(eg)
                        return(sum(GenomicRanges::width(exons)))
                      }
  )

  return(tpmlength)

}

#tpmlength.hg19 <- gettpmlength(version = 'hg19')

#saveRDS(tpmlength.hg19, 'tpmlength.hg19.rds')

#tpmlength.hg38 <- gettpmlength(version = 'hg38')

#saveRDS(tpmlength.hg38, 'tpmlength.hg38.rds')

#tpmlength.mm10 <- gettpmlength(version = 'mm10')

#saveRDS(tpmlength.mm10, 'tpmlength.mm10.rds')


caltpms <- function(countmat = refcounts,
                    version = 'hg19',
                    genekey = 'SYMBOL',
                    TPMorFPKM = 'TPM'){

  #library(scater)

  #library(TxDb.Hsapiens.UCSC.hg19.knownGene)

  tpmlengthdat <- paste0('C:/Users/liuy47/Desktop/Transfer/codestransfer/deconv/files/tpmlength.',
                         version, '.rds')
  tpmlength <- readRDS(tpmlengthdat)

  #tpmlengthdat <- paste0('tpmlength.', version)
  #tpmlength <- get(tpmlengthdat)


  orgmat <- orggenes(oriref = countmat,
                     version = version,
                     genekey = genekey)

  if(genekey == 'ENTREZID'){

    countmat <- orgmat

    sharedentrezids <- intersect(row.names(countmat), names(tpmlength))
    tpmlength <- tpmlength[sharedentrezids]

  }else{

    countmat <- orgmat$newref
    entrezids <- orgmat$finals

    sharedentrezids <- intersect(entrezids$ENTREZID, names(tpmlength))
    sharedentrezids <- subset(entrezids, ENTREZID %in% sharedentrezids)
    tpmlength <- tpmlength[sharedentrezids$ENTREZID]
    sharedentrezids <- sharedentrezids[,1]
  }


  countmat <- countmat[sharedentrezids,]
  names(tpmlength) <- sharedentrezids

  sharedkeys <- intersect(row.names(countmat), names(tpmlength))
  tpmlength <- tpmlength[sharedkeys]
  countmat <- countmat[sharedkeys,]

  if(TPMorFPKM == 'FPKM'){

    reftpms <- scater::calculateFPKM(x = countmat,
                                     lengths = tpmlength)

  }else{

    reftpms <- scater::calculateTPM(x = countmat,
                                    lengths = tpmlength)
  }

  return(reftpms)

}

removenonsds <- function(dat){

  if(ncol(dat) == 1){
    return(dat)
  }

  featuresds <- apply(X = dat, MARGIN = 1, FUN = sd)
  features <- names(featuresds[featuresds != 0])

  dat <- dat[features,]

  return(dat)

}


genepc1 <- function(expdata = expdata,
                    celltypemarkers = celltypemarkers){

  expcelltypemarkers <- intersect(colnames(expdata),
                                  celltypemarkers)
  expsub <- expdata[, expcelltypemarkers, drop = FALSE]


  pcs <- prcomp(expsub, scale. = TRUE)
  pc1 <- pcs$x[,1]

  #Note, The signs of the columns of the rotation matrix are
  #arbitrary, and so may differ between different programs for
  #PCA, and even between different builds of R.
  #Hence, need to check the correlation between pc1 and the
  #row means of expsub, if it is less than 0, it means the
  #sign of pc1 is reversed, and need to reverse it back
  pc1cor <- cor(rowMeans(expsub), pc1)
  if(pc1cor < 0){
    pc1 <- -pc1
  }


  cellpc1 <- data.frame(Cell = pc1, stringsAsFactors = FALSE)
  row.names(cellpc1) <- names(pc1)

  return(cellpc1)

}


#Organize cell markers from processSeuratobj via
#1) remove some unusual '.' from the original gene ids
#2) add manual markers obtained from publications, if provided
#3) add genes with a high correlation to the PC1 of the above markers in a
#target data set, if provided
orgmarkers <- function(orimarkers = Seuratobjlist$cellmarkers,
                       refgenes = row.names(reftpms),
                       targetdat = targetdat,
                       manualmarkerlist = NULL,
                       markerremovecutoff = 0.6){

  markergroup <- orimarkers

  datcolsnames <- c('gene', 'celltype', 'p_val_adj', 'avg_log2FC', 'p_val', 'pct.1', 'pct.2')

  if(sum(!(datcolsnames %in% names(markergroup))) > 0){
    datcolsnames <- c('marker', 'celltype', 'P.val.adj', 'log2FC', 'P.val')
  }

  markergroup <- markergroup[datcolsnames]
  newgenename <- gsub(pattern = '\\..*$', replacement = '', x = markergroup[,1])
  markergroup$gene <- newgenename

  markergroupgenes <- unique(markergroup$gene)
  markerreservegenes <- intersect(refgenes, markergroupgenes)
  markergroup <- subset(markergroup, gene %in% markerreservegenes)


  celltypenames <- unique(markergroup$celltype)

  uniqmarkers <- function(pool = markergroup, cell){
    cellsub <- subset(pool, celltype == cell)
    othersub <- subset(pool, celltype != cell)

    cellmarkers <- unique(cellsub$gene)
    othermarkers <- unique(othersub$gene)
    cellmarkers <- setdiff(cellmarkers, othermarkers)

    return(cellmarkers)
  }

  markerlist <- list()
  for(i in 1:length(celltypenames)){

    celltypename <- celltypenames[i]
    markerlist[[i]] <- uniqmarkers(pool = markergroup, cell = celltypename)
    names(markerlist)[i] <- celltypename

    if(!is.null(manualmarkerlist)){

      celltypemanualmarkers <- manualmarkerlist[[celltypename]]
      markerlist[[i]] <- unique(c(markerlist[[i]],
                                  celltypemanualmarkers))

    }

  }

  selectedmarkerlist <- list()
  i <- 1
  for(i in 1:length(celltypenames)){
    celltypename <- celltypenames[i]

    celltypemarkers <- markerlist[[i]]

    if(!is.null(targetdat)){

      if(is.matrix(targetdat)){

        if(ncol(targetdat) > 2){

          expdata <- t(targetdat)

          cellpc1 <- genepc1(expdata = expdata,
                             celltypemarkers = celltypemarkers)

          names(cellpc1) <- celltypename


          genecor <- cor(expdata, cellpc1)
          genecor <- genecor[genecor > markerremovecutoff, , drop = FALSE]
          corgenes <- row.names(genecor)
          corgenes <- unique(corgenes)

          corgenes <- unique(c(corgenes,
                               intersect(colnames(expdata), celltypemarkers)))

          celltypemarkers <- corgenes

        }
      }
    }

    selectedmarkerlist[[i]] <- unique(celltypemarkers)
    names(selectedmarkerlist)[i] <- celltypename

  }

  return(selectedmarkerlist)


}

#Order the genes in the reference TPM table according to the difference between
#the cell type with the highest expression of a gene and with the second
#highest expression of this gene, as well as the difference between the cell
#type with the highest expression of this gene and with the third highest
#expression of this gene, and then get the top 5%, 10%, 15%, etc genes so that
#condition number of these top gene groups can be calculated later
preparecondition <- function(tpms = reftpms,
                             conditioning = TRUE,
                             manualmarkerlist = NULL){


  if(!is.null(manualmarkerlist)){

    cellmanualmarkers <- unique(as.vector(unlist(manualmarkerlist)))

  }else{
    cellmanualmarkers <- NULL
  }


  cellnames <- colnames(tpms)
  cellnames <- gsub(pattern = '_[0-9]*$', replacement = '',
                    cellnames)
  cellnames <- unique(cellnames)

  cellsubs <- list()
  refmeanslist <- list()
  i <- 1
  for(i in 1:length(cellnames)){
    cellname <- cellnames[i]
    cellname <- gsub(pattern = '\\+', replacement = '\\\\+', x = cellname)
    cellsub <- tpms[,grep(pattern = paste0('^', cellname, '_'),
                          x = colnames(tpms),
                          fixed = FALSE)]

    cellsubs[[i]] <- cellsub
    names(cellsubs)[i] <- cellname

    cellmeans <- rowMeans(cellsub)
    cellmeans <- as.matrix(cellmeans)
    colnames(cellmeans) <- cellname

    refmeanslist[[i]] <- cellmeans
    names(refmeanslist)[i] <- cellname

  }

  refmeans <- do.call(cbind, refmeanslist)
  refmeans <- as.data.frame(refmeans, stringsAsFactors = FALSE)


  i <- 1
  for(i in 1:nrow(refmeans)){
    gene <- row.names(refmeans)[i]
    geneorder <- refmeans[i,]
    geneorder <- colnames(geneorder)[order(-unlist(geneorder))]

    top1 <- geneorder[1]
    top2 <- geneorder[2]
    top3 <- geneorder[3]

    diff12 <- wilcox.test(cellsubs[[top1]][gene,], cellsubs[[top2]][gene,])
    diff13 <- wilcox.test(cellsubs[[top1]][gene,], cellsubs[[top3]][gene,])

    diff12p <- diff12$p.value
    diff13p <- diff13$p.value

    diff12 <- mean(cellsubs[[top1]][gene,]) - mean(cellsubs[[top2]][gene,])
    diff13 <- mean(cellsubs[[top1]][gene,]) - mean(cellsubs[[top3]][gene,])

    diffframe <- data.frame(gene = gene, diff12p = diff12p, diff13p = diff13p,
                            diff12 = diff12, diff13 = diff13,
                            stringsAsFactors = FALSE)

    if(i == 1){
      diffframes <- diffframe
    }else{
      diffframes <- rbind(diffframes, diffframe)
    }


  }

  diffframes$diff12p.adj <- p.adjust(diffframes$diff12p, method = 'BH')
  diffframes$diff13p.adj <- p.adjust(diffframes$diff13p, method = 'BH')

  diffframes <- diffframes[order(diffframes$diff12p.adj,
                                 diffframes$diff13p.adj,
                                 diffframes$diff12p,
                                 diffframes$diff13p,
                                 -diffframes$diff12,
                                 -diffframes$diff13),]
  row.names(diffframes) <- 1:nrow(diffframes)
  diffframes <- diffframes[c('gene', 'diff12p', 'diff13p',
                             'diff12p.adj', 'diff13p.adj',
                             'diff12', 'diff13')]

  diff12frame <- diffframes[c('gene', 'diff12p.adj', 'diff12p', 'diff12')]
  diff13frame <- diffframes[c('gene', 'diff13p.adj', 'diff13p', 'diff13')]



  diff12frame <- diff12frame[order(diff12frame$diff12p.adj, diff12frame$diff12p,
                                   -diff12frame$diff12),]
  diff13frame <- diff13frame[order(diff13frame$diff13p.adj, diff13frame$diff13p,
                                   -diff13frame$diff13),]
  row.names(diff12frame) <- 1:nrow(diff12frame)
  row.names(diff13frame) <- 1:nrow(diff13frame)

  if(!is.null(cellmanualmarkers)){

    manualmarkerdiff12frame <- subset(diff12frame, gene %in% cellmanualmarkers)
    manualmarkerdiff13frame <- subset(diff13frame, gene %in% cellmanualmarkers)

    otherdiff12frame <- subset(diff12frame, !(gene %in% cellmanualmarkers))
    otherdiff13frame <- subset(diff13frame, !(gene %in% cellmanualmarkers))

    diff12frame <- rbind(manualmarkerdiff12frame, otherdiff12frame)
    diff13frame <- rbind(manualmarkerdiff13frame, otherdiff13frame)

  }

  row.names(diff12frame) <- 1:nrow(diff12frame)
  row.names(diff13frame) <- 1:nrow(diff13frame)

  genenum <- nrow(diff12frame)

  quantities <- seq(0.05, 1, 0.05)

  if(conditioning == FALSE){
    quantities <- 1
  }

  tpmsublist <- list()
  i <- 1
  cellmanualmarkernum <- length(cellmanualmarkers)
  for(i in 1:length(quantities)){
    subgenenum <- round(genenum*quantities[i])
    subgenenum <- max(subgenenum, cellmanualmarkernum)
    diff12sub <- diff12frame[1:subgenenum,]
    diff13sub <- diff13frame[1:subgenenum,]

    diffgenes <- unique(c(diff12sub$gene, diff13sub$gene))
    tpmsub <- tpms[diffgenes,]
    tpmsublist[[i]] <- tpmsub

  }

  return(tpmsublist)


}


removetop <- function(refdat,
                      cutoff = 1,
                      kepts = NULL){

  if(cutoff == 1){
    return(refdat)
  }

  removegenenum <- ceiling(nrow(refdat)*(1 - cutoff))
  removegenes <- c()
  i <- 1
  for(i in 1:ncol(refdat)){

    cellvals <- refdat[,i]
    names(cellvals) <- row.names(refdat)
    cellvals <- cellvals[order(-cellvals)]

    removegenes <- c(removegenes, names(cellvals)[1:removegenenum])

  }

  removegenes <- unique(removegenes)

  keptgenes <- row.names(refdat)[(!row.names(refdat) %in% setdiff(removegenes, kepts))]

  keptdat <- refdat[keptgenes, , drop = FALSE]

  return(keptdat)

}




calcorfactors <- function(cell1vals, cell2vals){

  denominator <- sd(cell1vals)*sd(cell2vals)

  numerators <- (cell1vals - mean(cell1vals))*(cell2vals - mean(cell2vals))

  corfactors <- numerators/denominator

  colnames(corfactors) <- 'corfactor'

  return(corfactors)

}

reducecor <- function(ref,
                      adjustcutoff = 0.5,
                      minrefgenenum = 1000){

  cormat <- cor(ref)

  abscormat <- abs(cormat)
  abscormat[lower.tri(abscormat, diag = TRUE)] <- 0

  coords <- which(abscormat > adjustcutoff, arr.ind = TRUE, useNames = TRUE)
  coords <- coords[coords[,1] != coords[,2], , drop = FALSE]

  if(nrow(coords) == 0){
    return(ref)
  }

  row.names(coords) <- 1:nrow(coords)

  corfactorlist <- list()
  i <- 1
  for(i in 1:nrow(coords)){

    cell1 <- row.names(abscormat)[coords[i,1]]
    cell2 <- colnames(abscormat)[coords[i,2]]

    cell1vals <- ref[, colnames(ref) == cell1, drop = FALSE]
    cell2vals <- ref[, colnames(ref) == cell2, drop = FALSE]

    corfactors <- calcorfactors(cell1vals = cell1vals, cell2vals = cell2vals)
    corfactorlist[[i]] <- corfactors

  }

  corfactorlist <- do.call(cbind, corfactorlist)
  genecorfactors <- rowMeans(corfactorlist)

  sign <- cormat[abscormat > adjustcutoff]
  sign <- sum(sign)

  if(sign > 0){
    genecorfactors <- names(genecorfactors[order(-genecorfactors)])
  }else{
    genecorfactors <- names(genecorfactors[order(genecorfactors)])
  }

  cutgenenum <- length(genecorfactors) - minrefgenenum
  if(cutgenenum <= 0){
    candidategenes <- genecorfactors
  }else{
    candidategenes <- genecorfactors[1:cutgenenum]
  }

  ref <- ref[genecorfactors, , drop = FALSE]

  while(nrow(ref) > minrefgenenum & sum(abscormat[upper.tri(abscormat)] > adjustcutoff) > 0){

    ref <- ref[2:nrow(ref), , drop = FALSE]
    abscormat <- abs(cor(ref))

  }

  ref <- ref[row.names(corfactorlist)[row.names(corfactorlist) %in% row.names(ref)],]

  return(ref)

}



#If a target data set is provided, combine it with the reference, for the
#top 5%, 10%, 15% gene groups generated by preparecondition, respectively.
#And calculate the condition numbers of the top gene groups and select the
#top gene group with the smallest condition number.
combatconditon <- function(reftpmsublist,
                           targetexp,
                           targetlogged,
                           conditioning,
                           refbatch = NULL,
                           savefile = TRUE,
                           minrefgenenum,
                           cutoff = 0.99,
                           adjustcutoff = 0.5,
                           suffix = '',
                           combat = TRUE){

  pseudocount <- 1

  if(!is.null(targetexp) & targetlogged == TRUE){

    if(sum(2^targetexp - 1 < 0) > 0){

      pseudocount <- 0

    }

    targetexp <- 2^targetexp - pseudocount

    targetexp[targetexp < 0] <- 0

  }


  #library(sva)

  ref.crt.list <- list()

  if(!is.null(targetexp)){

    targetexp.list <- list()

  }


  i <- 1
  for(i in 1:length(reftpmsublist)){

    reftpmsub <- reftpmsublist[[i]]
    reftpmsub[reftpmsub < 0] <- 0

    if(!is.null(targetexp) & combat == TRUE){
      targetexpsub <- targetexp[row.names(reftpmsub), , drop = FALSE]

      #reftpmsub <- log2(reftpmsub + pseudocount)
      #reftpmsub[is.infinite(reftpmsub)] <- 0
      #targetexpsub <- log2(targetexpsub + pseudocount)
      #targetexpsub[is.infinite(targetexpsub)] <- 0


      tmp <- cbind(targetexpsub, reftpmsub)
      batchinfo <- c(rep(1, ncol(targetexpsub)), rep(2, ncol(reftpmsub)))

      tmp.corrected <- sva::ComBat(tmp, batchinfo, c(), ref.batch = refbatch)

      targetexpsub.crt <- tmp.corrected[,1:ncol(targetexpsub), drop = FALSE]
      reftpmsub.crt <- tmp.corrected[,(ncol(targetexpsub)+1):ncol(tmp.corrected),
                                     drop = FALSE]


      #reftpmsub <- 2^reftpmsub - pseudocount
      #targetexpsub <- 2^targetexpsub - pseudocount

      reftpmsub[reftpmsub < 0] <- 0
      targetexpsub[targetexpsub < 0] <- 0

      #targetexpsub.crt <- 2^targetexpsub.crt - pseudocount
      #reftpmsub.crt <- 2^reftpmsub.crt - pseudocount

      targetexpsub.crt[targetexpsub.crt < 0] <- 0
      reftpmsub.crt[reftpmsub.crt < 0] <- 0

    }else{

      reftpmsub.crt <- reftpmsub

      if(!is.null(targetexp)){
        targetexpsub.crt <- targetexp[row.names(reftpmsub), , drop = FALSE]
      }

    }


    celltypelist <- gsub(pattern = '_[0-9]*$', replacement = '', x = colnames(reftpmsub.crt))
    celltypes <- unique(celltypelist)

    j <- 1
    for(j in 1:length(celltypes)){

      celltype <- celltypes[j]
      sub <- reftpmsub.crt[,celltypelist == celltype]

      submedian <- apply(X = sub, MARGIN = 1, FUN = median)
      submedian <- data.frame(medvalue = submedian, stringsAsFactors = FALSE)
      names(submedian) <- celltype
      if(j == 1){
        refmedian <- submedian
      }else{
        refmedian <- cbind(refmedian, submedian)
      }

    }

    ref.crt <- refmedian
    #backup <- ref.crt
    #targetexpsubbackup <- targetexpsub.crt

    ref.crt.list[[i]] <- ref.crt

    if(!is.null(targetexp)){
      targetexp.list[[i]] <- targetexpsub.crt
    }


  }

  #Remove top expressed outlier genes
  keptref <- removetop(refdat = ref.crt.list[[length(ref.crt.list)]],
                       cutoff = cutoff,
                       kepts = NULL)


  if(conditioning == TRUE){

    i <- 1
    for(i in 1:length(ref.crt.list)){

      ref.crt.list[[i]] <- ref.crt.list[[i]][row.names(ref.crt.list[[i]]) %in% row.names(keptref), ]

      ref.crt.sub <- ref.crt.list[[i]]

      if(!is.matrix(ref.crt.sub)){
        ref.crt.sub <- as.matrix(ref.crt.sub)
      }

      genenum <- nrow(ref.crt.sub)
      kappaval <- kappa(ref.crt.sub)

      if(i == 1){
        genenums <- genenum
        kappavals <- kappaval
      }else{
        genenums <- c(genenums, genenum)
        kappavals <- c(kappavals, kappaval)
      }


    }


    refidx <- seq(1, length(genenums), 1)

    if(minrefgenenum > max(genenums)){
      minrefgenenum <- 0
    }

    refidx <- refidx[genenums >= minrefgenenum]
    qualifiedkappavals <- kappavals[genenums >= minrefgenenum]
    qualifiedgenenums <- genenums[genenums >= minrefgenenum]


    bestrefidx <- match(min(qualifiedkappavals), qualifiedkappavals)
    bestrefidx <- refidx[bestrefidx]

    #library(ggplot2)
    #library(scales)
    conddata <- data.frame(genenums = genenums,
                           kappavals = kappavals,
                           stringsAsFactors = FALSE)
    p <- ggplot2::ggplot(data = conddata,
                         ggplot2::aes(x = genenums, y = kappavals))

    print(
      p + ggplot2::geom_point(color = scales::hue_pal()(1), size = 2) +
        ggplot2::geom_line(color = scales::hue_pal()(1), size = 1) +
        ggplot2::xlab('Gene Number') + ggplot2::ylab('Condition Number') +
        ggplot2::ggtitle(paste0('Reference matrix gene number/condition number'),
                         subtitle = paste0('(Best gene number = ', genenums[bestrefidx],
                                           '; Condition number = ',
                                           signif(min(qualifiedkappavals), 3), ')')) +
        ggplot2::theme_bw()
    )



  }else{

    ref.crt.list[[length(ref.crt.list)]] <-
      ref.crt.list[[length(ref.crt.list)]][row.names(ref.crt.list[[length(ref.crt.list)]]) %in% row.names(keptref), ]

    bestrefidx <- length(ref.crt.list)
  }

  ref.crt <- ref.crt.list[[bestrefidx]]

  ref.crt <- as.matrix(ref.crt)

  ref.crt <- reducecor(ref = ref.crt,
                       adjustcutoff = adjustcutoff,
                       minrefgenenum = minrefgenenum)

  if(!is.null(targetexp)){

    targetexp.crt <- targetexp.list[[bestrefidx]]
    targetexp.crt <- targetexp.crt[row.names(ref.crt),]

  }

  ref.crt <- as.matrix(ref.crt)

  if(savefile == TRUE){

    if(conditioning == TRUE){
      tag <- '.cond'
    }else{
      tag <- ''
    }

    stamp <- Sys.time()
    stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
    stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
    stamp <- paste0('.', stamp)


    saveRDS(ref.crt, file = paste0('ref.final.nolog', tag, stamp, suffix, '.rds'))

    if(!is.null(targetexp)){

      saveRDS(targetexp.crt, file = paste0('targetexp.final.nolog', tag, stamp, suffix, '.rds'))

    }


  }

  if(!is.null(targetexp)){

    res <- list(ref = ref.crt,
                targetnolog = targetexpsub.crt)

  }else{
    res <- list(ref = ref.crt)
  }

  return(res)


}


#'Make cell deconvolution reference from scRNA-seq data
#'
#'Make cell deconvolution reference matrix from scRNA-seq data.
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
#'@param pseudobulknum At the beginning of making the cell reference matrix,
#'  the scRNA-seq cell counts contained in \code{Seuratobj} will be sampled
#'  and used to generate some pseudo-bulk RNA-seq samples, for each cell type.
#'  The parameter \code{pseudobulknum} here defines how many pseudo-bulk
#'  RNA-seq data for each cell type need to be generated. Default is 10.
#'@param samplebalance During generating the pseudo-bulk RNA-seq data, the
#'  number of single cells can be sampled is always different for each cell
#'  type. If want to adjust this bias and make the single cell numbers used to
#'  make pseudo-bulk RNA-seq data same for different cell types, set this
#'  parameter as TRUE. Then, the cell types with too many candidate cells will
#'  be down-sampled while the ones with much fewer cells will be over-sampled.
#'  The down-sampling is performed using bootstrapping, and the over-sampling
#'  is conducted with SMOTE (Synthetic Minority Over-sampling Technique). This
#'  is a time-consuming step and the default value of this parameter is FALSE.
#'@param pseudobulkpercent If the parameter \code{samplebalance} is FALSE, for
#'  the pseudo-bulk sampling for each cell type, a percent of single cells for
#'  each cell type will be randomly sampled and this parameter is used to set
#'  this percent value and should be a number between 0 and 1, but if the
#'  parameter \code{samplebalance} is set as TRUE, bootstrapping and SMOTE
#'  will be performed to do the sampling and this parameter will be omitted.
#'@param pseudobulkdat If the scRNA-seq data transferred via \code{Seuratobj}
#'  is large, the pseudo-bulk RNA-seq data generation step will become time-
#'  consuming, and if this same scRNA-seq data needs to be used repeatedly for
#'  deconvolving different bulk datasets, to save time, it is recommended to
#'  use the function \code{prepseudobulk} to generate and save the pseudo-bulk
#'  RNA-seq data at the first time, and then the data can be transferred to
#'  this parameter \code{pseudobulkdat}, so that \code{scRef} can always skip
#'  its own pseudo-bulk data generation step and directly use the data here to
#'  further generate the final RNA deconvolution reference. The default value
#'  of this parameter is NULL, and in this case, the synthesis step will not
#'  be skipped and \code{scRef} will synthesize the pseudo-bulk data itself.
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
#'  gene IDs and column names are sample IDs. The default value of it is NULL.
#'  In this case, the reference matrix generation step will only base on the
#'  scRNA-seq data provided by the previous parameter \code{Seuratobj}, but if
#'  provide a matrix to be deconvolved to this parameter, both the reference
#'  matrix and this cell mixture matrix will be further processed, including
#'  combining the 2 matrices to remove their batch difference, selecting more
#'  genes into the reference matrix based on the correlation between the genes
#'  and the selected marker genes in the cell mixture , etc. It is recommended
#'  to provide the cell mixture matrix via this parameter, especially when the
#'  cell mixture is from RNA microarray, rather than RNA-seq data, so that the
#'  combination process will be performed to reduce the platform difference
#'  between RNA microarray and scRNA-seq.
#'@param targetlogged Whether the gene expression values in \code{targetdat}
#'  are log2 transformed values or not.
#'@param manualmarkerlist During making the reference matrix, for each cell
#'  type, the genes specially expressed in it with a high level will be deemed
#'  as markers and further used to generate the reference. However, it cannot
#'  be ensured that some known classical markers can be selected, and so if
#'  want to make sure these markers can be used to make the reference, a list
#'  can be used as an input to this parameter, with its element names as the
#'  cell type names and the elements as vectors with the gene IDs of these
#'  classical markers. It should be noted that before the final reference is
#'  determined, all the marker genes need to go through several filter steps,
#'  such as extremely highly expressed genes and collinearity contributing
#'  genes removal, to improve the reference quality, so that the classical
#'  genes provided via this parameter will be definitely used for reference
#'  generation, but may also be filtered out before the final one is returned.
#'  The default value of this parameter is NULL.
#'@param markerremovecutoff When a gene expression matrix is provided to the
#'  parameter \code{targetdat}, the gene expression values in it will be used
#'  to calculate the correlation with the scRNA-seq selected markers in this
#'  cell mixture matrix and the ones with a high Pearson correlation to the
#'  first principle component of these marker genes will also be used to make
#'  the reference. The cutoff of the Pearson correlation coefficient is set by
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
#'@param savefile Whether need to save the finally generated reference matrix,
#'  and the adjusted cell mixture matrix (if provided to \code{targetdat}),
#'  as rds file(s) in the working directory automatically. Default is FALSE.
#'@param threads Number of threads need to be used to do the computation. Its
#'  default value is 1.
#'@param cutoff To improve the robustness of the deconvolution result, some
#'  extremely highly expressed genes in the reference need to be filtered out
#'  due to their large variance. This cutoff is used to set the percent of
#'  genes can be kept in the reference while the other genes with a higher
#'  expression level will be filtered. The default value is 0.95, meaning the
#'  top 5% most highly expressed genes will be removed from the reference.
#'@param adjustcutoff For some similar cell types, their gene expressions in
#'  the reference matrix have a large correlation, which makes the downstream
#'  deconvolution difficult. To relive this problem, for each similar cell
#'  pair, some genes largely contributing to their correlation will be found
#'  and removed, so that their correlation in the reference can be reduced.
#'  This parameter \code{adjustcutoff} is used to set the cutoff of the cell
#'  pair correlation, and if a cell pair has a Pearson correlation coefficient
#'  greater than this value, the contributing gene filter process will be used
#'  to reduce the coefficient until it becomes smaller than this value. The
#'  default value is 0.4.
#'@return A list with the final reference matrix as its element, and if the
#'  cell mixture data matrix to be deconvolved is provided to the parameter
#'  \code{targetdat}, a adjusted one will also be returned as an element of
#'  this list. The gene values in this adjusted matrix are non-log transformed
#'  values.
#'@export
scRef <- function(Seuratobj,
                  targetcelltypes = NULL,
                  celltypecolname = "annotation",
                  pseudobulknum = 10,
                  samplebalance = FALSE,
                  pseudobulkpercent = 0.9,
                  pseudobulkdat = NULL,
                  geneversion = "hg19",
                  genekey = "SYMBOL",
                  targetdat = NULL,
                  targetlogged = FALSE,
                  manualmarkerlist = NULL,
                  markerremovecutoff = 0.6,
                  minrefgenenum = 500,
                  savefile = FALSE,
                  threads = 1,
                  cutoff = 0.95,
                  adjustcutoff = 0.4){

  TPMorFPKM <- 'TPM'

  Seuratobjlist <- processSeuratobj(Seuratobj = Seuratobj,
                                    targetcells = targetcelltypes,
                                    celltypecol = celltypecolname)

  if(!is.null(targetcelltypes)){

    #Extract and organize cell meta data
    meta <- Seuratobj[[]]
    meta <- meta[c('orig.ident', 'nCount_RNA', 'nFeature_RNA',
                   celltypecolname)]
    names(meta)[ncol(meta)] <- 'celltype'

    targetcelltypenum <- length(unique(targetcelltypes[targetcelltypes %in% meta$celltype]))

    if(length(unique(Seuratobjlist$cellmarkers$celltype)) < ceiling(targetcelltypenum*0.5)){

      Seuratobjlist <- processSeuratobj(Seuratobj = Seuratobj,
                                        targetcells = NULL,
                                        celltypecol = celltypecolname)


    }

  }

  #Seuratobjlist$meta$celltype[Seuratobjlist$meta$celltype == 'fFB1'] <- 'F'
  #Seuratobjlist$cellmarkers$celltype[Seuratobjlist$cellmarkers$celltype == 'fFB1'] <- 'F'

  if(is.null(pseudobulkdat)){

    refcounts <- makeref.cv(scmatrix = Seuratobjlist$counts,
                            metainfo = Seuratobjlist$meta,
                            nround = pseudobulknum,
                            samplepercent = pseudobulkpercent,
                            balance = samplebalance,
                            threads = threads)


  }else{

    refcounts <- pseudobulkdat

  }


  tag <- ''

  if(savefile == TRUE){

    if(samplebalance == TRUE){
      tag <- '.balance'
    }else{
      tag <- '.nobalance'
    }

    stamp <- Sys.time()
    stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
    stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
    stamp <- paste0('.', stamp)

    if(is.null(pseudobulkdat)){
      saveRDS(refcounts, file = paste0('orirefcounts', stamp, tag, '.rds'))
    }

  }

  #LUAD & LUSC
  #refcounts <- readRDS('orirefcounts.2021-09-19_21-01-45.nobalance.rds')
  #PBMC
  #refcounts <- readRDS('orirefcounts.2021-11-16_20-57-31.nobalance.rds')
  #PBMC filterred
  #refcounts <- readRDS('orirefcounts.2022-02-05_17-48-01.nobalance.rds')
  #PBMC ATAC
  #refcounts <- readRDS('orirefcounts.2021-11-16_20-57-31.nobalance.rds')

  refcounts <- refcounts[rowSums(refcounts) != 0, , drop = FALSE]

  reftpms <- caltpms(countmat = refcounts,
                     version = geneversion,
                     genekey = genekey,
                     TPMorFPKM = TPMorFPKM)

  reflogtpms <- log2(reftpms + 1)

  if(!is.null(targetdat)){

    if(targetlogged == FALSE){
      targetdat <- log2(targetdat + 1)
      targetlogged <- TRUE
    }

    targetdat <- orggeneids(oriref = targetdat)

    sharedgenes <- intersect(row.names(targetdat),
                             row.names(reftpms))

    #targetdatuniqgenes <- setdiff(row.names(targetdat),
    #                              row.names(reftpms))
    #refuniqgenes <- setdiff(row.names(reftpms),
    #                        row.names(targetdat))

    targetdat <- targetdat[sharedgenes, , drop = FALSE]
    reftpms <- reftpms[sharedgenes, , drop = FALSE]

    targetdat <- removenonsds(dat = targetdat)
    reftpms <- removenonsds(dat = reftpms)

    sharedgenes <- intersect(row.names(targetdat),
                             row.names(reftpms))
    targetdat <- targetdat[sharedgenes, , drop = FALSE]
    reftpms <- reftpms[sharedgenes, , drop = FALSE]


    reflogtpms <- reflogtpms[sharedgenes,]

    Seuratobjlist$cellmarkers <- subset(Seuratobjlist$cellmarkers,
                                        gene %in% sharedgenes)


  }

  cellmarkerlist <- orgmarkers(orimarkers = Seuratobjlist$cellmarkers,
                               refgenes = row.names(reftpms),
                               targetdat = targetdat,
                               manualmarkerlist = manualmarkerlist,
                               markerremovecutoff = markerremovecutoff)

  markers <- do.call(c, cellmarkerlist)
  markers <- unique(markers)

  if(!is.null(targetdat)){

    markers <- intersect(markers, row.names(targetdat))
    targetdat <- targetdat[markers, , drop = FALSE]

  }

  reftpms <- reftpms[markers, ,drop = FALSE]
  reflogtpms <- reflogtpms[markers, ,drop = FALSE]

  #manualmarkerlist <- NULL

  reftpmsublist <- preparecondition(tpms = reftpms,
                                    conditioning = FALSE,
                                    manualmarkerlist = manualmarkerlist)

  res <- combatconditon(reftpmsublist = reftpmsublist,
                        targetexp = targetdat,
                        targetlogged = targetlogged,
                        conditioning = FALSE,
                        savefile = FALSE,
                        minrefgenenum = minrefgenenum,
                        cutoff = cutoff,
                        adjustcutoff = adjustcutoff,
                        suffix = tag)

  if(!is.null(targetcelltypes)){
    sharedcelltypes <- intersect(colnames(res$ref), targetcelltypes)
    res$ref <- res$ref[,colnames(res$ref) %in% sharedcelltypes, drop = FALSE]
  }


  if(savefile == TRUE){

    stamp <- Sys.time()
    stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
    stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
    stamp <- paste0('.', stamp)


    saveRDS(res$ref, file = paste0('ref.final.nolog', stamp, tag, '.rds'))

    if(!is.null(targetdat)){

      saveRDS(res$targetnolog, file = paste0('targetexp.final.nolog', stamp, tag, '.rds'))

    }


  }


  return(res)

}


