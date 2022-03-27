
#epDeconv####

#Find methylation probes highly correlated to the cell contents
highrelatedprobes <- function(methylmat,
                              cellcontents,
                              corprobecutoff,
                              removedups = FALSE){


  #Select probes highly correlated to cell contents#######
  corprobes <- list()

  k <- 1
  for(k in 1:ncol(cellcontents)){
    cellname <- colnames(cellcontents)[k]
    cellvalue <- cellcontents[,k]
    probecor <- cor(methylmat, cellvalue)

    colnames(probecor) <- cellname
    colprobes <- probecor[abs(probecor) > corprobecutoff,]
    probenames <- names(colprobes)

    print(length(probenames))

    if(length(probenames) != 0){
      corprobes[[cellname]] <- probenames
    }

  }

  corprobes <- unlist(corprobes)


  if(removedups == TRUE){

    dups <- unique(corprobes[duplicated(corprobes)])

    corprobes <- corprobes[!(corprobes %in% dups)]
  }



  corprobes <- unique(corprobes)

  if(length(corprobes) > 10){
    corbetas <- methylmat[,corprobes]
  }else{
    corbetas <- methylmat
  }

  corbetaslist <- list()

  for(l in 1:ncol(corbetas)){
    corprobename <- colnames(corbetas)[l]
    corbetaslist[[l]] <- corbetas[,l]
    names(corbetaslist)[l] <- corprobename
  }

  return(corbetaslist)

}

#Calculate the methylation cell reference for one probe from a single probe
#methylation data across multiple samples (`expmat`) and cell contents of
#the samples (`refmat`)
methylconstraintreg <- function(expmat,
                                refmat){

  tst <- lsfit(x = refmat, y = expmat, intercept = FALSE)
  coeffs <- tst$coefficients
  return(coeffs)

}

#Parallelize the methylation cell reference calculation for multiple probes
#from multiple probe methylation data list (`corbetaslist`) and cell contents
#of the samples (`bootsolutions`)
makemethylref <- function(corbetaslist,
                          bootsolutions,
                          threads = 1){

  methylconstraintreg <- function(expmat,
                                  refmat){

    tst <- lsfit(x = refmat, y = expmat, intercept = FALSE)
    coeffs <- tst$coefficients
    return(coeffs)

  }

  if(threads == 1){

    methylrefs <- list()
    for(o in 1:length(corbetaslist)){

      corbetas <- corbetaslist[[o]]
      methylref <- methylconstraintreg(expmat = corbetas,
                                       refmat = bootsolutions)

      methylrefs[[o]] <- methylref
    }

  }else{

    #library(doParallel)

    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))

    doParallel::registerDoParallel(cl)

    #date()
    `%dopar%` <- foreach::`%dopar%`
    methylrefs <- foreach::foreach(expmat = corbetaslist,
                                   #.export = ls(name = globalenv())) %dopar% {
                                   .export = NULL) %dopar% {
                                     methylconstraintreg(expmat,
                                                         refmat = bootsolutions)
                                   }
    #.export character vector of variables to export. This can be useful when
    #accessing a variable that isn't defined in the current environment. The
    #default value is NULL
    #The function methylconstraintreg is used by foreach, if methylconstraintreg
    #is defined within this function makemethylref, and .export = NULL, when run
    #makemethylref, there will be no error. If methylconstraintreg is defined
    #out of this function makemethylref, .export should be set as
    #ls(name = globalenv()), otherwise an error will occur when run the whole
    #function makemethylref
    #date()

    parallel::stopCluster(cl)

    unregister_dopar()

  }

  names(methylrefs) <- names(corbetaslist)
  methylrefsdat <- do.call(rbind, methylrefs)

  return(methylrefsdat)


}

methylpres <- function(methylrefs,
                       elasticnetmodel,
                       targetmethyldat,
                       targetmethyldatlogged = FALSE,
                       #deconvmod,
                       resscale = FALSE,
                       adjustminus = FALSE,
                       baselearner = 'lasso'){

  if(baselearner == 'lasso'){

    targetmethyldat <- targetmethyldat[, row.names(elasticnetmodel$beta[[1]]),
                                       drop = FALSE]

    pres <- tryCatch({

      predict(object = elasticnetmodel, newx = targetmethyldat)

    }, error = function(err){

      NULL

    })

    if(is.null(pres)){

      library(glmnet)

      pres <- predict(object = elasticnetmodel, newx = targetmethyldat)

    }


  }else{

    #Predict methyl cell contents#####
    targetmethyldat <- t(targetmethyldat)

    #targetdat should be a matrix with genes/probes as rows and samples as columns
    pres <- refDeconv(ref = methylrefs,
                      targetdat = targetmethyldat,
                      targetlogged = targetmethyldatlogged,
                      #modmethod = deconvmod,
                      resscale = resscale)

  }

  pres <- as.data.frame(pres)
  names(pres) <- gsub(pattern = '\\.s0', replacement = '', x = names(pres))

  if(adjustminus == TRUE){
    pres[pres < 0] <- 0
  }

  if(resscale == TRUE){

    pres <- pres/rowSums(pres)

  }

  return(pres)

}

calbaseweight <- function(rnares,
                          methylres){

  uniquesamples <- row.names(methylres)
  changednames <- uniquesamples[!(uniquesamples %in% row.names(rnares))]
  orinames <- gsub(pattern = '\\.[0-9]*$', replacement = '', x = changednames)

  uniquesamples[!(uniquesamples %in% row.names(rnares))] <- orinames

  rnares <- rnares[uniquesamples,]

  rnacellcor <- cor(rnares)
  methylcellcor <- cor(methylres)

  corcor <- cor(rnacellcor, methylcellcor)

  corcorrsqr <- diag(corcor)

  if(sum(is.na(corcorrsqr)) > 0){
    return(NULL)
  }

  if(any(corcorrsqr <= sqrt(0.5))){
    return(NULL)
  }

  corcorrsqr <- corcorrsqr^2
  weight <- 0.5*log(corcorrsqr/(1 - corcorrsqr))
  return(weight)

}


#'Predict cell contents for new methylation data using an ensemble model
#'
#'Predict cell contents for new methylation data with the model trained from
#'the function \code{epDeconv}.
#'
#'@param model Result returned from the function \code{epDeconv}. It contains
#'  the deconvoltion ensemble model. Default value is NULL, and in this case,
#'  the other two parameters \code{normweights} and \code{modellist} will be
#'  used instead to finish the prediction, if they are not NULL.
#'@param normweight The slot "normweigth" of the result of \code{epDeconv}. It
#'  provides the base learner weights of the deconvolution ensemble model. If
#'  \code{model} is NULL, this parameter will be used to restore the model and
#'  finish the prediction. Default value is NULL.
#'@param modellist The slot "modellist" of the result of \code{epDeconv}, and
#'  it provides the base learners of the deconvolution ensemble model. If the
#'  parameter \code{model} is NULL, it will be used to restore the model, and
#'  finish the prediction. Default value is NULL.
#'@param targetmethyldat The target cell mixture methylation data need to be
#'  deconvolved. Should be a matrix with each column representing one sample
#'  and each row for one feature. Row names are feature names and column names
#'  are sample IDs. It is recommended to adjust the batch difference between
#'  this dataset and the paired methylation data used by \code{epDeconv} with
#'  \code{ComBat} in advance, and using the paired data as the reference batch
#'  when adjusting, so that the cell deconvolution model trained can be used
#'  on these target data with the influence from batch difference minimized.
#'@param resscale For each sample, whether its cell contents result should be
#'  scaled so that the sum of different cell types is 1. Default is FALSE.
#'@param adjustminus In some extreme situations, the cell contents predicted
#'  may be minus values and this parameter indicates whether these minus ones
#'  should be changed to 0 in the final result. Default is FALSE.
#'@return A matrix recording the cell composition result for the samples.
#'@export
methylpredict <- function(model = NULL,
                          normweights = NULL,
                          modellist = NULL,
                          targetmethyldat,
                          resscale = FALSE,
                          adjustminus = FALSE){

  methylreflist <- NULL

  if(!is.null(model)){

    if('normweights' %in% names(model)){
      normweights <- model$normweights
    }else if(!is.null(normweights)){
      normweights <- normweights
    }else{
      return(NULL)
    }

    if('methylreflist' %in% names(model)){

      methylreflist <- model$methylreflist
      baselearner <- 'methylref'

    }else if('modellist' %in% names(model)){

      modellist <- model$modellist
      baselearner <- 'lasso'

    }else if(!is.null(methylreflist)){

      methylreflist <- methylreflist
      baselearner <- 'methylref'


    }else if(!is.null(modellist)){

      modellist <- modellist
      baselearner <- 'lasso'

    }else{
      return(NULL)
    }

  }else if(!is.null(normweights) & !is.null(methylreflist)){

    normweights <- normweights
    methylreflist <- methylreflist
    baselearner <- 'methylref'

  }else if(!is.null(normweights) & !is.null(modellist)){

    normweights <- normweights
    modellist <- modellist
    baselearner <- 'lasso'

  }else{
    return(NULL)
  }

  methylcellcontlist <- list()

  if(baselearner == 'lasso'){

    for(i in 1:length(modellist)){

      singlemodel <- modellist[[i]]

      methylcellcont <- methylpres(elasticnetmodel = singlemodel,
                                   targetmethyldat = t(targetmethyldat),
                                   baselearner = 'lasso',
                                   resscale = FALSE,
                                   adjustminus = adjustminus)

      methylcellcont <- as.matrix(methylcellcont)

      methylcellcontlist[[i]] <- methylcellcont

    }


  }else{

    for(i in 1:length(methylreflist)){

      singlemethylref <- methylreflist[[i]]

      methylcellcont <- methylpres(methylrefs = singlemethylref,
                                   targetmethyldat = t(targetmethyldat),
                                   #deconvmod = deconvmod,
                                   resscale = resscale,
                                   adjustminus = adjustminus)

      methylcellcont <- as.matrix(methylcellcont)

      methylcellcontlist[[i]] <- methylcellcont


    }

  }

  methylcellconts <- ensemblecellcontents(cellcontlist = methylcellcontlist,
                                          normweights = normweights)

  if(resscale == TRUE){

    methylcellconts <- methylcellconts/rowSums(methylcellconts)

  }

  return(methylcellconts)

}


protectnames <- function(changednames,
                         orinames){

  orinchars <- nchar(orinames)
  changednchars <- nchar(changednames)

  suffix <- changednames[orinchars != changednchars]
  suffix <- gsub(pattern = '^.*\\.', replacement = '', x = suffix)

  newnames <- orinames
  newnames[orinchars != changednchars] <- paste0(newnames[orinchars != changednchars],
                                                 '.', suffix)

  return(newnames)

}


#'Deconvolve bulk DNA methylation data using RNA reference
#'
#'Deconvolve bulk DNA methylation data with RNA reference and also by a paired
#'bulk RNA-bulk DNA methylation dataset
#'
#'@param rnaref The RNA reference recording the signature of each cell type.
#'  Each row is one gene, and each column is one cell type. Each entry should
#'  be a gene TPM value. Column names are cell type names and row names are
#'  gene names. The default is NULL and in this case, it can be synthesized
#'  from the scRNA-seq data transferred to the parameter \code{Seuratobj}.
#'@param Seuratobj An object of class Seurat generated with the \code{Seurat}
#'  R package from scRNA-seq data, should contain read count data, normalized
#'  data, and cell meta data. The meta data should contain a column recording
#'  the cell type name of each cell. When \code{rnaref} is set as NULL, but
#'  this parameter is provided with the matching data, it will be used to make
#'  the RNA reference for the downstream deconvolution.
#'@param targetcelltypes When use \code{Seuratobj} to make the RNA reference,
#'  this parameter defines the cell types should be coverred by the reference.
#'  If NULL, all the cell types included in \code{Seuratobj} will be included.
#'  Default is NULL.
#'@param celltypecolname When use \code{Seuratobj} to make the RNA reference,
#'  this parameter indicates which column in its "meta.data" slot records the
#'  cell type information for each cell and the name of this column should be
#'  transferred to this parameter. Default value is "annotation".
#'@param samplebalance When use \code{Seuratobj} to make the RNA reference, at
#'  the beginning, the scRNA-seq cell counts data in \code{Seuratobj} will be
#'  sampled and used to make 100 pseudo-bulk RNA-seq samples, for each cell
#'  type, and during synthesizing, the number of single cells can be sampled
#'  is always different for each cell type. If want to adjust the bias and
#'  make the single cell numbers used to generate pseudo-bulk RNA-seq data
#'  same for different cell types, set this parameter as TRUE. Then, the cell
#'  types with too many candidate cells will be down-sampled while the ones
#'  with much fewer cells will be over-sampled. The down-sampling is performed
#'  with bootstrapping, and the over-sampling is via SMOTE (Synthetic Minority
#'  Over-sampling Technique). This is a time-consuming step and the default is
#'  FALSE, so that no such adjustment will be performed during generating the
#'  pseudo-bulk samples.
#'@param pseudobulkdat If the scRNA-seq data transferred via \code{Seuratobj}
#'  is large, the pseudo-bulk RNA-seq data generation step will become time-
#'  consuming, and if this same scRNA-seq data needs to be used repeatedly for
#'  deconvolving different bulk datasets, to save time, it is recommended to
#'  use the function \code{prepseudobulk} to generate and save the pseudo-bulk
#'  RNA-seq data at the first time, and then the data can be transferred to
#'  this parameter \code{pseudobulkdat}, so that \code{epDeconv} can skip its
#'  own pseudo-bulk data generation step and use the data here to generate the
#'  final RNA deconvolution reference. The default value of this parameter is
#'  NULL, and in this case, the synthesis step will not be skipped.
#'@param geneversion To calculate the TPM value of the genes when generating
#'  the reference matrix, the effective length of the genes will be needed.
#'  This parameter is used to define from which genome version the effective
#'  gene length will be extracted. For human genes, "hg19" or "hg38" can be
#'  used, for mouse, "mm10" can be used. Default is "hg19".
#'@param genekey The type of the gene IDs used in the \code{Seuratobj}, it is
#'  "SYMBOL" in most cases, and the default value of this parameter is also
#'  "SYMBOL", but sometimes it may be "ENTREZID", "ENSEMBL", or other types.
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
#'@param rnamat The RNA data of the paired bulk RNA-bulk methylation dataset.
#'  Its sample cell contents will be first deconvolved via the RNA reference
#'  provided to the parameter \code{rnaref}, or generated by \code{Seuratobj},
#'  then downstream steps will be started to fulfill the bulk methylation data
#'  deconvolution. Should be a matrix with each column represennting a sample
#'  and each row for one gene. Row names are gene names and column names are
#'  sample IDs. If the reference matrix is transferred via \code{rnaref} and
#'  generated with the function \code{scRef}, and both the scRNA-seq and this
#'  paired RNA dataset were transferred to it, the result reference matrix can
#'  be transferred to \code{rnaref} and the adjusted paired RNA data returned
#'  by \code{scRef} can be transferred to this parameter.
#'@param methylmat The DNA methylaiton data of the paired bulk RNA-bulk DNA
#'  methylaiton dataset. Should be a matrix with each column represennting a
#'  sample and each row representing a feature. Row names are feature names
#'  and column names are sample IDs. The sample IDs should be the same as the
#'  ones in \code{rnamat}, because they are data for paired samples.
#'@param leanernum The base leaner number for the bagging model to deconvolve
#'  DNA methylation data. Default is 10.
#'@param rnamatlogged A logical value indicating whether the gene values in
#'  \code{rnamat} are log2 transformed or not.
#'@param resscale For each sample, whether its cell contents result should be
#'  scaled so that the sum of different cell types is 1. Default is FALSE.
#'@param threads Number of threads need to be used to do the computation. Its
#'  default value is 1.
#'@param lassoerrortype The base learners of the bagging model to deconvolve
#'  the DNA methylation data are LASSO models and the lambda value for each
#'  of them (regularization coefficient) is selected from a grid search. This
#'  parameter is used to determine whether the lambda value should be the one
#'  giving the minimum cross-validation error (set it as "min"), or the one
#'  giving an error within 1 standard error of the minimum (set it as "1se").
#'  Default is "min".
#'@param targetmethyldat The target cell mixture methylation data need to be
#'  deconvolved. Should be a matrix with each column representing one sample
#'  and each row for one feature. Row names are feature names and column names
#'  are sample IDs. It is recommended to adjust the batch difference between
#'  this dataset and \code{methylmat} with \code{ComBat} in advance, and using
#'  \code{methylmat} as the reference batch when adjusting, so that the cell
#'  deconvolution model trained from \code{methylmat} can be transferred to
#'  these data with the influence from batch difference minimized. The default
#'  value of this parameter is NULL, and it won't influence the deconvolution
#'  model training, and the model returned by this function can still be used
#'  on other cell mixture data via the function \code{methylpredict}.
#'@param plot Whether generate box plots, heatmaps, and scatter plots for the
#'  deconvolution results for the paired RNA data, paired methylaiton data,
#'  and target methylation data. Default is FALSE.
#'@param pddat If set \code{plot} as TRUE, this parameter can be used to show
#'  the sample group information of the paired bulk RNA-bulk DNA methylation
#'  data, so that their box plots will also compare the group difference for
#'  each cell type, and heatmaps with this comparison will also be generated.
#'  It should be a data frame recording the sample groups, and must include 2
#'  columns. One is named as "sampleid", recording the sample IDs same as the
#'  column names of \code{rnamat} and \code{methylmat}, the other column is
#'  "Samplegroup", recording the sample group to which each sample belongs. It
#'  can also be NULL, meaning all the samples are from the same group.
#'@param targetmethylpddat If \code{plot} is TRUE, and \code{targetmethyldat}
#'  is also provided, this parameter can be used to indicate the sample group
#'  information of the target DNA methylaiton data, its format requirment and
#'  effect are similar to \code{pddat} on the paired dataset.
#'@return A list containing several slots recording the deconvolution results
#'  for the paired RNA and paired DNA methylation data (slots "rnacellconts"
#'  and "methylcellconts"), the base leaners of the cell deconvolution model
#'  (slots "modellist" and "modelcoeflist"), the weights of the base learners
#'  (slots "normweights" and "weights"), the gene subsets used by each RNA
#'  data deconvolution base learner (slot "rnageneidxlist"), and the paired
#'  RNA-methylation sample cell contents correlation (expressed as R square)
#'  deconvolved by each base learner (slot "rnamethylsqrs"). If the target DNA
#'  methylation data is provided to the parameter \code{targetmethyldat}, a
#'  slot recording its cell contents result predicted by the model will also
#'  be returned (slot "methyltargetcellcounts").
#'@export
epDeconv <- function(rnaref = NULL,

                     Seuratobj = NULL,
                     targetcelltypes = NULL,
                     celltypecolname = 'annotation',
                     samplebalance = FALSE,
                     pseudobulkdat = NULL,
                     geneversion = 'hg19',
                     genekey = "SYMBOL",
                     manualmarkerlist = NULL,

                     rnamat,
                     methylmat,
                     learnernum = 10,
                     rnamatlogged,
                     resscale = FALSE,
                     threads = 1,
                     lassoerrortype = 'min',

                     targetmethyldat = NULL,
                     plot = FALSE,
                     pddat = NULL,
                     targetmethylpddat = NULL){

  corprobecutoff <- 0.5
  weightfrom <- 'boot'
  lassoprescreen = FALSE


  if(is.null(rnaref) & !is.null(Seuratobj)){

    refres <- scRef(Seuratobj = Seuratobj,
                    targetcelltypes = targetcelltypes,
                    celltypecolname = celltypecolname,
                    pseudobulknum = 100,
                    samplebalance = samplebalance,
                    pseudobulkpercent = 0.9,
                    pseudobulkdat = pseudobulkdat,
                    geneversion = geneversion,
                    genekey = genekey,
                    targetdat = rnamat,
                    targetlogged = rnamatlogged,
                    manualmarkerlist = manualmarkerlist,
                    markerremovecutoff = 0.6,
                    #conditioning = conditioning,
                    minrefgenenum = 1000,
                    savefile = TRUE,
                    threads = threads)

    rnaref <- refres$ref
    rnamat <- refres$targetnolog
    rnamatlogged <- FALSE

  }else if(is.null(rnaref) & is.null(Seuratobj)){

    cat('Either a RNA reference matrix or a scRNA-seq Seurat object should be provided/n')
    return(NULL)
  }

  rnasharedgenes <- intersect(row.names(rnaref),
                              row.names(rnamat))

  rnaref <- rnaref[rnasharedgenes,]
  rnamat <- rnamat[rnasharedgenes,]

  probenum <- nrow(methylmat)
  subprobenum <- round(probenum/3)

  modellist <- list()
  modelcoeflist <- list()


  weightlist <- list()
  weightrsqrlist <- list()
  cellcontlist <- list()
  rnageneidxlist <- list()

  n <- 1
  i <- 1

  while(n <= learnernum){

    i <- i + 1

    set.seed(i)

    bootgenes <- sample(x = 1:nrow(rnamat), size = nrow(rnamat), replace = TRUE)
    bootrnamat <- rnamat[bootgenes,]
    bootref <- rnaref[bootgenes,]

    #When deconvolve RNA data, only genes are bootstrapped, samples are still
    #the same as original ones

    constraintsolutions <- refDeconv(ref = bootref,
                                     targetdat = bootrnamat,
                                     targetlogged = rnamatlogged,
                                     resscale = resscale,

                                     plot = FALSE,
                                     pddat = NULL,

                                     threads = threads)

    row.names(constraintsolutions) <- colnames(bootrnamat)



    methylmat <- methylmat[,row.names(constraintsolutions)]

    set.seed(i)

    #When deconvolve methylation data, samples are bootstrapped, and probes are
    #randomly selected
    bootsamples <- sample(x = 1:ncol(methylmat), size = ncol(methylmat), replace = TRUE)
    bootmethylmat <- methylmat[,bootsamples]

    set.seed(i)
    bootprobes <- sample(x = 1:nrow(bootmethylmat), size = subprobenum, replace = FALSE)
    bootmethylmat <- bootmethylmat[bootprobes,]

    bootmethylmat <- t(bootmethylmat)
    bootsolutions <- constraintsolutions[rownames(bootmethylmat),]

    #Add oob sample data####
    oobsamples <- setdiff(colnames(methylmat), unique(colnames(methylmat)[bootsamples]))
    oobmethylmat <- methylmat[,oobsamples]
    oobmethylmat <- oobmethylmat[bootprobes,]
    oobmethylmat <- t(oobmethylmat)
    oobsolutions <- constraintsolutions[rownames(oobmethylmat),]

    #Add complete sample data####
    completemethylmat <- t(methylmat)
    completemethylmat <- completemethylmat[,colnames(bootmethylmat)]
    completesolutions <- constraintsolutions[rownames(completemethylmat),]




    if(lassoprescreen == TRUE){

      corbetaslist <- highrelatedprobes(methylmat = bootmethylmat,
                                        cellcontents = bootsolutions,
                                        corprobecutoff = 0.5,
                                        removedups = TRUE)
      corbetasprobes <- unique(names(corbetaslist))
      bootmethylmat <- bootmethylmat[,corbetasprobes]

      rm(corbetaslist)
      rm(corbetasprobes)

    }

    lassoreslist <- singlelasso(bootmethylmat = bootmethylmat,
                                bootsolutions = bootsolutions,
                                threads = threads,
                                seednum = i,
                                errortype = lassoerrortype)



    elasticnetmodel <- lassoreslist$elasticnetmodel
    modelcoef <- lassoreslist$modelcoef


    if(weightfrom == 'oob'){

      pres <- methylpres(elasticnetmodel = elasticnetmodel,
                         targetmethyldat = oobmethylmat,
                         baselearner = 'lasso',
                         adjustminus = TRUE)

      rnapres <- oobsolutions

    }else if(weightfrom == 'complete'){

      pres <- methylpres(elasticnetmodel = elasticnetmodel,
                         targetmethyldat = completemethylmat,
                         baselearner = 'lasso',
                         adjustminus = TRUE)

      rnapres <- completesolutions

    }else{

      pres <- methylpres(elasticnetmodel = elasticnetmodel,
                         targetmethyldat = bootmethylmat,
                         baselearner = 'lasso',
                         adjustminus = TRUE)

      rnapres <- bootsolutions

    }

    row.names(pres) <- protectnames(changednames = row.names(pres),
                                    orinames = row.names(rnapres))


    modelcoeflist[[n]] <- modelcoef
    modellist[[n]] <- elasticnetmodel



    pretruecormat <- cor(pres, rnapres)
    pretruersqr <- diag(pretruecormat)

    weightrsqr <- pretruersqr
    weightpres <- pres

    #Select base learners with a high correlation between its predictions and
    #the true values (RNA predicted cell contents), and only the base learners
    #predicting ALL the cell type contents well would be kept, so the base
    #learners are selected according to their fitness to the RNA predicted
    #cell contents
    if(sum(is.na(weightrsqr)) > 0){
      next()
    }

    if(any(weightrsqr <= sqrt(0.5))){
      next()
    }

    #Calculate base learner weights according to the CORRELATION BETWEEN
    #RNA predicted cell type content CORRELATION MATRIX AND methylation
    #predicted cell type content CORRELATION MATRIX, and for each base
    #learner, it will get several base learner weights for all the cell
    #types deconvolved
    weight <- calbaseweight(rnares = constraintsolutions,
                            methylres = weightpres)

    if(is.null(weight)){
      next()
    }


    weightlist[[n]] <- weight

    weightrsqrlist[[n]] <- weightrsqr

    cellcontlist[[n]] <- constraintsolutions

    rnageneidxlist[[n]] <- bootgenes

    n <- n + 1

  }

  #Correlation between methylaiton predicted values and RNA true values used
  #to select base learners
  rnamethylsqrs <- do.call(rbind, weightrsqrlist)
  row.names(rnamethylsqrs) <- 1:nrow(rnamethylsqrs)

  #Base learner weights
  weights <- do.call(rbind, weightlist)
  row.names(weights) <- 1:nrow(weights)

  normweights <- weights
  for(i in 1:ncol(normweights)){
    normweights[,i] <- normweights[,i]/sum(normweights[,i])
  }


  rnacellconts <- ensemblecellcontents(cellcontlist = cellcontlist,
                                       normweights = normweights)



  methylcellconts <- methylpredict(normweights = normweights,
                                   modellist = modellist,
                                   targetmethyldat = methylmat,
                                   resscale = resscale,
                                   adjustminus = TRUE)

  res <- list(rnacellconts = rnacellconts,
              methylcellconts = methylcellconts,
              modellist = modellist,
              modelcoeflist = modelcoeflist,
              normweights = normweights,
              weights = weights,
              rnamethylsqrs = rnamethylsqrs,
              rnageneidxlist = rnageneidxlist)


  if(!is.null(targetmethyldat)){

    methyltargetcellcounts <- methylpredict(model = res,
                                            targetmethyldat = targetmethyldat,
                                            #deconvmod = deconvmod,
                                            resscale = resscale,
                                            adjustminus = TRUE)

    res$methyltargetcellcounts <- methyltargetcellcounts

  }

  if(plot == TRUE){

    celldeconvplots(cellcontres = res$rnacellconts,
                    pddat = pddat, title = 'RNA')

    celldeconvplots(cellcontres = res$methylcellconts,
                    pddat = pddat, title = 'Methylation')

    if(!is.null(targetmethyldat)){

      celldeconvplots(cellcontres = res$methyltargetcellcounts,
                      pddat = targetmethylpddat, title = NULL)

    }

    cellscatterplots(dat1 = res$rnacellconts,
                     dat2 = res$methylcellconts,
                     title = 'RNA and Methylation',
                     colorful = TRUE)



  }

  return(res)


}



