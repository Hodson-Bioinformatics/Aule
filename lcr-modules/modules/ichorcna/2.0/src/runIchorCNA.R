# file:   ichorCNA.R
# authors: Gavin Ha, Ph.D.
#          Fred Hutch
# contact: <gha@fredhutch.org>
#
#         Justin Rhoades
#          Broad Institute
# contact: <rhoades@broadinstitute.org>

# ichorCNA: https://github.com/broadinstitute/ichorCNA
# date:   July 24, 2019
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

library(optparse)

option_list <- list(
  make_option(c("--WIG"), type = "character", help = "Path to tumor WIG file. Required."),
  make_option(c("--NORMWIG"), type = "character", default=NULL, help = "Path to normal WIG file. Default: [%default]"),
  make_option(c("--gcWig"), type = "character", help = "Path to GC-content WIG file; Required"),
  make_option(c("--mapWig"), type = "character", default=NULL, help = "Path to mappability score WIG file. Default: [%default]"),
  make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals. Default: [%default]"),
  make_option(c("--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions. Default: [%default]"),
  make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
  make_option(c("--centromere"), type="character", default=NULL, help = "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package. Default: [%default]"),
  make_option(c("--minMapScore"), type = "numeric", default=0.9, help="Include bins with a minimum mappability score of this value. Default: [%default]."),
  make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
  make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
  make_option(c("--scStates"), type="character", default="NULL", help = "Subclonal states to consider. Default: [%default]"),
  make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage. Default: [%default]"),
  make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data. Default: [%default]"),
  make_option(c("--lambdaScaleHyperParam"), type="numeric", default=3, help="Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [%default]"),
  #	make_option(c("--kappa"), type="character", default=50, help="Initial state distribution"),
  make_option(c("--ploidy"), type="character", default="2", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
  make_option(c("--maxCN"), type="numeric", default=7, help = "Total clonal CN states. Default: [%default]"),
  make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal. Default: [%default]"),
  make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),
  make_option(c("--estimatePloidy"), type="logical", default=TRUE, help = "Estimate tumour ploidy. Default: [%default]"),
  make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value. Default: [%default]"),
  make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value. Default: [%default]"),
  make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
  make_option(c("--altFracThreshold"), type="numeric", default=0.05, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [%default]"),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases. Default: [%default]"),
  make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze. Default: [%default]"),
  make_option(c("--genomeBuild"), type="character", default="hg19", help="Geome build. Default: [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--normalizeMaleX"), type="logical", default=TRUE, help = "If male, then normalize chrX by median. Default: [%default]"),
  make_option(c("--minTumFracToCorrect"), type="numeric", default=0.1, help = "Tumor-fraction correction of bin and segment-level CNA if sample has minimum estimated tumor fraction. [Default: %default]"), 
  make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
  make_option(c("--includeHOMD"), type="logical", default=FALSE, help="If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [%default]"),
  make_option(c("--txnE"), type="numeric", default=0.9999999, help = "Self-transition probability. Increase to decrease number of segments. Default: [%default]"),
  make_option(c("--txnStrength"), type="numeric", default=1e7, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
  make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots. Default: [%default]"),
  make_option(c("--plotYLim"), type="character", default="c(-2,2)", help = "ylim to use for chromosome plots. Default: [%default]"),
  make_option(c("--outDir"), type="character", default="./", help = "Output Directory. Default: [%default]"),
  make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

library(HMMcopy)
library(GenomicRanges)
library(GenomeInfoDb)
library(foreach)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

patientID <- opt$id
tumour_file <- opt$WIG
normal_file <- opt$NORMWIG
gcWig <- opt$gcWig
mapWig <- opt$mapWig
normal_panel <- opt$normalPanel
exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
minMapScore <- opt$minMapScore
flankLength <- opt$rmCentromereFlankLength
normal <- eval(parse(text = opt$normal))
scStates <- eval(parse(text = opt$scStates))
lambda <- eval(parse(text = opt$lambda))
lambdaScaleHyperParam <- opt$lambdaScaleHyperParam
estimateNormal <- opt$estimateNormal
estimatePloidy <- opt$estimatePloidy
estimateScPrevalence <- opt$estimateScPrevalence
maxFracCNASubclone <- opt$maxFracCNASubclone
maxFracGenomeSubclone <- opt$maxFracGenomeSubclone
minSegmentBins <- opt$minSegmentBins
altFracThreshold <- opt$altFracThreshold
ploidy <- eval(parse(text = opt$ploidy))
coverage <- opt$coverage
maxCN <- opt$maxCN
txnE <- opt$txnE
txnStrength <- opt$txnStrength
normalizeMaleX <- as.logical(opt$normalizeMaleX)
includeHOMD <- as.logical(opt$includeHOMD)
minTumFracToCorrect <- opt$minTumFracToCorrect
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale
chrXMedianForMale <- -0.1
outDir <- opt$outDir
libdir <- opt$libdir
plotFileType <- opt$plotFileType
plotYLim <- eval(parse(text=opt$plotYLim))
gender <- NULL
outImage <- paste0(outDir,"/", patientID,".RData")
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
chrs <- as.character(eval(parse(text = opt$chrs)))
chrTrain <- as.character(eval(parse(text=opt$chrTrain))); 
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize))); 
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle
seqlevelsStyle(chrTrain) <- genomeStyle

## load ichorCNA library or source R scripts
if (!is.null(libdir) && libdir != "None"){
  source(paste0(libdir,"/R/utils.R"))
  source(paste0(libdir,"/R/segmentation.R"))
  source(paste0(libdir,"/R/EM.R"))
  source(paste0(libdir,"/R/output.R"))
  source(paste0(libdir,"/R/plotting.R"))
} else {
  library(ichorCNA)
}



# Modified HMMsegment & runEM functions to fix errors with lambda

HMMsegment <- function (x, validInd = NULL, dataType = "copy", param = NULL, 
                        chrTrain = c(1:22), maxiter = 50, estimateNormal = TRUE, 
                        estimatePloidy = TRUE, estimatePrecision = TRUE, estimateVar = TRUE, 
                        estimateSubclone = TRUE, estimateTransition = TRUE, estimateInitDist = TRUE, 
                        logTransform = FALSE, likChangeConvergence = 0.001, verbose = TRUE) {
  chr <- as.factor(seqnames(x[[1]]))
  dataMat <- as.matrix(as.data.frame(lapply(x, function(y) {
    mcols(y)[, dataType]
  })))
  if (logTransform) {
    dataMat <- apply(dataMat, 2, function(x) {
      log(x/median(x, na.rm = TRUE))
    })
  }
  else {
    dataMat <- log(2^dataMat)
  }
  for (i in 1:length(x)) {
    mcols(x[[i]])[, dataType] <- dataMat[, i]
  }
  if (!is.null(chrTrain)) {
    chrInd <- chr %in% chrTrain
  }
  else {
    chrInd <- !logical(length(chr))
  }
  if (!is.null(validInd)) {
    chrInd <- chrInd & validInd
  }
  if (is.null(param)) {
    param <- getDefaultParameters(dataMat[chrInd])
  }
  convergedParams <- runEM(dataMat, chr, chrInd, param, maxiter, 
                           verbose, estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, 
                           estimateVar = estimateVar, estimateSubclone = estimateSubclone, 
                           estimatePrecision = estimatePrecision, estimateTransition = estimateTransition, 
                           estimateInitDist = estimateInitDist, likChangeConvergence = likChangeConvergence)
  viterbiResults <- runViterbi(convergedParams, chr)
  segs <- segmentData(x, validInd, viterbiResults$states, convergedParams)
  names <- c("HOMD", "HETD", "NEUT", "GAIN", "AMP", "HLAMP", 
             paste0("HLAMP", 2:1000))
  cnaList <- list()
  S <- length(x)
  for (s in 1:S) {
    id <- names(x)[s]
    copyNumber <- param$jointCNstates[viterbiResults$state, 
                                      s]
    subclone.status <- param$jointSCstatus[viterbiResults$state, 
                                           s]
    cnaList[[id]] <- data.frame(cbind(sample = as.character(id), 
                                      chr = as.character(seqnames(x[[s]])), start = start(x[[s]]), 
                                      end = end(x[[s]]), copy.number = copyNumber, event = names[copyNumber + 
                                                                                                   1], logR = round(log2(exp(dataMat[, s])), digits = 4), 
                                      subclone.status = as.numeric(subclone.status)))
    cnaList[[id]] <- transform(cnaList[[id]], start = as.integer(as.character(start)), 
                               end = as.integer(as.character(end)), copy.number = as.numeric(copy.number), 
                               logR = as.numeric(as.character(logR)), subclone.status = as.numeric(subclone.status))
    chrOrder <- unique(chr)
    cnaList[[id]] <- cnaList[[id]][order(match(cnaList[[id]][, 
                                                             "chr"], chrOrder)), ]
    cnaList[[id]] <- cnaList[[id]][cnaList[[id]][, "chr"] %in% 
                                     chrOrder, ]
    segs[[s]]$median <- log2(exp(segs[[s]]$median))
    segs[[s]]$subclone.status <- param$jointSCstatus[segs[[s]]$state, 
                                                     s]
  }
  convergedParams$segs <- segs
  return(list(cna = cnaList, results = convergedParams, viterbiResults = viterbiResults))
}

runEM <- function (copy, chr, chrInd, param, maxiter, verbose = TRUE, 
                   estimateNormal = TRUE, estimatePloidy = TRUE, estimateVar = TRUE, 
                   estimatePrecision = TRUE, estimateTransition = TRUE, estimateInitDist = TRUE, 
                   estimateSubclone = TRUE, likChangeConvergence = 0.001) {
  if (nrow(copy) != length(chr) || nrow(copy) != length(chrInd)) {
    stop("runEM: Length of inputs do not match for one of: copy, chr, chrInd")
  }
  if (is.null(param$ct) || is.null(param$lambda) || is.null(param$nu) || 
      is.null(param$kappa)) {
    stop("runEM: Parameter missing, ensure all parameters exist as columns in", 
         "data frame: ct, lambda, nu, kappa")
  }
  S <- param$numberSamples
  K <- length(param$ct)
  Z <- sum(param$ct.sc)
  KS <- nrow(param$jointStates)
  N <- nrow(copy)
  rho <- matrix(0, KS, N)
  lambdas <- array(0, dim = c(length(na.omit(param$lambda)), ### change nrow to length
                              S, maxiter))
  vars <- array(0, dim = c(KS, S, maxiter))
  phi <- matrix(NA, length(param$phi_0), maxiter)
  n <- matrix(NA, S, maxiter)
  sp <- matrix(NA, S, maxiter)
  piG <- matrix(0, KS, maxiter)
  A <- param$A
  converged <- FALSE
  Zcounts <- matrix(0, KS, KS)
  loglik <- rep(0, maxiter)
  mus <- array(0, dim = c(KS, S, maxiter))
  py <- matrix(0, KS, N)
  ptmTotal <- proc.time()
  chrs <- levels(chr)
  chrsI <- vector("list", length(chrs))
  for (i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
  }
  if (verbose) {
    message("runEM: Initialization")
  }
  i <- 1
  piG[, i] <- normalize(param$kappa)
  n[, i] <- param$n_0
  sp[, i] <- param$sp_0
  phi[, i] <- param$phi_0
  lambdas[, , i] <- na.omit(param$lambda)
  lambdasKS <- as.matrix(na.omit(expand.grid(as.data.frame(lambdas[, 
                                                                   , i]))))
  if (param$likModel == "Gaussian") {
    vars[, , i] <- as.matrix(na.omit(expand.grid(as.data.frame(param$var))))
    varsKS <- as.matrix(vars[, , i])
  }
  mus[, , i] <- as.matrix(get2and3ComponentMixture(param$jointCNstates, 
                                                   param$jointSCstatus, n[, i], sp[, i], phi[, i]))
  if (param$likModel == "t") {
    for (ks in 1:KS) {
      probs <- tdistPDF(copy, mus[ks, , i], lambdasKS[ks, 
      ], param$nu)
      py[ks, ] <- apply(probs, 1, prod)
    }
  }
  else if (param$likModel == "Gaussian") {
    message("Using ", param$likModel, " emission model.")
    py <- getNormLik(copy, mus[, , i, drop = FALSE], varsKS, 
                     param$sw)
  }
  loglik[i] <- -Inf
  while (!converged && (i < maxiter)) {
    ptm <- proc.time()
    i <- i + 1
    if (verbose) {
      message("runEM iter", i - 1, ": Expectation")
    }
    if (getDoParWorkers() > 1) {
      fwdBack <- foreach(j = 1:length(chrsI), .combine = rbind, 
                         .noexport = c("copy")) %dopar% {
                           if (verbose == TRUE) {
                             message(j, " ", appendLF = FALSE)
                           }
                           I <- intersect(chrsI[[j]], which(chrInd))
                           if (length(I) > 0) {
                             .Call("forward_backward", piG[, i - 1], A, 
                                   py[, I], PACKAGE = "HMMcopy")
                           }
                         }
      if (verbose == TRUE) {
        message("")
      }
      if (length(chrsI) > 1) {
        loglik[i] <- sum(do.call(rbind, fwdBack[, "loglik"]))
        rho[, chrInd] <- do.call(cbind, fwdBack[, "rho"])
        tmp <- lapply(as.list(fwdBack[, "xi"]), function(x) {
          rowSums(x, dims = 2)
        })
        Zcounts <- Reduce("+", tmp)
      }
      else {
        loglik[i] <- fwdBack[["loglik"]]
        rho[, chrInd] <- fwdBack[["rho"]]
        Zcounts <- rowSums(fwdBack[["xi"]], dims = 2)
      }
    }
    else {
      Zcounts <- matrix(0, KS, KS)
      for (j in 1:length(chrsI)) {
        I <- intersect(chrsI[[j]], which(chrInd))
        if (length(I) > 0) {
          output <- .Call("forward_backward", piG[, i - 
                                                    1], A, py[, I], PACKAGE = "HMMcopy")
          rho[, I] <- output$rho
          loglik[i] <- loglik[i] + output$loglik
          Zcounts <- Zcounts + t(colSums(aperm(output$xi, 
                                               c(3, 2, 1))))
        }
      }
    }
    elapsedTime <- proc.time() - ptm
    message("runEM iter ", i - 1, ": E-step ", format(elapsedTime[3]/60, 
                                                      digits = 2), " min.")
    if (verbose) {
      message("runEM iter", i - 1, ": Maximization")
    }
    ptm.mstep <- proc.time()
    if (param$likModel == "Gaussian") {
      output <- estimateGaussianParamsMap(copy[chrInd, 
      ], n[, i - 1], sp[, i - 1], phi[, i - 1], varsKS, 
      piG[, i - 1], A, param, rho[, chrInd], Zcounts, 
      estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, 
      estimateVar = estimateVar, estimateTransition = estimateTransition, 
      estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
      vars[, , i] <- output$var
      if (verbose == TRUE) {
        for (s in 1:S) {
          message("Sample", s, " n=", signif(output$n[s], 
                                             digits = 4), ", sp=", signif(output$sp[s], 
                                                                          digits = 4), ", phi=", signif(output$phi[s], 
                                                                                                        digits = 4))
        }
      }
    }
    else {
      output <- estimateTParamsMap(copy[chrInd, ], n[, 
                                                     i - 1], sp[, i - 1], phi[, i - 1], lambdas[, 
                                                                                                , i - 1], piG[, i - 1], A, param, rho[, chrInd], 
                                   Zcounts, loglik[i], estimateNormal = estimateNormal, 
                                   estimatePloidy = estimatePloidy, estimatePrecision = estimatePrecision, 
                                   estimateTransition = estimateTransition, estimateInitDist = estimateInitDist, 
                                   estimateSubclone = estimateSubclone)
      lambdas[, , i] <- output$lambda
      lambdasKS <- as.matrix(na.omit(expand.grid(as.data.frame(lambdas[, 
                                                                       , i]))))
      if (verbose == TRUE) {
        for (s in 1:S) {
          message("Sample", s, " n=", signif(output$n[s], 
                                             digits = 4), ", sp=", signif(output$sp[s], 
                                                                          digits = 4), ", phi=", signif(output$phi[s], 
                                                                                                        digits = 4), ", lambda=", paste0(signif(output$lambda, 
                                                                                                                                                digits = 5), collapse = ","), ", F=", signif(output$F, 
                                                                                                                                                                                             digits = 4))
        }
      }
    }
    elapsedTime <- proc.time() - ptm.mstep
    message("runEM iter ", i - 1, ": M-step ", format(elapsedTime[3]/60, 
                                                      digits = 2), " min.")
    n[, i] <- output$n
    sp[, i] <- output$sp
    phi[, i] <- output$phi
    piG[, i] <- output$piG
    A <- output$A
    estF <- output$F
    varsKS <- as.matrix(vars[, , i])
    mus[, , i] <- as.matrix(get2and3ComponentMixture(param$jointCNstates, 
                                                     param$jointSCstatus, n[, i], sp[, i], phi[, i]))
    if (param$likModel == "t") {
      for (ks in 1:KS) {
        probs <- tdistPDF(copy, mus[ks, , i], lambdasKS[ks, 
        ], param$nu)
        py[ks, ] <- apply(probs, 1, prod)
      }
    }
    else if (param$likModel == "Gaussian") {
      py <- getNormLik(copy, mus[, , i, drop = FALSE], 
                       varsKS, param$sw)
    }
    prior <- priorProbs(n[, i], sp[, i], phi[, i], lambdas[, 
                                                           , i], varsKS, piG[, i], A, param, estimateNormal = estimateNormal, 
                        estimatePloidy = estimatePloidy, estimateVar = estimateVar, 
                        estimatePrecision = estimatePrecision, estimateTransition = estimateTransition, 
                        estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
    loglik[i] <- loglik[i] + prior$prior
    elapsedTime <- proc.time() - ptm
    if (verbose) {
      message(paste("runEM iter", i - 1, " Log likelihood:", 
                    loglik[i]))
      message("runEM iter", i - 1, " Time: ", format(elapsedTime[3]/60, 
                                                     digits = 2), " min.")
    }
    if ((abs(loglik[i] - loglik[i - 1])/abs(max(loglik[i - 
                                                       1], -1e+30))) < likChangeConvergence) {
      message("runEM iter", i - 1, " EM Converged")
      converged = 1
    }
    if (loglik[i] < loglik[i - 1]) {
      message("Using previous iteration ", i - 2)
      i <- i - 1
      converged = 1
    }
  }
  if (converged) {
    if (verbose) {
      message("runEM iter", i - 1, ": Re-calculating responsibilties from converged parameters.")
    }
    for (j in 1:length(chrsI)) {
      I <- chrsI[[j]]
      output <- .Call("forward_backward", piG[, i], A, 
                      py[, I], PACKAGE = "HMMcopy")
      rho[, I] <- output$rho
    }
  }
  if (verbose) {
    totalTime <- proc.time() - ptmTotal
    message("runEM: Using optimal parameters from iter", 
            i - 1)
    message("runEM: Total elapsed time: ", format(totalTime[3]/60, 
                                                  digits = 2), "min.")
  }
  n <- n[, 1:i, drop = FALSE]
  sp <- sp[, 1:i, drop = FALSE]
  phi <- phi[, 1:i, drop = FALSE]
  mus <- mus[, , 1:i, drop = FALSE]
  vars <- vars[, , 1:i, drop = FALSE]
  lambdas <- lambdas[, , 1:i, drop = FALSE]
  piG <- piG[, 1:i, drop = FALSE]
  loglik = loglik[1:i]
  output <- vector("list", 0)
  output$n <- n
  output$sp <- sp
  output$phi <- phi
  output$mus <- mus
  output$vars <- varsKS
  output$lambdas <- lambdas
  output$pi <- piG
  output$A <- A
  output$loglik <- loglik
  output$rho <- rho
  output$param <- param
  output$py <- py
  output$iter <- i
  return(output)
}



## load seqinfo 
# seqinfo <- getSeqInfo(genomeBuild, genomeStyle)
seqinfo <- NULL

if (substr(tumour_file,nchar(tumour_file)-2,nchar(tumour_file)) == "wig") {
  wigFiles <- data.frame(cbind(patientID, tumour_file))
} else {
  wigFiles <- read.delim(tumour_file, header=F, as.is=T)
}

## FILTER BY EXONS IF PROVIDED ##
## add gc and map to GRanges object ##
if (is.null(exons.bed) || exons.bed == "None" || exons.bed == "NULL"){
  targetedSequences <- NULL
}else{
  targetedSequences <- read.delim(exons.bed, header=T, sep="\t")  
}

## load PoN
if (is.null(normal_panel) || normal_panel == "None" || normal_panel == "NULL"){
  normal_panel <- NULL
}

if (is.null(centromere) || centromere == "None" || centromere == "NULL"){ # no centromere file provided
  centromere <- system.file("extdata", "GRCh37.p13_centromere_UCSC-gapTable.txt", 
                            package = "ichorCNA")
}
centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
save.image(outImage)
## LOAD IN WIG FILES ##
numSamples <- nrow(wigFiles)

tumour_copy <- list()
for (i in 1:numSamples) {
  id <- wigFiles[i,1]
  ## create output directories for each sample ##
  dir.create(paste0(outDir, "/", id, "/"), recursive = TRUE)
  ### LOAD TUMOUR AND NORMAL FILES ###
  message("Loading tumour file:", wigFiles[i,1])
  tumour_reads <- wigToGRanges(wigFiles[i,2])
  
  ## LOAD GC/MAP WIG FILES ###
  # find the bin size and load corresponding wig files #
  binSize <- as.data.frame(tumour_reads[1,])$width 
  message("Reading GC and mappability files")
  if (is.null(gcWig) || gcWig == "None" || gcWig == "NULL"){
    stop("GC wig file is required")
  }
  gc <- wigToGRanges(gcWig)
  if (is.null(mapWig) || mapWig == "None" || mapWig == "NULL"){
    message("No mappability wig file input, excluding from correction")
    map <- NULL
  } else {
    map <- wigToGRanges(mapWig)
  }
  message("Correcting Tumour")
  
  counts <- loadReadCountsFromWig(tumour_reads, chrs = chrs, gc = gc, map = map, 
                                  centromere = centromere, flankLength = flankLength, 
                                  targetedSequences = targetedSequences, chrXMedianForMale = chrXMedianForMale,
                                  genomeStyle = genomeStyle, fracReadsInChrYForMale = fracReadsInChrYForMale,
                                  chrNormalize = chrNormalize, mapScoreThres = minMapScore)
  tumour_copy[[id]] <- counts$counts #as(counts$counts, "GRanges")
  gender <- counts$gender
  ## load in normal file if provided 
  if (!is.null(normal_file) && normal_file != "None" && normal_file != "NULL"){
    message("Loading normal file:", normal_file)
    normal_reads <- wigToGRanges(normal_file)
    message("Correcting Normal")
    counts <- loadReadCountsFromWig(normal_reads, chrs=chrs, gc=gc, map=map, 
                                    centromere=centromere, flankLength = flankLength, targetedSequences=targetedSequences,
                                    genomeStyle = genomeStyle, chrNormalize = chrNormalize, mapScoreThres = minMapScore)
    normal_copy <- counts$counts #as(counts$counts, "GRanges")
    gender.normal <- counts$gender
  }else{
    normal_copy <- NULL
  }
  
  ### DETERMINE GENDER ###
  ## if normal file not given, use chrY, else use chrX
  message("Determining gender...", appendLF = FALSE)
  gender.mismatch <- FALSE
  if (!is.null(normal_copy)){
    if (gender$gender != gender.normal$gender){ #use tumour # use normal if given
      # check if normal is same gender as tumour
      gender.mismatch <- TRUE
    }
  }
  message("Gender ", gender$gender)
  
  ## NORMALIZE GENOME-WIDE BY MATCHED NORMAL OR NORMAL PANEL (MEDIAN) ##
  tumour_copy[[id]] <- normalizeByPanelOrMatchedNormal(tumour_copy[[id]], chrs = chrs, 
                                                       normal_panel = normal_panel, normal_copy = normal_copy, 
                                                       gender = gender$gender, normalizeMaleX = normalizeMaleX)
  
  ### OUTPUT FILE ###
  ### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
  outMat <- as.data.frame(tumour_copy[[id]])
  #outMat <- outMat[,c(1,2,3,12)]
  outMat <- outMat[,c("seqnames","start","end","copy")]
  colnames(outMat) <- c("chr","start","end","log2_TNratio_corrected")
  outFile <- paste0(outDir,"/",id,".correctedDepth.txt")
  message(paste("Outputting to:", outFile))
  write.table(outMat, file=outFile, row.names=F, col.names=T, quote=F, sep="\t")
  
} ## end of for each sample

chrInd <- as.character(seqnames(tumour_copy[[1]])) %in% chrTrain
## get positions that are valid
valid <- tumour_copy[[1]]$valid
if (length(tumour_copy) >= 2) {
  for (i in 2:length(tumour_copy)){ 
    valid <- valid & tumour_copy[[i]]$valid 
  } 
}
save.image(outImage)

### RUN HMM ###
## store the results for different normal and ploidy solutions ##
ptmTotalSolutions <- proc.time() # start total timer
results <- list()
loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
                               dimnames = list(c(), c("init", "n_est", "phi_est", "BIC", 
                                                      "Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
counter <- 1
compNames <- rep(NA, nrow(loglik))
mainName <- rep(NA, length(normal) * length(ploidy))
#### restart for purity and ploidy values ####
for (n in normal){
  for (p in ploidy){
    if (n == 0.95 & p != 2) {
      next
    }
    logR <- as.data.frame(lapply(tumour_copy, function(x) { x$copy })) # NEED TO EXCLUDE CHR X #
    param <- getDefaultParameters(logR[valid & chrInd, , drop=F], maxCN = maxCN, includeHOMD = includeHOMD, 
                                  ct.sc=scStates, ploidy = floor(p), e=txnE, e.same = 50, strength=txnStrength)
    param$phi_0 <- rep(p, numSamples)
    param$n_0 <- rep(n, numSamples)
    
    ############################################
    ######## CUSTOM PARAMETER SETTINGS #########
    ############################################
    # 0.1x cfDNA #
    if (is.null(lambda)){
      logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
      param$lambda <- rep(logR.var, length(param$ct))
      param$lambda[param$ct %in% c(2)] <- logR.var 
      param$lambda[param$ct %in% c(1,3)] <- logR.var 
      param$lambda[param$ct >= 4] <- logR.var / 5
      param$lambda[param$ct == max(param$ct)] <- logR.var / 15
      param$lambda[param$ct.sc.status] <- logR.var / 10
    }else{
      param$lambda[param$ct %in% c(2)] <- lambda[2]
      param$lambda[param$ct %in% c(1)] <- lambda[1]
      param$lambda[param$ct %in% c(3)] <- lambda[3]
      param$lambda[param$ct >= 4] <- lambda[4]
      param$lambda[param$ct == max(param$ct)] <- lambda[2] / 15
      param$lambda[param$ct.sc.status] <- lambda[2] / 10
    }
    param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct))  
    # 1x bulk tumors #
    #param$lambda[param$ct %in% c(2)] <- 2000
    #param$lambda[param$ct %in% c(1)] <- 1750
    #param$lambda[param$ct %in% c(3)] <- 1750
    #param$lambda[param$ct >= 4] <- 1500
    #param$lambda[param$ct == max(param$ct)] <- 1000 / 25
    #param$lambda[param$ct.sc.status] <- 1000 / 75
    #param$alphaLambda[param$ct.sc.status] <- 4
    #param$alphaLambda[param$ct %in% c(1,3)] <- 5
    #param$alphaLambda[param$ct %in% c(2)] <- 5
    #param$alphaLambda[param$ct == max(param$ct)] <- 4
    
    #############################################
    ################ RUN HMM ####################
    #############################################
    hmmResults.cor <- HMMsegment(tumour_copy, valid, dataType = "copy", 
                                 param = param, chrTrain = chrTrain, maxiter = 50,
                                 estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                 estimateSubclone = estimateScPrevalence, verbose = TRUE)
    
    for (s in 1:numSamples){
      iter <- hmmResults.cor$results$iter
      id <- names(hmmResults.cor$cna)[s]
      
      ## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
      ## check if there is an altered segment that has at least a minimum # of bins
      segsS <- hmmResults.cor$results$segs[[s]]
      segsS <- segsS[segsS$chr %in% chrTrain, ]
      segAltInd <- which(segsS$event != "NEUT")
      maxBinLength = -Inf
      if (sum(segAltInd) > 0){
        maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
        maxSegRD <- GRanges(seqnames=segsS$chr[segAltInd[maxInd]], 
                            ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
        hits <- findOverlaps(query=maxSegRD, subject=tumour_copy[[s]][valid, ])
        maxBinLength <- length(subjectHits(hits))
      }
      ## check if there are proportion of total bins altered 
      # if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
      cnaS <- hmmResults.cor$cna[[s]]
      altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
      altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
      if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
        hmmResults.cor$results$n[s, iter] <- 1.0
      }
      
      # correct integer copy number based on estimated purity and ploidy
      correctedResults <- correctIntegerCN(cn = hmmResults.cor$cna[[s]],
                                           segs = hmmResults.cor$results$segs[[s]], 
                                           purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
                                           cellPrev = 1 - hmmResults.cor$results$sp[s, iter], 
                                           maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = minTumFracToCorrect, 
                                           gender = gender$gender, chrs = chrs, correctHOMD = includeHOMD)
      hmmResults.cor$results$segs[[s]] <- correctedResults$segs
      hmmResults.cor$cna[[s]] <- correctedResults$cn
      
      ## plot solution ##
      outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_", "n", n, "-p", p)
      mainName[counter] <- paste0(id, ", n: ", n, ", p: ", p, ", log likelihood: ", signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4))
      plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType=plotFileType, 
                     logR.column = "logR", call.column = "Corrected_Call",
                     plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, seqinfo=seqinfo, main=mainName[counter])
    }
    iter <- hmmResults.cor$results$iter
    results[[counter]] <- hmmResults.cor
    loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
    subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
    fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
    fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
    fracAltSub <- lapply(fracAltSub, function(x){if (is.na(x)){0}else{x}})
    loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits=2), collapse=",")
    loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits=2), collapse=",")
    loglik[counter, "init"] <- paste0("n", n, "-p", p)
    loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
    loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")
    
    counter <- counter + 1
  }
}
## get total time for all solutions ##
elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

### SAVE R IMAGE ###
save.image(outImage)
#save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))

### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
loglik <- loglik[!is.na(loglik$init), ]
if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
  fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
                     loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
  if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
    ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
  }else{ # otherwise just take largest likelihood
    ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
  }
}else{#sort by likelihood only
  ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
}

#new loop by order of solutions (ind)
outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_all_sols")
for(i in 1:length(ind)) {
  hmmResults.cor <- results[[ind[i]]]
  turnDevOff <- FALSE
  turnDevOn <- FALSE
  if (i == 1){
    turnDevOn <- TRUE
  }
  if (i == length(ind)){
    turnDevOff <- TRUE
  }
  plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
                 logR.column = "logR", call.column = "Corrected_Call",
                 plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
                 seqinfo = seqinfo,
                 turnDevOn = turnDevOn, turnDevOff = turnDevOff, main=mainName[ind[i]])
}

hmmResults.cor <- results[[ind[1]]]
hmmResults.cor$results$loglik <- as.data.frame(loglik)
hmmResults.cor$results$gender <- gender$gender
hmmResults.cor$results$chrYCov <- gender$chrYCovRatio
hmmResults.cor$results$chrXMedian <- gender$chrXMedian
hmmResults.cor$results$coverage <- coverage

outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
          results = hmmResults.cor$results, patientID = patientID, outDir=outDir)
outFile <- paste0(outDir, "/", patientID, ".params.txt")
outputParametersToFile(hmmResults.cor, file = outFile)

## plot solutions for all samples 
plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir, counts, numSamples=numSamples,
              logR.column = "logR", call.column = "Corrected_Call",
              plotFileType=plotFileType, plotYLim=plotYLim, seqinfo = seqinfo,
              estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)