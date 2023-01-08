#=======================================================
#' @title Simulate Genotypes
#'
#' @description Simulate SNP genotype data from a pedigree, with optional
#'   missingess and errors.
#'
#' @details Please ensure the pedigree is a valid pedigree, for example by first
#'   running \code{\link{PedPolish}}. For founders, i.e. individuals with no
#'   known parents, genotypes are drawn according to the provided MAF and
#'   assuming Hardy-Weinberg equilibrium. Offspring genotypes are generated
#'   following Mendelian inheritance, assuming all loci are completely
#'   independent. Individuals with one known parent are allowed: at each locus,
#'   one allele is inherited from the known parent, and the other drawn from the
#'   genepool according to the provided MAF.
#'
#'   Genotyping errors are generated following a user-definable 3x3 matrix with
#'   probabilities that actual genotype \eqn{i} (rows) is observed as genotype
#'   \eqn{j} (columns). This is specified as \code{ErrorFM}, which is a function
#'   of \code{SnpError}. By default (\code{ErrorFM} = "version2.0"),
#'   \code{SnpError} is interpreted as a locus-level error rate (rather than
#'   allele-level), and equals the probability that a homozygote is observed as
#'   heterozygote, and the probability that a heterozygote is observed as either
#'   homozygote (i.e., the probability that it is observed as AA = probability
#'   that observed as aa = \code{SnpError}/2). The probability that one
#'   homozygote is observed as the other is (\code{SnpError}/2\eqn{)^2}.
#'
#'   Note that this differs from versions up to 1.1.1, where a proportion of
#'   \code{SnpError}*3/2 of genotypes were replaced with random genotypes. This
#'   corresponds to \code{ErrorFM} = "Version111".
#'
#'   Error rates differ between SNPs, but the same error pattern is used across
#'   all SNPs, even when inheritance patterns vary. When two or more different
#'   error patterns are required, SimGeno should be run on the different SNP
#'   subsets separately, and results combined.
#'
#'   Variation in call rates is assumed to follow a highly skewed (beta)
#'   distribution, with many samples having call rates close to 1, and a
#'   narrowing tail of lower call rates. The first shape parameter defaults to 1
#'   (but see \code{\link{MkGenoErrors}}), and the second shape parameter is
#'   defined via the mean as \code{CallRate}. For 99.9\% of SNPs to have a call
#'   rate of 0.8 (0.9; 0.95) or higher, use a mean call rate of 0.969 (0.985;
#'   0.993).
#'
#'   Variation in call rate between samples can be specified by providing a
#'   named vector to \code{CallRate}, which supersedes PropLQ in versions up to
#'   1.1.1. Otherwise, variation in call rate and error rate between samples
#'   occurs only as side-effect of the random nature of which individuals are
#'   hit by per-SNP errors and drop-outs. Finer control is possible by first
#'   generating an error-free genotype matrix, and then calling
#'   \code{\link{MkGenoErrors}} directly on subsets of the matrix.
#'
#' @param Pedigree  dataframe, pedigree with the first three columns being id -
#'   dam - sire. Column names are ignored, as are additional columns, with the
#'   exception of a 'Sex' column when Inherit is not 'autosomal'.
#' @param nSnp  number of SNPs to simulate.
#' @param ParMis  single number or vector length two with proportion of parents
#'   with fully missing genotype. Ignored if CallRate is a named vector.
#' @param MAF  minimum minor allele frequency, and allele frequencies will be
#'   sampled uniformly between this minimum and 0.5, OR a vector with minor
#'   allele frequency at each locus. In both cases, this is the MAF among
#'   pedigree founders, the MAF in the sample will deviate due to drift.
#' @param CallRate either a single number for the mean call rate (genotyping
#'   success), OR a vector with the call rate at each SNP, OR a named vector
#'   with the call rate for each individual. In the third case, ParMis is
#'   ignored, and individuals in the pedigree (as id or parent) not included in
#'   this vector are presumed non-genotyped.
#' @param SnpError  mean per-locus genotyping error rate across SNPs, and a
#'   beta-distribution will be used to simulate the number of missing cases per
#'   SNP, OR a vector with the genotyping error for each SNP.
#' @param ErrorFM  function taking the error rate (scalar) as argument and
#'   returning a 3x3 matrix with probabilities that actual genotype i (rows) is
#'   observed as genotype j (columns). Inbuilt ones are as used in sequoia
#'   'version2.0', 'version1.3', or 'version1.1'. See details.
#' @param ReturnStats in addition to the genotype matrix, return the input
#'   parameters and mean & quantiles of MAF, error rate and call rates.
#' @param OutFile  file name for simulated genotypes. If NA (default), return
#'   results within R.
#' @param Inherit  inheritance pattern, scalar or vector of length nSnp,
#'   Defaults to 'autosomal'. An excel file included in the package has
#'   inheritance patterns for the X and Y chromosome and mtDNA, and allows
#'   custom inheritance patterns. Note that these are experimental, and NOT
#'   currently supported by the pedigree reconstruction with
#'   \code{\link{sequoia}} !
#' @param InheritFile  file name of file with inheritance patterns, with
#'   extension csv, txt, xls or xlsx (the latter two require library
#'   \pkg{openxlsx}).
#' @param quiet suppress messages.
#'
#' @return If \code{ReturnStats=FALSE} (the default), a matrix with genotype
#'   data in sequoia's input format, encoded as 0/1/2/-9.
#'
#'   If \code{ReturnStats=TRUE}, a named list with three elements: list
#'   'ParamsIN', matrix 'SGeno', and list 'StatsOUT':
#'   \item{AF}{Frequency in 'observed' genotypes of '1' allele}
#'   \item{AF.act}{Allele frequency in 'actual' (without genotyping errors &
#'     missingness)}
#'   \item{SnpError}{Error rate per SNP (actual /= observed AND observed /=
#'     missing)}
#'   \item{SnpCallRate}{Non-missing per SNP}
#'   \item{IndivError}{Error rate per individual}
#'   \item{IndivCallRate}{Non-missing per individual}
#'
#' @seealso The wrapper \code{\link{EstConf}} for repeated simulation and
#'   pedigree reconstruction; \code{\link{MkGenoErrors}} for fine control over
#'   the distribution of genotyping errors in simulated data.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @section Disclaimer: This simulation is highly simplistic and assumes that
#'   all SNPs segregate completely independently, that the SNPs are in
#'   Hardy-Weinberg equilibrium in the pedigree founders. It assumes that
#'   genotyping errors are not due to heritable mutations of the SNPs, and that
#'   missingness is random and not e.g. due to heritable mutations of SNP
#'   flanking regions. Results based on this simulated data will provide an
#'   minimum estimate of the number of SNPs required, and an optimistic estimate
#'   of pedigree reconstruction performance.
#'
#' @examples
#' GenoM <- SimGeno(Pedigree = Ped_HSg5, nSnp = 100, ParMis = c(0.2, 0.7))
#'
#' \dontrun{
#' # Alternative genotyping error model
#' EFM <- function(E) {   # Whalen, Gorjanc & Hickey 2018
#'  matrix(c(1-E*3/4, E/4, E/4,
#'           E/4, 1/2-E/4, 1/2-E/4, E/4,
#'           E/4, E/4, 1-E*3/4),
#'           3,3, byrow=TRUE)  }
#' EFM(0.01)
#' GenoM <- SimGeno(Pedigree = Ped_HSg5, nSnp = 100, ParMis = 0.2,
#'  SnpError = 5e-3, ErrorFM = EFM)
#'
#' # combination of high & low quality SNPs
#' Geno.HQ <- SimGeno(Ped_HSg5, nSnp=50, MAF=0.3, CallRate=runif(50, 0.7, 1))
#' Geno.LQ <- SimGeno(Ped_HSg5, nSnp=20, MAF=0.1, CallRate=runif(20, 0.1, 5))
#' Geno.HQLQ <- merge(Geno.HQ, Geno.LQ, by="row.names")
#' }
#'
#' @importFrom stats rbinom runif rbeta
#' @importFrom plyr laply
#'
#' @export

SimGeno <- function(Pedigree,
                    nSnp = 400,
                    ParMis = 0.4,
                    MAF = 0.3,
                    CallRate = 0.99,
                    SnpError = 5e-4,
                    ErrorFM = "version2.0",
					          ReturnStats = FALSE,
					          OutFile = NA,
					          Inherit = "autosomal",
					          InheritFile = NA,
					          quiet = FALSE)
{
  if (is.null(OutFile)) stop("'OutFile' must be filename or NA")
  if (missing(Pedigree)) stop("please provide a pedigree to simulate from")


  # unavoidable partial matching when 'Err' is specified instead of 'SnpError'
  if (is.numeric(ErrorFM)) {
    SnpError <- ErrorFM
    ErrorFM <- "version2.0"   # problem doesn't occur when 'ErrorFM' specified; restore default
  }

  #================================
  # check input
  if(!(is.numeric(nSnp)) || nSnp<=0)  stop("nSnp should be a number greater than 0")
  if(length(ParMis)==1) ParMis <- rep(ParMis, 2)

  params <- list(ParMis1 = ParMis[1], ParMis2 = ParMis[2], MAF = MAF,
                 SnpError = SnpError, CallRate = CallRate)
  for (p in seq_along(params)) {
    if (p==1 & length(params[[p]]) != 1) {
      stop("Length ", names(params)[p], " should be 1")
    } else if (p<4 & !length(params[[p]]) %in% c(1, nSnp)) {
      stop("Length ", names(params)[p], " should be 1 or nSnp")
    } else if (is.null(params[[p]]) | all(is.na(params[[p]]))) {
      stop("Please provide ", names(params)[p])
    } else if ( ! (is.numeric(params[[p]]) & all(params[[p]]>=0) & all(params[[p]]<=1)) ) {
      stop(names(params)[p], " must be a number between 0 and 1")
    }
  }

  if (length(CallRate) > 1) {
    if (is.null(names(CallRate))) {
      if (length(CallRate) != nSnp)  stop("CallRate should be length 1, or length nSnp, or a named vector")
    } else {
      if (length(intersect(names(CallRate), Pedigree[,1]))==0) stop("names of CallRate vector do not match pedigree")
    }
  }


  # Error matrix: rows = actual, columns = observed
	ErFunc <- ErrToM(Err=0.1, flavour=ErrorFM, Return="function")  # 0.1 is test value only

  if(!is.na(OutFile) & ReturnStats)  stop("Cannot write Return Stats to OutFile")

  if (interactive() & !quiet & !is.na(OutFile)) {
    if (file.exists(OutFile)) {
      ANS <- readline(prompt = paste("WARNING: ", OutFile,
                                     "will be overwritten.",
                                     "Press <N> to abort, or any other key to continue."))
    } else {
      ANS <- readline(prompt = paste("Genotypes will be written to ", OutFile,
                                     ". Press <N> to abort, or any other key to continue."))
    }
    if (substr(ANS, 1, 1) %in% c("N", "n")) stop()
  }

  ParamsIN <- as.list(environment())

  #================================
  # minor allele frequencies (among founders)
  if (length(MAF)==1) {
    Q <- round(runif(nSnp, min=MAF, max=0.5),3)
  } else {
    Q <- as.numeric(MAF)
  }

  #================================
  # check & prep
  Ped <- PedPolish(Pedigree, ZeroToNA=TRUE)
  nInd <- nrow(Ped)
  if (any(round(Q*nInd) %in% c(0,1)))  warning("some simulated SNPs have fixed alleles")

  #================================
  # simulate genotypes
  # founders: random draw of alleles under HWE
  # non-founders: following Mendelian inheritance

  # rownumber of dam & sire
  Ped$damIDx <- sapply(Ped[,2], function(x) ifelse(is.na(x), 0, which(Ped[,1]==x)))
  Ped$sireIDx <- sapply(Ped[,3], function(x) ifelse(is.na(x), 0, which(Ped[,1]==x)))

  #~~~~~~~~~~~~~~
  # divide pedigree into `generations': the parents of an individual must come
  # from earlier cohorts than itself, or from the founder population (gen 0)
  Ped$gen <- getGenerations(Ped, StopIfInvalid=TRUE)
  nGen <- max(Ped$gen)


  #~~~~~~~~~~~~~~
  if (all(Inherit == "autosomal")) {

    # founders
    SGeno <- matrix(NA, nInd, nSnp,
                    dimnames=list(Ped[,1], NULL))
    for (i in which(Ped$gen==0)) {
      SGeno[i, ] <- rbinom(nSnp, 2, prob=Q)
    }

    # non-founders
    getHaplo <- function(j, G, Q) {
      if (j >0) {
        rbinom(ncol(G), 1, prob=G[j,]/2)
      } else {
        rbinom(length(Q), 1, prob=Q)
      }
    }

    for (g in 1:nGen) {
      for (i in which(Ped$gen==g)) {
        SGeno[i, ] <- rowSums( cbind( getHaplo(Ped$damIDx[i], SGeno, Q),
                                      getHaplo(Ped$sireIDx[i], SGeno, Q) ))
      }
    }

  #~~~~~~~~~~~~~~
  } else {

    if (!"Sex" %in% names(Ped))  stop("Inherit other than autosomal requires 'Sex' column in Ped")
    Ped$Sex[is.na(Ped$Sex)] <- 3
    if (any(Ped$Sex == 4)) stop("non-autosomal inheritance not yet implemented for hermaphrodites (Sex=4)")
    Ped$Sex[!Ped$Sex %in% 1:2] <- 3

    # inheritance patterns
    if (is.na(InheritFile) | is.null(InheritFile)) {
      Inherit_patterns <- NULL
      utils::data(Inherit_patterns, package='sequoia')
      INHA <- Inherit_patterns
    } else {
      INHA <- ReadSpecialInherit(InheritFile, quiet)  # array: inherit - off sex - geno off - dam - sire
    }

    # founders
    founderProp <- apply(INHA, 1:3, sum)
    SGeno.4 <- matrix(NA, nInd, nSnp, dimnames=list(Ped[,1], NULL))
    if (length(Inherit)==1)  Inherit <- rep(Inherit, nSnp)
    for (i in which(Ped$gen==0)) {
      SGeno.4[i, ] <- sapply(1:nSnp, function(l) sample.int(4, size=1,
                                prob=c(Q[l]^2, Q[l]*(1-Q[l]), Q[l]*(1-Q[l]), (1-Q[l])^2) *
                founderProp[Inherit[l], Ped$Sex[i], ]) )   # CHECK
    }

    # non-founders
    Gprob <- matrix(NA, nSnp, 4,
                    dimnames=list(1:nSnp, c("aa", "aA", "Aa", "AA")))
    Inherit <- factor(Inherit, levels = dimnames(INHA)[[1]])
    for (g in 1:nGen) {
      for (i in which(Ped$gen==g)) {
        if (Ped$damIDx[i] != 0 & Ped$sireIDx[i] != 0) {
          for (x in 1:4) {  #
            Gprob[,x] <- INHA[,Ped$Sex[i],x,,][cbind(Inherit, SGeno.4[Ped$damIDx[i],], SGeno.4[Ped$sireIDx[i],]) ]
          }
        } else {
          stop("Non-autosomal single parents not implemented yet!!")
        }
        SGeno.4[i,] <- apply(Gprob, 1, function(p) sample.int(4, size=1, prob=p))
      }
    }

    SGeno <- apply(SGeno.4, 1, function(v) c(0,1,1,2)[v])   # CHECK
  }


  #================================
  # genotyping errors & missing values:

  SGeno.actual <- SGeno
  SGeno <- MkGenoErrors(SGeno, CallRate, SnpError, ErFunc)

  #================================
  # Non-genotyped parents

  NotSampled <- which(apply(SGeno, 1, function(x) all(x==-9)))

  if (any(ParMis>0) & is.null(names((CallRate)))==1) {
    if (length(na.exclude(intersect(Ped[,2], Ped[,3]))) >0) {
      if (ParMis[1] != ParMis[2]) {
        stop("With hermaphrodites, 'ParMis' must be equal for dams & sires")
#      } else if (!quiet) {
#        message("detected hermaphrodites ... ")
      }
      IsParent <- which(Ped[,1] %in% Ped[,2] | Ped[,1] %in% Ped[,3])
      if (round(length(IsParent)*ParMis[1]) > 0) {
        NotSampled <- c(NotSampled,
                        sample(IsParent, round(length(IsParent)*ParMis[1]),
                               replace=FALSE) )
      }

    } else {

      for (p in 1:2) {
        if (ParMis[p]>0) {
          IsParent <- which(Ped[,1] %in% Ped[,p+1])
        }
        if (round(length(IsParent)*ParMis[p]) > 0) {
          NotSampled <- c(NotSampled,
                          sample(IsParent, round(length(IsParent)*ParMis[p]),
                                 replace=FALSE) )
        }
      }
    }
  }
  if (length(NotSampled)>0) {
    SGeno <- SGeno[-NotSampled, ]
    SGeno.actual <- SGeno.actual[-NotSampled, ]
    nInd.g <- nInd -length(NotSampled)
  } else {
    nInd.g <- nInd
  }

  #================================
  # output

  if (!is.na(OutFile)) {
    utils::write.table(SGeno, OutFile, quote=FALSE, col.names=FALSE)

  } else if (ReturnStats) {

    Params <- c(ParamsIN[c("nSnp", "ParMis", "MAF", "CallRate", "SnpError",
                           "Inherit", "ErrorFM", "OutFile", "InheritFile")],
                list(nInd = nrow(SGeno)))

    StatsOUT <- list(AF = apply(SGeno, 2, function(x) sum(x[x!=-9])/(2*sum(x!=-9))),
                     AF.act = apply(SGeno.actual, 2, function(x) sum(x[x!=-9])/(2*sum(x!=-9))),
                     SnpError = sapply(1:nSnp, function(l) {
                       sum(SGeno[,l] != SGeno.actual[,l] & SGeno[,l]!=-9)/nInd.g }),
                     SnpCallRate = apply(SGeno, 2, function(x) sum(x!=-9))/nInd.g,
                     IndivError = sapply(1:nInd.g, function(i) {
                       sum(SGeno[i,] != SGeno.actual[i,] & SGeno[i,]!=-9)/nSnp }),
                     IndivCallRate = apply(SGeno, 1, function(x) sum(x!=-9))/nSnp )

    return( list(ParamsIN = ParamsIN, StatsOUT = StatsOUT, SGeno = SGeno) )  # Ped = Ped,

  } else {
    return( SGeno )
  }
}


#=============================================================================
#' @title Simulate Genotyping Errors
#'
#' @description Generate errors and missing values in a (simulated) genotype
#'   matrix.
#'
#' @param SGeno  matrix with genotype data in Sequoia's format: 1 row per
#'   individual, 1 column per SNP, and genotypes coded as 0/1/2.
#' @param CallRate either a single number for the mean call rate (genotyping
#'   success), OR a vector with the call rate at each SNP, OR a named vector
#'   with the call rate for each individual. In the third case, ParMis is
#'   ignored, and individuals in the pedigree (as id or parent) not included in
#'   this vector are presumed non-genotyped.
#' @param SnpError  mean per-locus genotyping error rate across SNPs, and a
#'   beta-distribution will be used to simulate the number of missing cases per
#'   SNP, OR a vector with the genotyping error for each SNP.
#' @param ErrorFM  function taking the error rate (scalar) as argument and
#'   returning a 4x4 or 3x3 matrix with probabilities that actual genotype i
#'   (rows) is observed as genotype j (columns).
#' @param Error.shape first shape parameter (alpha) of beta-distribution of
#'   per-SNP error rates. A higher value results in a flatter distribution.
#' @param CallRate.shape as Error.shape, for per-SNP call rates.
#'
#' @return  The input genotype matrix, with some genotypes replaced, and some
#'   set to missing (-9).
#'
#' @examples
#' GenoM <- SimGeno(Ped = Ped_HSg5, nSnp = 100, ParMis = 0.2,
#'                  SnpError=0, CallRate=1)
#' GenoM.actual <- GenoM
#' LowQ <- sample.int(nrow(GenoM), 42)  # low-quality samples
#' GenoM[LowQ, ] <- MkGenoErrors(GenoM[LowQ, ], SnpError = 0.05)
#' GenoM[-LowQ, ] <- MkGenoErrors(GenoM[-LowQ, ], SnpError = 0.001)
#' ErrorCount <- sapply(1:nrow(GenoM), function(i) {
#'   sum(GenoM.actual[i,] != GenoM[i,] & GenoM[i,] != -9) } )
#' mean(ErrorCount[LowQ])
#' mean(ErrorCount[-LowQ])
#'
#' @export

MkGenoErrors <- function(SGeno,
                         CallRate = 0.99,
                         SnpError = 5e-4,
                         ErrorFM = function(E) {
                           matrix(c(1-E-(E/2)^2, E, (E/2)^2,
                                    E/2, 1-E, E/2,
                                    (E/2)^2, E, 1-E-(E/2)^2),
                                  3,3, byrow=TRUE) },
                         Error.shape=0.5,
                         CallRate.shape=1)
{
  nSnp <- ncol(SGeno)
  nInd <- nrow(SGeno)
#  if (! all(SGeno %in% c(0,1,2)) ) stop("SGeno may only contain 0, 1, or 2")

  #~~~~~~~~~
  if (any(SnpError >0)) {
    if (length(SnpError)==1) {
      El <- rbeta(nSnp, shape1=Error.shape, shape2=Error.shape*(1/SnpError -1))
    } else if (length(SnpError) == ncol(SGeno)) {
      El <- SnpError
    } else {
      stop("length of SnpError should equal 1 or number of SNPs")
    }

    shrinkET <- function(M) {
      N <- matrix(NA, 3,3)
      N[c(1,3), c(1,3)] <- M[c(1,4), c(1,4)]
      N[2, c(1,3)] <- M[2, c(1,4)]+M[3, c(1,4)]
      N[c(1,3), 2] <- M[c(1,4), 2] + M[c(1,4), 2]
      N[2, 2] <- sum(M[2:3, 2:3])
      return( N )
    }

    if (all(dim(ErrorFM(0.1))==3)) {
      RealToObs <- laply(El, ErrorFM)
    } else if (all(dim(ErrorFM(0.1))==4)) {
      RealToObs <- laply(El, function(e) shrinkET(ErrorFM(e)))
    } else {
      stop("ErrorFM(E) should return a 4x4 or 3x3 matrix")
    }

    # for (l in 1:nSnp) {
    #   SGeno[,l] <- sapply(SGeno[,l], function(x) sample.int(3, 1, prob=RealToObs[l,x+1,]) -1 )
    # } # rather slow; implemented in Fortran instead:
    SGeno <- DoErrors(SGeno, RealToObs)
  }

  #~~~~~~~~~
  if (any(CallRate <1)) {
    CRtype <- ifelse(length(CallRate)==1, "mean",
                ifelse(!is.null(names(CallRate)), "Indiv", "SNP"))
    MisX <- matrix(FALSE, nInd, nSnp)

    if (CRtype == "Indiv") {
      IndivCallRate <- setNames(CallRate[rownames(SGeno)], rownames(SGeno))
      IndivCallRate[is.na(IndivCallRate)] <- 0
      lmis <- round((1-IndivCallRate) *nSnp) # no. missing SNPs per indiv
      for (i in 1:nInd) {
        MisX[i, sample.int(nSnp, lmis[i])] <- TRUE   # ERROR
      }
    } else {
      if (CRtype == "mean") {
        imis <- round(rbeta(nSnp, CallRate.shape, CallRate.shape*(1/(1-CallRate) -1)) *nInd)
        # no. missing indiv per SNP
      } else if (CRtype == "SNP") {
        imis <- round((1-CallRate) *nInd)
      }
      for (l in 1:nSnp) {
        MisX[sample.int(nInd, imis[l]), l] <- TRUE
      }
    }
    SGeno[MisX] <- -9
  }

  #~~~~~~~~~
  return( SGeno )
}


#=============================================================================
#' @title Fortran Simulate Genotyping Errors
#'
#' @description Wrapper for Fortran function to simulate genotyping errors.
#'
#' @param SGeno matrix with genotype data, size nInd x nSnp.
#' @param RealToObs array with conditional probability of observing genotype i
#'   conditional on actual genotype j, size nSnp x 3 x 3.
#'
#' @return \code{SGeno} with errors.
#'
#' @useDynLib sequoia, .registration = TRUE
#'
#' @keywords internal

DoErrors <- function(SGeno, RealToObs) {
   dnames <- dimnames(SGeno)
   TMP <- .Fortran(mkerrors,
                   nind = as.integer(nrow(SGeno)),
                   nsnp = as.integer(ncol(SGeno)),
                   genofr = as.integer(SGeno),
                   eprobfr = as.double(RealToObs))
   return( matrix(TMP$genofr, nrow(SGeno), ncol(SGeno), dimnames=dnames) )
}


#=============================================================================

ReadSpecialInherit <- function(InheritFile, quiet) {
  inherit.L <- list()
  if (tools::file_ext(InheritFile) %in% c("xls", "xlsx")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      if (interactive() & !quiet) {
        ANS <- readline(prompt = paste("library 'openxlsx' not found. Install Y/N? "))
        if (!substr(ANS, 1, 1) %in% c("Y", "y")) stop()
      }
      utils::install.packages("openxlsx")
    }
    TypeNames <- openxlsx::getSheetNames(InheritFile)
    TmpL <- list()
    for (type in TypeNames) {
      TmpL[[type]] <-  openxlsx::read.xlsx(xlsxFile = InheritFile, sheet = type,
                                           colnames = TRUE, detectDates = FALSE)
    }

  } else {
    if (tools::file_ext(InheritFile)==".csv") {
      TmpDF <- utils::read.csv(InheritFile, stringsAsFactors=FALSE)
    } else {
      TmpDF <- utils::read.table(InheritFile, header=TRUE, stringsAsFactors=FALSE)
    }
    if (ncol(TmpDF)!=8)  stop("InheritFile must have 8 columns (possible read error)")
    TmpL <- plyr::dlply(TmpDF, "Mode", function(df) df[, -1])
  }

  for (y in seq_along(TmpL)) {
    if (is.null(TmpL[[y]]))  next
    INH <- array(0, dim = c(3,4,4,4),  # off sex, geno off - dam - sire
                 dimnames = c(list(c("fem","male", "unk")),
                              lapply(1:3, function(i) c("aa", "aA", "Aa", "AA"))))
    if (all(TmpL[[y]]$Sex == 3)) {
      for (i in 1:4) {
        INH[3,i,,] <- t(matrix(as.matrix(TmpL[[y]][,4:7])[,i], 4,4))
      }
      for (s in 1:2) {
        INH[s,,,] <- INH[3,,,]
      }

    } else if (all(TmpL[[y]]$Sex %in% 1:2)){
      for (s in 1:2) {
        if (!any(TmpL[[y]]$Sex == s)) next
        for (i in 1:4) {
          INH[s,i,,] <- t(matrix(as.matrix(TmpL[[y]][TmpL[[y]]$Sex==s, 4:7])[,i], 4,4))
        }
      }
      INH[3,,,] <- apply(INH[1:2,,,], c(2:4), mean)

    } else {
      stop("mix of known & unknown sex in INHERIT not implemented")
    }

    inherit.L[[TypeNames[y]]] <- INH
  }

  INHA <- plyr::laply(inherit.L, function(x) x)
  dimnames(INHA) <- c(list(names(inherit.L)), dimnames(INH))
  return( INHA )
}
