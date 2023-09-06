#=======================================================
#' @title Simulate Genotypes
#'
#' @description Simulate SNP genotype data from a pedigree, with optional
#'   missingness, genotyping errors, and non-genotyped parents.
#'
#' @details For founders, i.e. individuals with no
#'   known parents, genotypes are drawn according to the provided MAF and
#'   assuming Hardy-Weinberg equilibrium. Offspring genotypes are generated
#'   following Mendelian inheritance, assuming all loci are completely
#'   independent. Individuals with one known parent are allowed: at each locus,
#'   one allele is inherited from the known parent, and the other drawn from the
#'   genepool according to the provided MAF.
#'
#' @param Pedigree  dataframe, pedigree with the first three columns being id -
#'   dam - sire, additional columns are ignored.
#' @param nSnp  number of SNPs to simulate.
#' @param ParMis  single number or vector length two with proportion of parents
#'   with fully missing genotype. Ignored if CallRate is a named vector.
#' @param MAF  either a single number with minimum minor allele frequency, and
#'   allele frequencies will be sampled uniformly between this minimum and 0.5,
#'   OR a vector with minor allele frequency at each locus. In both cases, this
#'   is the MAF among pedigree founders, the MAF in the sample will deviate due
#'   to drift.
#' @param CallRate either a single number for the mean call rate (genotyping
#'   success), OR a vector with the call rate at each SNP, OR a named vector
#'   with the call rate for each individual. In the third case, ParMis is
#'   ignored, and individuals in the pedigree (as id or as parent) not included
#'   in this vector are presumed non-genotyped.
#' @param SnpError  either a single value which will be combined with
#'   \code{ErrorFM}, or a length 3 vector with probabilities (observed given
#'   actual) hom|other hom, het|hom, and hom|het; OR a vector or 3XnSnp matrix
#'   with the genotyping error rate(s) for each SNP.
#' @param ErrorFM  function taking the error rate (scalar) as argument and
#'   returning a 3x3 matrix with probabilities that actual genotype i (rows) is
#'   observed as genotype j (columns). See below for details
#' @param ReturnStats in addition to the genotype matrix, return the input
#'   parameters and mean & quantiles of MAF, error rate and call rates.
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
#' @section Genotyping errors:
#'   If \code{SnpError} is a length 3 vector, genotyping errors are generated
#'   following a length 3 vector with probabilities that 1) an actual homozygote
#'   is observed as the other homozygote, 2) an actual homozygote is observed as
#'   a heterozygote, and 3) an heterozygote is observed as an homozygote. The
#'   only assumption made is that the two alleles can be treated equally, i.e.
#'   observing actual allele $A$ as $a$ is as likely as observing actual $a$ as
#'   $A$.
#'
#'   If \code{SnpError} is a single value, by default this is interpreted as a
#'   locus-level error rate (rather than allele-level), and equals the
#'   probability that a homozygote is observed as heterozygote, and the
#'   probability that a heterozygote is observed as either homozygote (i.e., the
#'   probability that it is observed as AA = probability that observed as aa =
#'   \code{SnpError}/2). The probability that one homozygote is observed as the
#'   other is (\code{SnpError}/2\eqn{)^2}. How this single value is rendered
#'   into a 3x3 error matrix is fully flexible and specified via \code{ErrorFM};
#'   see \code{link{ErrToM}} for details.
#'
#'   The default values of \code{SnpError=5e-4} and \code{ErrorFM='version2.0'}
#'   correspond to the length 3 vector \code{c((5e-4/2)^2, 5e-4/2,
#'   5e-4*(1-5e-4/2))}.
#'
#'   A beta-distribution is used to simulate variation in the error rate between
#'   SNPs, the shape parameter of this distribution can be specified via
#'   \code{\link{MkGenoErrors}}. It is also possible to specify the error rate
#'   per SNP.
#'
#'
#' @section Call Rate:
#'   Variation in call rates across SNPs is assumed to follow a highly skewed
#'   (beta) distribution, with many SNPs having call rates close to 1, and a
#'   narrowing tail of lower call rates. The first shape parameter defaults to 1
#'   (but see \code{\link{MkGenoErrors}}), and the second shape parameter is
#'   defined via the mean as \code{CallRate}. For 99.9\% of SNPs to have a call
#'   rate of 0.8 (0.9; 0.95) or higher, use a mean call rate of 0.969 (0.985;
#'   0.993).
#'
#'   Variation in call rate between samples can be specified by providing a
#'   named vector to \code{CallRate}. Otherwise, variation in call rate and
#'   error rate between samples occurs only as side-effect of the random nature
#'   of which individuals are hit by per-SNP errors and drop-outs. Finer control
#'   is possible by first generating an error-free genotype matrix, and then
#'   calling \code{\link{MkGenoErrors}} directly on (subsets of) the matrix.
#'
#'
#' @seealso The wrapper \code{\link{EstConf}} for repeated simulation and
#'   pedigree reconstruction; \code{\link{MkGenoErrors}} for fine control over
#'   the distribution of genotyping errors in simulated data;
#'   \code{\link{EstEr}} to estimate genotyping error rate.
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
#' Geno_A <- SimGeno(Pedigree = Ped_griffin, nSnp=200, ParMis=c(0.1, 0.6),
#'                   MAF = 0.25, SnpError = 0.001)
#'
#' Geno_B <- SimGeno(Pedigree = Ped_HSg5, nSnp = 100, ParMis = 0.2,
#'                  SnpError = c(0.01, 0.04, 0.1))
#'
#' # genotype matrix with duplicated samples:
#' Dups_grif <- data.frame(ID1 = c('i006_2001_M', 'i021_2002_M', 'i064_2004_F'))
#' Dups_grif$ID2 <- paste0(Dups_grif$ID1, '_2')
#' Err <- c(0.01, 0.04, 0.1)
#' Geno_act <- SimGeno(Ped_griffin, nSnp=500, ParMis=0, CallRate=1, SnpError=0)
#' Geno_sim <- MkGenoErrors(Geno_act, SnpError=Err, CallRate=0.99)
#' Geno_dups <- MkGenoErrors(Geno_act[Dups_grif$ID1, ], SnpError=Err,
#'                           CallRate=0.99)
#' rownames(Geno_dups) <- Dups_grif$ID2
#' Geno_sim <- rbind(Geno_sim, Geno_dups)
#'
#'
#' @importFrom stats rbinom runif rbeta
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
					          quiet = FALSE)
{
  if (missing(Pedigree)) stop("please provide a pedigree to simulate from")

  # unavoidable partial matching when 'Err' is specified instead of 'SnpError'
  if (is.numeric(ErrorFM)) {
    SnpError <- ErrorFM
    ErrorFM <- "version2.0"   # problem doesn't occur when 'ErrorFM' specified; restore default
  }

  #================================
  # check input
  if(!(is.numeric(nSnp)) || nSnp<=0)  stop("nSnp should be a number greater than 0")

  params <- list('ParMis' = ParMis, 'MAF' = MAF,
                 'SnpError' = SnpError, 'CallRate' = CallRate)
  ValidLength <- list('ParMis' = c(1,2),
  'MAF' = c(1, nSnp),
  'SnpError' = c(1,3,nSnp,3*nSnp),
  'CallRate' = c(1, nSnp))

  for (p in seq_along(params)) {
    if (is.null(params[[p]]) | all(is.na(params[[p]]))) {
      stop("Please provide ", names(params)[p])
    } else if (!length(params[[p]]) %in% ValidLength[[p]]) {
      if (p < 4) {
        stop("Length ", names(params)[p], " should be one of ", ValidLength[[p]])
      } else { # callrate
        if (is.null(names(CallRate))) {
          if (length(CallRate) != nSnp)  stop("CallRate should be length 1, or length nSnp, or a named vector")
        } else {
          if (length(intersect(names(CallRate), Pedigree[,1]))==0) stop("names of CallRate vector do not match pedigree")
        }
      }
    } else if ( ! (is.numeric(params[[p]]) & all(params[[p]]>=0) & all(params[[p]]<=1)) ) {
      stop(names(params)[p], " must be a number between 0 and 1")
    }
  }

  if(length(ParMis)==1) ParMis <- rep(ParMis, 2)

  ParamsIN <- as.list(environment())

  #================================
  # minor allele frequencies (among founders)
  if (length(MAF)==1) {
    Q <- round(runif(nSnp, min=MAF, max=0.5),3)
  } else {
    Q <- as.numeric(MAF)
  }

  #================================
  # check & prep ===
  Ped <- sequoia::PedPolish(Pedigree, ZeroToNA=TRUE)
  nInd <- nrow(Ped)
  if (any(round(Q*nInd) %in% c(0,1)))  warning("some simulated SNPs have fixed alleles")


  #================================
  # simulate genotypes ===
  # founders: random draw of alleles under HWE
  # non-founders: following Mendelian inheritance

  # rownumber of dam & sire
  Ped$damIDx <- sapply(Ped[,2], function(x) ifelse(is.na(x), 0, which(Ped[,1]==x)))
  Ped$sireIDx <- sapply(Ped[,3], function(x) ifelse(is.na(x), 0, which(Ped[,1]==x)))

  #~~~~~~~~~~~~~~
  # divide pedigree into `generations': the parents of an individual must come
  # from earlier cohorts than itself, or from the founder population (gen 0)
  Ped$gen <- sequoia::getGenerations(Ped, StopIfInvalid=TRUE)
  nGen <- max(Ped$gen)


  #~~~~~~~~~~~~~~
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


  #================================
  # genotyping errors & missing values ===

  SGeno.actual <- SGeno
  # Error matrix: rows = actual, columns = observed
  SGeno <- MkGenoErrors(SGeno, CallRate, SnpError)


  #================================
  # Non-genotyped parents ===

  NotSampled <- which(apply(SGeno, 1, function(x) all(x<0)))

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

  if (ReturnStats) {

    Params <- c(ParamsIN[c("nSnp", "ParMis", "MAF", "CallRate", "SnpError")],
                list(nInd = nrow(SGeno)))

    StatsOUT <- list(AF = apply(SGeno, 2, function(x) sum(x[x>=0])/(2*sum(x>=0))),
                     AF.act = apply(SGeno.actual, 2, function(x) sum(x[x>=0])/(2*sum(x>=0))),
                     SnpError = sapply(1:nSnp, function(l) {
                       sum(SGeno[,l] != SGeno.actual[,l] & SGeno[,l]>=0)/nInd.g }),
                     SnpCallRate = apply(SGeno, 2, function(x) sum(x>=0))/nInd.g,
                     IndivError = sapply(1:nInd.g, function(i) {
                       sum(SGeno[i,] != SGeno.actual[i,] & SGeno[i,]>=0)/nSnp }),
                     IndivCallRate = apply(SGeno, 1, function(x) sum(x>=0))/nSnp )

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
#' @param SnpError  either a single value which will be combined with \code{ErrorFM},
#'    or a length 3 vector with  hom->other hom, hom->het, het-hom error rates;
#'     OR a vector or 3XnSnp matrix with the genotyping error rate(s) for each SNP.
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
#' @export

MkGenoErrors <- function(SGeno,
                         CallRate = 0.99,
                         SnpError = c((5e-4/2)^2, 5e-4/2, 5e-4*(1-5e-4/2)),
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

  #~~~~~~~~~
  if (any(SnpError >0)) {
    if (length(SnpError) %in% c(1, nSnp)) {
      if (length(SnpError)==1) {
        El <- rbeta(nSnp, shape1=Error.shape, shape2=Error.shape*(1/SnpError -1))
      } else if (length(SnpError) == ncol(SGeno)) {
        El <- SnpError
      }
      Act2Obs <- plyr::laply(El, ErrorFM)

    } else if (length(SnpError) %in% c(3, 3*nSnp)) {
      if (length(SnpError)==3) {
        Err_per_SNP <- sapply(SnpError, function(e) rbeta(nSnp, shape1=Error.shape, shape2=Error.shape*(1/e -1)))
      } else if (length(SnpError) == 3*ncol(SGeno)) {
        Err_per_SNP <- SnpError
      }
      Act2Obs <- plyr::aaply(Err_per_SNP, .margins=1, .fun=ErV2M)

    } else {
      stop("length of SnpError should equal 1, 3, nSnp, or be a matrix of 3 by number of SNPs")
    }

    # for (l in 1:nSnp) {
    #   SGeno[,l] <- sapply(SGeno[,l], function(x) sample.int(3, 1, prob=Act2Obs[l,x+1,]) -1 )
    # } # rather slow; implemented in Fortran instead:
    SGeno <- DoErrors(SGeno, Act2Obs)
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
#' @param Act2Obs array with conditional probability of observing genotype i
#'   conditional on actual genotype j, size nSnp x 3 x 3.
#'
#' @return \code{SGeno} with errors.
#'
#' @useDynLib sequoia, .registration = TRUE
#'
#' @keywords internal

DoErrors <- function(SGeno, Act2Obs) {
   dnames <- dimnames(SGeno)
   # Generate random numbers to determine which SNPs are erroneous (r < p)
   # random number generation by Fortran not allowed: F90 not always supported  
   # + may interfere with other pkgs
   
   randomV <- runif(n=nrow(SGeno)*ncol(SGeno), min=0, max=1)
   
   TMP <- .Fortran(mkerrors,
                   nind = as.integer(nrow(SGeno)),
                   nsnp = as.integer(ncol(SGeno)),
                   genofr = as.integer(SGeno),
                   eprobfr = as.double(Act2Obs),
                   randomv = as.double(randomV))
   return( matrix(TMP$genofr, nrow(SGeno), ncol(SGeno), dimnames=dnames) )
}


#=============================================================================

# vector -> matrix Pr(observed genotype (columns) | actual genotype (rows)):
ErV2M <- function(ErV, CallRate=NULL)
{
  if (is.null(CallRate)) {
    ErrM <- matrix(NA, 3,3, dimnames = list(act=0:2, obs=0:2))
  } else {
    ErrM <- matrix(NA, 3,4, dimnames = list(act=0:2, obs=-1:2))
  }
  ErrM['0', c('0','1','2')] <- c(1-ErV[1]-ErV[2],  ErV[2],  ErV[1])
  ErrM['1', c('0','1','2')] <- c(ErV[3],  1-2*ErV[3],  ErV[3])
  ErrM['2', c('0','1','2')] <- c(ErV[1],  ErV[2],  1-ErV[1]-ErV[2])

  if (!is.null(CallRate)) {
    ErrM[,c('0','1','2')] <- ErrM[,c('0','1','2')] * CallRate
    ErrM[,'-1'] <- 1 - CallRate
  }
  return( ErrM )
}
