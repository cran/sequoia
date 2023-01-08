#=======================================================================
#' @title Check Genotype Matrix
#'
#' @description Check that the provided genotype matrix is in the correct
#'   format, and check for low call rate samples and SNPs.
#'
#' @section Thresholds: Appropriate call rate thresholds for SNPs and
#'   individuals depend on the total number of SNPs, distribution of call rates,
#'   genotyping errors, and the proportion of candidate parents that are SNPd
#'   (sibship clustering is more prone to false positives). Note that filtering
#'   first on SNP call rate tends to keep more individuals in.
#'
#' @param GenoM the genotype matrix.
#' @param quiet suppress messages.
#' @param Plot  display the plots of \code{\link{SnpStats}}.
#' @param Return  either 'GenoM' to return the cleaned-up genotype matrix, or
#'   'excl' to return a list with excluded SNPs and individuals (see Value).
#' @param DumPrefix length 2 vector, to check if these don't occur among
#'   genotyped individuals.
#'
#' @return If \code{Return='excl'} a list with, if any are found:
#'  \item{ExcludedSNPs}{SNPs scored for <10% of individuals; automatically
#'    excluded when running \code{\link{sequoia}}}
#'  \item{ExcludedSnps-mono}{monomorphic (fixed) SNPs; automatically excluded
#'    when running \code{\link{sequoia}}. This includes nearly-fixed SNPs with
#'    MAF \eqn{= 1/2N}. Column numbers are *after* removal of
#'    \code{ExcludedSNPs}, if any.}
#'  \item{ExcludedIndiv}{Individuals scored for <5% of SNPs; these cannot be
#'    reliably included during pedigree reconstruction. Individual call rate is
#'    calculated after removal of 'Excluded SNPs'}
#'  \item{Snps-LowCallRate}{SNPs scored for 10% -- 50% of individuals; strongly
#'    recommended to be filtered out}
#'  \item{Indiv-LowCallRate}{individuals scored for <50% of SNPs; strongly
#'    recommended to be filtered out}
#'
#' When \code{Return='excl'} the return is \code{\link{invisible}}, i.e. a check
#' is run and warnings or errors are always displayed, but nothing may be
#' returned.
#'
#' @seealso \code{\link{SnpStats}} to calculate SNP call rates;
#'   \code{\link{CalcOHLLR}} to count the number of SNPs scored in both focal
#'   individual and parent.
#'
#' @examples
#' GenoM <- SimGeno(Ped_HSg5, nSnp=400, CallRate = runif(400, 0.2, 0.8))
#' # the quick way:
#' GenoM.checked <- CheckGeno(GenoM, Return="GenoM")
#'
#' # the user supervised way:
#' Excl <- CheckGeno(GenoM, Return = "excl")
#' GenoM.orig <- GenoM   # make a 'backup' copy
#' if ("ExcludedSnps" %in% names(Excl))
#'   GenoM <- GenoM[, -Excl[["ExcludedSnps"]]]
#' if ("ExcludedSnps-mono" %in% names(Excl))
#'   GenoM <- GenoM[, -Excl[["ExcludedSnps-mono"]]]
#' if ("ExcludedIndiv" %in% names(Excl))
#'   GenoM <- GenoM[!rownames(GenoM) %in% Excl[["ExcludedIndiv"]], ]
#'
#' # warning about  SNPs scored for <50% of individuals ?
#' SnpCallRate <- apply(GenoM, MARGIN=2,
#'                      FUN = function(x) sum(x!=-9)) / nrow(GenoM)
#' hist(SnpCallRate, breaks=50, col="grey")
#' GenoM <- GenoM[, SnpCallRate > 0.6]
#'
#' # to be on the safe side, filter out low call rate individuals
#' IndivCallRate <- apply(GenoM, MARGIN=1,
#'                        FUN = function(x) sum(x!=-9)) / ncol(GenoM)
#' hist(IndivCallRate, breaks=50, col="grey")
#' GoodSamples <- rownames(GenoM)[ IndivCallRate > 0.8]
#'
#' @export

CheckGeno <- function(GenoM, quiet=FALSE, Plot=FALSE,
                      Return="GenoM", DumPrefix = c("F0", "M0"))
{

  # basic checks ----
  if (is.null(GenoM)) stop("please provide 'GenoM'")
  if (!is.matrix(GenoM)) stop("'GenoM' should be a numeric matrix")
	if (!all(GenoM %in% c(0,1,2,-9))) {
	  UniqueValues <- unique(c(GenoM))
	  InvalidValues <- UniqueValues[!UniqueValues %in% c(0,1,2,-9)]
	  stop(paste0("'GenoM' includes invalid values: ", "'",
	             paste(InvalidValues, collapse="', '"), "'"))
	}
  if (!Return %in% c("GenoM", "excl"))  stop("'Return' must be 'GenoM' or 'excl'")


  # check IDs ----
  if (is.null(rownames(GenoM))) stop("'GenoM' has no rownames, these should be the individual IDs")
  if (any(duplicated(rownames(GenoM))))  stop("'GenoM' has duplicate IDs. Please exclude or rename these samples,",
                                              " or run GenoConvert with UseFID=TRUE.")
  if (any(grepl(" ", rownames(GenoM))))  stop("GenoM rownames must not include spaces")

  DP <- sapply(DumPrefix, function(x) ifelse(nchar(x)==1, paste0(x,"0"), x))
  for (i in seq_along(DumPrefix)) {
    if (any(substr(rownames(GenoM),1,nchar(DP[i]))==DP[i])) {
      stop("DummyPrefix must not occur in GenoM rownames")
    }
  }

  warn <- function(...)
    warning(paste(strwrap(paste(...), prefix=" "), "\n"),    # strwrap() destroys whitespace
            call. = FALSE, immediate.=TRUE)


  # Exclude low quality SNPs ----
  Excl <- list()
  sstats <- SnpStats(GenoM, Plot = Plot)
  SNPs.TooMuchMissing <- sstats[,"Mis"] >= 0.9*nrow(GenoM)
  if (any(SNPs.TooMuchMissing)) {
    warn("There are ", sum(SNPs.TooMuchMissing)," SNPs scored for <10% of individuals, these will be excluded")
		GenoM <- GenoM[, !SNPs.TooMuchMissing]
		Excl[["ExcludedSnps"]] <- which(SNPs.TooMuchMissing)
  }
  SNPs.mono <- sstats[,"AF"] <= 1/(2*nrow(GenoM)) | sstats[,"AF"] >= 1 - 1/(2*nrow(GenoM))
  if (any(SNPs.mono)) {
    warn("There are ", sum(SNPs.mono)," monomorphic (fixed) SNPs, these will be excluded")
		GenoM <- GenoM[, !SNPs.mono]
		Excl[["ExcludedSnps-mono"]] <- which(SNPs.mono)
  }
  Nscored.2 <- apply(GenoM, 2, function(x) sum(x!=-9))
  if (any(Nscored.2 < nrow(GenoM)/2)) {
    warn(ifelse("ExcludedSnps" %in% names(Excl), "In addition, there", "There"),
                  "are ", sum(Nscored.2 < nrow(GenoM)/2)," SNPs scored for <50% of individuals,",
                  "it is strongly advised to exclude those")
    Excl[["Snps-LowCallRate"]] <- which(Nscored.2 < nrow(GenoM)/2)
  }


  # Exclude low quality individuals ----
  Lscored <- apply(GenoM, 1, function(x) sum(x!=-9))
  if (any(Lscored < ncol(GenoM)/20)) {
    warn("\n *********** \n There are ", sum(Lscored < ncol(GenoM)/20),
                  " individuals scored for <5% of SNPs,",
                  "these WILL BE IGNORED \n ***********")
		Excl[["ExcludedIndiv"]] <- rownames(GenoM)[which(Lscored < ncol(GenoM)/20)]
		GenoM <- GenoM[Lscored >= ncol(GenoM)/20, ]
  }

  if (any(Lscored < ncol(GenoM)/2)) {
    warn("There are ", sum(Lscored < ncol(GenoM)/2)," individuals scored for <50% of SNPs,",
                  "it is strongly advised to exclude those")
		Excl[["Indiv-LowCallRate"]] <- rownames(GenoM)[which(Lscored < ncol(GenoM)/2)]
#		GenoM <- GenoM[Lscored >= ncol(GenoM)/2, ]
  }


  if (!quiet) {
    MSG <- paste("There are ", nrow(GenoM), " individuals and ", ncol(GenoM), " SNPs.\n")
    if (any(c("ExcludedSnps", "ExcludedSnps-mono", "ExcludedIndiv") %in% names(Excl))) {
      message("After exclusion, ", MSG)
    } else if (length(Excl)>0) {
      message(MSG)
    } else {
      message("Genotype matrix looks OK! ", MSG)
    }
  }

  if (Return == "GenoM") {
    return( GenoM )
  } else {
    invisible( Excl )
  }
}

