#============================================================================
#============================================================================
#' @title preparations
#'
#' @description Check parameter values and save as named vector.
#'
#' @details  Please do not increasing the number of SNPs or individuals beyond
#'  the numbers present in the datasets, as this may cause R to crash.
#'
#' @param GenoM matrix with genotype data, size nInd x nSnp
#' @param LifeHistData Dataframe with 3 columns:
#'  \itemize{
#'  \item{ID: }{max. 30 characters long,}
#'  \item{Sex: }{1 = females, 2 = males, other numbers = unkown,}
#'  \item{Birth Year: }{(or hatching year) Zero and negative numbers are
#'    interpreted as missing values.}}
#' @param nAgeClasses Number of age classes (= no. rows in AgePriors)
#' @param MaxSibIter Maximum number of iterations of sibship clustering
#'   (up to 42).
#' @param Err Estimated genotyping error rate.
#' @param MaxMismatch Maximum number of loci at which candidate parent and
#'  offspring are allowed to be opposite homozygotes, or be excluded.
#' @param Tfilter Threshold log-likelihood ratio between a proposed
#'   relationship versus unrelated, to select candidate relatives. Typically a
#'   negative value, related to the fact that unconditional likelihoods are
#'   calculated during the filtering steps. More negative values may decrease
#'   non-assignment, but will increase computational time.
#' @param Tassign Minimum log-likelihood ratio required for acceptance of
#'  proposed relationship, relative to next most likely relationship. Higher
#'  values result in more conservative assignments.
#' @param MaxSibshipSize  Maximum number of offspring for a single individual
#'  (a generous safety margin is advised).
#' @param DummyPrefix character vector of length 2 with prefixes for dummy
#'   dams (mothers) and sires (fathers); maximum 20 characters each.
#' @param Complexity Either "full" (default), "simp" (no explicit consideration
#'   of inbred relationships), "mono" (monogamous breeding system), or "herm"
#'   (hermaphrodites)
#' @param FindMaybeRel  Identify pairs of non-assigned likely relatives after
#'   pedigree reconstruction. Can be time-consuming in large datasets.
#' @param CalcLLR  Calculate log-likelihood ratios for all assignments. Can be
#'  time-consuming in large datasets.
#' @return A 1-row dataframe with parameter values
#'
#' @keywords internal

SeqPrep <- function(GenoM = NULL,
                     LifeHistData = NULL,
                     nAgeClasses = 1,
                     MaxSibIter = 5,
                     Err = 0.0001,
                     MaxMismatch = 3,
                     Tfilter = -2.0,
                     Tassign = 0.5,
            				 MaxSibshipSize = 100,
            				 DummyPrefix = c("F", "M"),
            				 Complexity = "full",
            				 UseAge = "yes",
            				 FindMaybeRel = TRUE,
            				 CalcLLR = TRUE)

{
  nIndG <- nrow(GenoM)
  IDs_geno <- rownames(GenoM)
  nSnp <- ncol(GenoM)

  if (Err<0 || Err>1 || !is.double(Err)) stop("invalid value for 'Err'")
  if (MaxMismatch<0 || !is.wholenumber(MaxMismatch))  stop("invalid value for 'MaxMismatch'")
  if (!is.double(Tfilter))  stop("invalid value for 'Tfilter'")
  if (!is.double(Tassign) || Tassign<0)  stop("invalid value for 'Tassign'")
  if (MaxSibshipSize<0 || !is.wholenumber(MaxSibshipSize)) {
    stop("invalid value for MaxSibshipSize")
  }
  if (!Complexity %in% c("full", "simp", "mono", "herm")) stop("invalid value for 'Complexity'")
  if (!UseAge %in% c("yes", "no", "extra")) stop("invalid value for 'UseAge'")
  if (!FindMaybeRel %in% c(TRUE, FALSE)) stop("invalid value for 'FindMaybeRel'")
  if (!CalcLLR %in% c(TRUE, FALSE)) stop("invalid value for 'CalcLLR'")
  if (!is.wholenumber(MaxSibIter))  stop("invalid value for 'MaxSibIter'")

  if (!is.null(LifeHistData)) {
    nIndLH <- nrow(LifeHistData)
    if (nIndLH >1 & any(LifeHistData$BY >= 0) & any(LifeHistData$ID %in% IDs_geno)) {
      nAgeClasses <- with(LifeHistData, diff(range(BY[BY >= 0 & ID %in% IDs_geno],
                                                 na.rm = TRUE))) + 1
    }
    if (nAgeClasses > 100) {
      nAgeClasses <- length(table(LifeHistData$BY[LifeHistData$BY >= 0]))
       if (nAgeClasses > 100) stop("Cannot handle >100 cohorts!")
    }
  }

  Specs <- data.frame(NumberIndivGenotyped = nIndG,
             NumberSnps = nSnp,
             GenotypingErrorRate = Err,
             MaxMismatch = MaxMismatch,
             Tfilter = Tfilter,
             Tassign = Tassign,
      			 nAgeClasses = nAgeClasses,
             MaxSibshipSize = MaxSibshipSize,
      			 MaxSibIter = MaxSibIter,
             DummyPrefixFemale = DummyPrefix[1],
             DummyPrefixMale = DummyPrefix[2],
             Complexity = Complexity,
      			 UseAge = UseAge,
      			 FindMaybeRel = FindMaybeRel,
      			 CalcLLR = CalcLLR,
      			 stringsAsFactors = FALSE)
  Specs
}



#=======================================================================
#' @title check GenoM
#'
#' @description Check that the provided genotype matrix is in the correct format
#'
#' @param GenoM the genotype matrix
#'
#' @return updates GenoM internally to remove low call rate individuals (<.5)
#'   and low call rate SNPs (<0.1), if necessary
#'
#' @keywords internal

CheckGeno <- function(GenoM) {
	if (is.null(GenoM)) stop("please provide 'GenoM'")
  if (!is.matrix(GenoM)) stop("'GenoM' should be a numeric matrix")
	if (!all(GenoM %in% c(0,1,2,-9))) stop("'GenoM' not in <seq> format, please use GenoConvert")

  if (is.null(rownames(GenoM))) stop("'GenoM' has no rownames, these should be the individual IDs")
	if (any(duplicated(rownames(GenoM))))  stop("'GenoM' has duplicate IDs. Please exclude or rename these samples, or run GenoConvert with UseFID=TRUE.")

  Excl <- list()
	Lscored <- apply(GenoM, 1, function(x) sum(x!=-9))
  if (any(Lscored < ncol(GenoM)/2)) {
    warning(paste("There are ", sum(Lscored < ncol(GenoM)/2)," individuals scored for <50% of SNPs, these will be excluded"),
            immediate.=TRUE)
		GenoM <<- GenoM[Lscored >= ncol(GenoM)/2, ]
		Excl[["ExcludedInd"]] <- rownames(GenoM)[which(Lscored < ncol(GenoM)/2)]
  }

  Nscored <- apply(GenoM, 2, function(x) sum(x!=-9))
  if (any(Nscored < nrow(GenoM)/10)) {
    warning(paste("There are ", sum(Nscored < nrow(GenoM)/10)," SNPs scored for <10% of individuals, these will be excluded"),
            immediate.=TRUE)
		GenoM <<- GenoM[, Nscored >= ncol(GenoM)/10]
		Excl[["ExcludedSnps"]] <- which(Nscored < nrow(GenoM)/10)
  }
  Excl
}


