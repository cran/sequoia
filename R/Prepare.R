#============================================================================
#' @title preparations
#'
#' @description Check parameter values and save as named vector.
#'
#' @details  Please do not increasing the number of SNPs or individuals beyond
#'   the numbers present in the datasets, as this may cause R to crash.
#'
#' @param GenoM matrix with genotype data, size nInd x nSnp
#' @param LifeHistData Dataframe with 3 columns:
#'  \describe{
#'  \item{ID}{max. 30 characters long,}
#'  \item{Sex}{1 = females, 2 = males, other numbers = unknown,}
#'  \item{Birth Year}{(or hatching year) Zero and negative numbers are
#'    interpreted as missing values.}}
#' @param nAgeClasses Number of age classes (= no. rows in AgePriors)
#' @param MaxSibIter Maximum number of iterations of sibship clustering
#'   (up to 42).
#' @param Err Estimated genotyping error rate.
#' @param MaxMismatchV Maximum number of loci at which (1) a duplicate sample
#'   mismatches; (2) candidate parent and offspring are allowed to be opposite
#'   homozygotes; (3) parent-parent-offspring trios can have Mendelian errors.
#' @param Tfilter Threshold log-likelihood ratio between a proposed relationship
#'   versus unrelated, to select candidate relatives. Typically a negative
#'   value, related to the fact that unconditional likelihoods are calculated
#'   during the filtering steps. More negative values may decrease
#'   non-assignment, but will increase computational time.
#' @param Tassign Minimum log-likelihood ratio required for acceptance of
#'   proposed relationship, relative to next most likely relationship. Higher
#'   values result in more conservative assignments.
#' @param MaxSibshipSize  Maximum number of offspring for a single individual (a
#'   generous safety margin is advised).
#' @param DummyPrefix character vector of length 2 with prefixes for dummy dams
#'   (mothers) and sires (fathers); maximum 20 characters each.
#' @param Complexity Either "full" (default), "simp" (no explicit consideration
#'   of inbred relationships), "mono" (monogamous breeding system), or "herm"
#'   (hermaphrodites)
#' @param FindMaybeRel  Identify pairs of non-assigned likely relatives after
#'   pedigree reconstruction. Can be time-consuming in large datasets.
#' @param CalcLLR  Calculate log-likelihood ratios for all assignments. Can be
#'   time-consuming in large datasets.
#'
#' @return A 1-row dataframe with parameter values
#'
#' @keywords internal

SeqPrep <- function(GenoM = NULL,
                    LifeHistData = NULL,
                    nAgeClasses = 1,
                    MaxSibIter = 5,
                    Err = 0.0001,
                    MaxMismatchV = NULL,
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
  #if (any(make.names(IDs_geno) != IDs_geno)) {
  if (any(grepl(" ", IDs_geno))) {
    stop("GenoM rownames must not include spaces") #includes syntactically invalid names, see ?make.names")
  }
  nSnp <- ncol(GenoM)

  if (Err<0 || Err>1 || !is.double(Err)) stop("'Err' must be a number between 0 and 1")
  if (any(MaxMismatchV<0 | !is.wholenumber(MaxMismatchV)))  stop("'MaxMismatchV' but be positive whole numbers")
  if (!is.double(Tfilter))  stop("'Tfilter' must be a number (typically negative)")
  if (!is.double(Tassign) || Tassign<0)  stop("'Tassign' must be a positive number")
  if (MaxSibshipSize<0 || !is.wholenumber(MaxSibshipSize)) {
    stop("'MaxSibshipSize' must be a positive whole number")
  }
  if (!Complexity %in% c("full", "simp", "mono", "herm")) stop("'Complexity' must be 'full', 'simp', 'mono', or 'herm'")
  if (is.logical(UseAge))   UseAge <- ifelse(UseAge, "yes", "no")
  if (!UseAge %in% c("yes", "no", "extra")) stop("'UseAge' must be 'yes', 'no', or 'extra'")
  if (!FindMaybeRel %in% c(TRUE, FALSE)) stop("'FindMaybeRel' must be TRUE or FALSE")
  if (!CalcLLR %in% c(TRUE, FALSE)) stop("'CalcLLR' must be TRUE or FALSE")
  if (!is.wholenumber(MaxSibIter))  stop("'MaxSibIter' must be a whole number")


  DP <- sapply(DummyPrefix, function(x) ifelse(nchar(x)==1, paste0(x,"0"), x))
  if (any(substr(rownames(GenoM),1,nchar(DP[1]))==DP[1]) |
      any(substr(rownames(GenoM),1,nchar(DP[2]))==DP[2])) {
    stop("DummyPrefix must not occur in GenoM rownames")
  }

  if (!is.null(LifeHistData)) {
    nIndLH <- nrow(LifeHistData)
    if (nIndLH >1 & any(LifeHistData$BirthYear >= 0 & LifeHistData$ID %in% IDs_geno)) {
      nAgeClasses <- with(LifeHistData, diff(range(BirthYear[BirthYear >= 0 & ID %in% IDs_geno],
                                                   na.rm = TRUE))) + 1
    }
    if (nAgeClasses > 100) {
      nAgeClasses <- length(table(LifeHistData$BirthYear[LifeHistData$BirthYear >= 0]))
      if (nAgeClasses > 100) stop("Cannot handle >100 cohorts!")
    }
  }

  Specs <- data.frame(NumberIndivGenotyped = nIndG,
                      NumberSnps = nSnp,
                      GenotypingErrorRate = Err,
                      MaxMismatchDUP = MaxMismatchV["DUP"],
                      MaxMismatchOH = MaxMismatchV["OH"],
                      MaxMismatchME = MaxMismatchV["ME"],
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
                      SequoiaVersion = utils::packageVersion("sequoia"),
                      stringsAsFactors = FALSE)
  rownames(Specs) <- "Specs"
  Specs
}

#=======================================================================
#=======================================================================
#' @title check LifeHistData
#'
#' @description Check that the provided LifeHistData is in the correct format
#'
#' @param LifeHistData the dataframe with ID - Sex - Birth year, and
#'   optionally BY.min - BY.max
#'
#' @return a dataframe with LifeHistData formatted for use by the Fortran
#'    part of the program
#'
#' @keywords internal

CheckLH <- function(LifeHistData) {
  names(LifeHistData)[1:3] <- c("ID", "Sex", "BirthYear")
	if (ncol(LifeHistData)==3) {
		LifeHistData$BY.min <- NA
		LifeHistData$BY.max <- NA
	} else if (any(colnames(LifeHistData) %in% c("BY.min", "BY.max"))) {
		if (colnames(LifeHistData)[4] == "BY.min" & colnames(LifeHistData)[5] == "BY.max") {
			LifeHistData <- LifeHistData[, 1:5]
		} else {
		  if (any(colnames(LifeHistData) == "BY.min")) {
			  BYmin <- LifeHistData$BY.min
		  } else {
		    BYmin <- NA
		  }
		  if (any(colnames(LifeHistData) == "BY.min")) {
			  BYmax <- LifeHistData$BY.max
		  } else {
		    BYmax <- NA
		  }
			LifeHistData <- LifeHistData[, 1:3]
			LifeHistData$BY.min <- BYmin
			LifeHistData$BY.max <- BYmax
		}
	} else {
		LifeHistData <- LifeHistData[, 1:5]
		names(LifeHistData)[4:5] <- c("BY.min", "BY.max")
	}

  LifeHistData$ID <- as.character(LifeHistData$ID)
  LifeHistData <- unique(LifeHistData[!is.na(LifeHistData$ID), ])

  for (x in c("Sex", "BirthYear", "BY.min", "BY.max")) {
    IsInt <- check.integer(LifeHistData[,x])
    if (any(!IsInt, na.rm=TRUE)) {
      if (sum(!IsInt, na.rm=TRUE) > sum(!is.na(LifeHistData[,x]))/2) {
        stop("LifeHistData column ", x, " should be integers (whole numbers)")
      } else {
        warning("Converting all values in LifeHistData column ", x, " to integers",
                immediate. = TRUE)
      }
    }
    LifeHistData[, x] <- ifelse(IsInt, suppressWarnings(as.integer(as.character(LifeHistData[, x]))), NA)
		if (x=="Sex") {
			if ((sum(LifeHistData$Sex %in% c(1,2,4)) < nrow(LifeHistData)/10) &
				(sum(is.na(LifeHistData$Sex)) + sum(LifeHistData$Sex==3, na.rm=TRUE) < nrow(LifeHistData)/2)) {
			stop("LifeHistData column 2 should contain Sex coded as 1=female, 2=male, 3/NA=unknown, 4=hermaphrodite")
			}
			LifeHistData$Sex[is.na(LifeHistData$Sex)] <- 3
			LifeHistData$Sex[!LifeHistData$Sex %in% 1:4] <- 3
		} else {
			LifeHistData[is.na(LifeHistData[,x]), x] <- -999
			LifeHistData[LifeHistData[,x] < 0, x] <- -999
		}
  }
  LifeHistData <- unique(LifeHistData)
  if (any(duplicated(LifeHistData$ID))) {
    stop("Some IDs occur multiple times in LifeHistData, with different sex and/or birth year")
  }
  if (any(grepl(" ", LifeHistData$ID))) {
    stop("LifeHistData IDs must not include spaces")
  }

  return( LifeHistData )
}
