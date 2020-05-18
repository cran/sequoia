# this is P(age|relationship) / P(Age), to be used in P(relationship|age) =
# P(age|relationship) * P(relationship) / P(age) where P(Age) is the emperical
# age distribution in the sample.

#======================================================================
#======================================================================

#' @title Age priors
#'
#' @description For various categories of pairwise relatives (R), calculate
#'   age-difference (A) based probability ratios \eqn{P(A|R) / P(A)} .
#'
#' @details The ratio \eqn{P(A|R) / P(A)} is the ratio between the observed
#'   counts of pairs with age difference A and relationship R (\eqn{N_{A,R}}),
#'   and the expected counts if age and relationship were independent
#'   (\eqn{N_{.,.}*p_A*p_R}).
#'
#'   During pedigree reconstruction, the ratios \eqn{P(A|R) / P(A)} calculated
#'   here are multiplied by the age-independent genetic-only \eqn{P(R|G)} to
#'   obtain a probability that the pair are relatives of type R conditional on
#'   both their age difference and their genotypes (i.e. using Bayes' theorem,
#'   \eqn{P(R|A, G) =P(A|R) / P(A) * P(R|G)}).
#'
#'   The age-difference prior is used for pairs of genotyped individuals, as
#'   well as for dummy individuals. This assumes that the propensity for a pair
#'   with a given age difference to both be sampled does not depend on their
#'   relationship, so that the ratio \eqn{P(A|R) / P(A)} does not differ between
#'   sampled and unsampled pairs.
#'
#' @section CAUTION: The small sample correction with \code{Smooth} and/or
#'   \code{Flatten} prevents errors in one dataset, but may introduce errors in
#'   another; a single solution that fits to the wide variety of life histories
#'   and datasets is impossible. Please do inspect the matrix, e.g. with
#'   \code{PlotAgePrior}.
#'
#' @section Single cohort: When no birth year information is given, or all
#'   individuals have the same birth year, it is assumed that a single cohort has
#'   been analysed and a matrix with 0's and 1's is returned. When
#'   \code{Discrete=FALSE}, avuncular pairs are assumed potentially present,
#'   while when \code{Discrete=TRUE} avuncular is not considered as a
#'   relationship possibility.
#'
#' @section Other time units: "Birth year" may be in any arbitrary time unit
#'   relevant to the species (day, month, decade), as long as parents are never
#'   born in the same time unit as their offspring, but always before their
#'   putative offspring (e.g. parent's BirthYear= 1 (or 2001) and offspring
#'   BirthYear=5 (or 2005)). Negative numbers and NA's are interpreted as unknown,
#'   and fractional numbers are not allowed.
#'
#' @section Maximum parental age: The number of rows in the output ageprior
#'   matrix equals the maximum parental age +1 (the first row is for age
#'   difference 0). The maximum parental age equals:
#'  \itemize{
#'  \item{}{the maximum age of parents if a pedigree is provided, or}
#'  \item{}{the (largest) value of \code{MaxAgeParent}, or}
#'  \item{}{1, if generations are discrete, or}
#'  \item{}{the maximum range of birth years in LifeHistData (including BY.min
#'  and BY.max, when provided)}}
#'  Exception is when \code{MaxAgeParent} is larger than the maximum age of
#'  parents in the provided skeleton pedigree, then \code{MaxAgeParent} is used.
#'  Thus, \code{MaxAgeParent} can be used when the birth year range in
#'  LifeHistData and/or the age distribution of assigned parents does not
#'  capture the absolutely maximum age of parents. Not adjusting this may hinder
#'  subsequent assignment of both dummy parents and grandparents.
#'
#' @param Pedigree dataframe with id - dam - sire in columns 1-3, and optional
#'   column with birth years. Other columns are ignored.
#' @param LifeHistData dataframe with 3 or 5 columns: id - sex (not used) -
#'   birth year (- BY.min - BY.max), with unknown birth years coded as negative
#'   numbers or NA. Column names are ignored, so the column order is important.
#'   "Birth year" may be in any arbitrary discrete time unit relevant to the
#'   species (day, month, decade), as long as parents are never born in the same
#'   time unit as their offspring. It may include individuals not in the
#'   pedigree, and not all individuals in the pedigree need to be in
#'   LifeHistData.
#' @param MaxAgeParent  maximum age of a parent, a single number (max across
#'   dams and sires) or a vector of length two (dams, sires). If NULL, it will
#'   be estimated from the data. If there are fewer than 20 parents of either
#'   sex assigned, \code{MaxAgeParent} is set to the maximum age difference in
#'   the birth year column of \code{Pedigree} or \code{LifeHistData}.
#' @param Discrete  Discrete generations? By default (NULL), discrete
#'   generations are assumed if all parent-offspring pairs have an age
#'   difference of 1, and all siblings an age difference of 0, and there are at
#'   least 20 pairs of each category (mother, father, maternal sibling, paternal
#'   sibling). Otherwise, overlapping generations are presumed. When
#'   \code{Discrete=TRUE} (explicitly or deduced), \code{Smooth} and
#'   \code{Flatten} are always automatically set to \code{FALSE}. Use
#'   \code{Discrete=FALSE} to enforce (potential for) overlapping generations.
#' @param Flatten To deal with small sample sizes for some or all relationships,
#'   calculate weighed average between the observed age difference distribution
#'   among relatives and a flat (0/1) distribution. When \code{Flatten=NULL}
#'   (the default) automatically set to TRUE when there are fewer than 20
#'   parents with known age of either sex assigned, or fewer than 20 maternal or
#'   paternal siblings with known age difference. Also advisable if the sampled
#'   relative pairs with known age difference are non-typical of the pedigree as
#'   a whole.
#' @param lambdaNW  Control weighing factors when \code{Flatten=TRUE}. Weights
#'   are calculated as \eqn{W(R) = 1 - exp(-lambdaNW * N(R))}, where \eqn{N(R)}
#'   is the number of pairs with relationship R for which the age difference is
#'   known. Large values (>0.2) put strong emphasis on the pedigree, small
#'   values (<0.0001) cause the pedigree to be ignored. Default results in
#'   \eqn{W=0.5} for \eqn{N=100}.
#' @param Smooth Smooth the tails of and any dips in the distribution? Sets dips
#'   (<10\% of average of neighbouring ages) to the average of the neighbouring
#'   ages, sets the age after the end (oldest observed age) to LR(end)/2, and
#'   assigns a small value (0.001) to the ages before the front (youngest
#'   observed age) and after the new end. Peaks are not smoothed out, as these
#'   are less likely to cause problems than dips, and are more likely to be
#'   genuine characteristics of the species. Is set to \code{FALSE} when
#'   generations do not overlap (\code{Discrete=TRUE}).
#' @param Plot  plot a heatmap of the results? Only when \code{Pedigree} is
#'   provided
#' @param Return  return only a matrix with the likelihood-ratio \eqn{P(A|R) /
#'   P(A)} (\code{"LR"}) or a list including also various intermediate statistics
#'   (\code{"all"}) ?
#' @param quiet suppress messages
#'
#' @return A matrix with the probability ratio of the age difference between two
#'   individuals conditional on them being a certain type of relative
#'   (\eqn{P(A|R)}) versus being a random draw from the sample (\eqn{P(A)}). For
#'   siblings and avuncular pairs, this is the absolute age difference.
#'
#'   The matrix has one row per age difference (0 - nAgeClasses) and five
#'   columns, one for each relationship type, with abbreviations:
#'   \item{M}{Mothers}
#'   \item{P}{Fathers}
#'   \item{FS}{Full siblings}
#'   \item{MS}{Maternal half-siblings}
#'   \item{PS}{Paternal half-siblings}
#'
#'   When \code{Return}='all', a list is returned with in addition to this
#'   matrix ('LR.RU.A') the following elements:
#'    \item{BirthYearRange}{vector length 2}
#'    \item{MaxAgeParent}{single number, estimated from the data or provided}
#'    \item{tblA.R}{matrix with the counts per age difference (0 -
#'   nAgeClasses) and the five relationship types as for 'LR.RU.A', plus a
#'   column 'X' with age differences across all pairs of individuals, including
#'   those in LifeHistData but not in Pedigree.}
#'    \item{Weights}{vector length 4, the weights used to flatten the
#'     distributions}
#'    \item{LR.RU.A.unweighed}{matrix with nAgeClasses+1 rows and 5 columns;
#'     LR.RU.A prior to flattening and smoothing}
#'    \item{Specs.AP}{the names of the input Pedigree and LifeHistData (or
#'    NULL), the 'effective' settings of Discrete, Smooth, and Flatten, and the
#'    value of lambdaNW}
#'
#' @seealso \code{\link{sequoia}} (and its argument \code{args.AP}),
#'   \code{\link{PlotAgePrior}} for visualisation. The age vignette gives
#'   further details, mathematical justification, and some examples.
#'
#' @examples
#' data(LH_HSg5, Ped_HSg5, package="sequoia")
#'
#' # no pedigree available:
#' MakeAgePrior(LifeHistData = LH_HSg5)
#' MakeAgePrior(LifeHistData = LH_HSg5, Discrete=TRUE)
#' MakeAgePrior(LifeHistData = LH_HSg5, MaxAgeParent = c(2,3))
#' \dontrun{
#' # with pedigree:
#' MakeAgePrior(Pedigree=Ped_HSg5[1:100,], LifeHistData = LH_HSg5)
#' MakeAgePrior(Ped_HSg5[1:100,], LH_HSg5, Discrete=FALSE)
#' # With 'Flatten', the value depens on the no. pairs per relationship:
#' MakeAgePrior(Ped_HSg5[1:100,], LH_HSg5, Flatten=TRUE)
#' AP.all <- MakeAgePrior(Ped_HSg5[1:200,], LH_HSg5, Flatten=TRUE)
#' AP.all$tblA.R
#' }
#'
#' @importFrom stats setNames
#'
#' @export

MakeAgePrior <- function(Pedigree = NULL,
                         LifeHistData = NULL,
                         MaxAgeParent = NULL,
                         Discrete = NULL,
                         Flatten = NULL,
                         lambdaNW = -log(0.5)/100,
                         Smooth = TRUE,
                         Plot = TRUE,
												 Return = "LR",
												 quiet = FALSE)
{
  call.AP <- as.list(match.call())
  # Input check
	if (is.null(LifeHistData)) {
    if (!is.null(Pedigree) && any(c("BY", "BirthYear", "birthyear", "Birthyear")
		 %in% colnames(Pedigree))) {
      LifeHistData <- data.frame(id = Pedigree[,1],
                                 Sex = 3,
                                 stringsAsFactors = FALSE)
			for (yy in c("BY", "BirthYear", "birthyear", "Birthyear")) {
				if (yy %in% colnames(Pedigree)) {
					LifeHistData$BirthYear <- Pedigree[, yy]
				}
			}
    } else {
      LifeHistData <- data.frame(id = NA,
                                 Sex = 3,
                                 BirthYear = NA)
    }
	}
  LifeHistData <- CheckLH(LifeHistData)
  LifeHistData$BirthYear[which(LifeHistData$BirthYear < 0)] <- NA
  if (any(colnames(LifeHistData) %in% c("BY.min", "BY.max"))) {
    LifeHistData$BY.min[which(LifeHistData$BY.min < 0)] <- NA
    LifeHistData$BY.max[which(LifeHistData$BY.max < 0)] <- NA
    BYrange <- suppressWarnings(range(unlist(LifeHistData[,
                        c("BirthYear", "BY.min", "BY.max")]), na.rm=TRUE))
  } else if (any(!is.na(LifeHistData$BirthYear))) {
    BYrange <- suppressWarnings(range(LifeHistData$BirthYear, na.rm=TRUE))
  } else {
    BYrange <- NA
  }

	if (is.null(lambdaNW) | is.na(lambdaNW))  lambdaNW <- -log(0.5)/100  # used for Flatten
	if (!Return %in% c("LR", "all"))  stop("Return must be 'LR' or 'all'")

	if (!is.null(Discrete) && is.na(Discrete))  Discrete <- NULL
	if (!(is.null(Discrete) || Discrete %in% c(TRUE, FALSE))) {
		stop("'Discrete' must be TRUE, FALSE, or NULL")
	}
	if (!is.null(Flatten) && is.na(Flatten))  Flatten <- NULL
	if (!(is.null(Flatten) || Flatten %in% c(TRUE, FALSE))) {
		stop("'Flatten' must be TRUE, FALSE, or NULL")
	}
	if (!Smooth %in% c(TRUE, FALSE))  stop("'Smooht' must be TRUE or FALSE")

	# ~~ maximum age of parents ~~
  if (is.null(MaxAgeParent) || all(is.na(MaxAgeParent))) {
    if (!is.null(Discrete) && Discrete) {
      MaxAgePO <- c(1,1)
    } else if (all(is.na(BYrange))) {
      stop("Must provide MaxAgeParent and/or birth years in LifeHistData when Discrete=FALSE")
    } else {
      MaxAgePO <- rep(max(1, min(diff(BYrange), 99)), 2)
    }
  } else {
    if (length(MaxAgeParent)==1) {
      MaxAgeParent <- rep(MaxAgeParent, 2)
    } else if (length(MaxAgeParent) > 2) {
      stop("MaxAgeParent must be NULL, a single number, or length 2 vector")
    }
    for (p in 1:2) {
      if (!is.na(MaxAgeParent[p]) && (MaxAgeParent[p]<0 || !is.wholenumber(MaxAgeParent[p]))) {
        stop("'MaxAgeParent' must be a positive whole number")
      }
    }
    if (!is.null(Discrete) && Discrete && (MaxAgeParent[1] != MaxAgeParent[2])) {
      stop("When Discrete=TRUE, MaxAgeParent must be identical for dams & sires (?)")
    } else {
      MaxAgePO <- rep(max(1, min(diff(BYrange), 99)), 2)
      for (p in 1:2) {
        if (!is.na(MaxAgeParent[p])) {
          MaxAgePO[p] <- MaxAgeParent[p]
        } # else: use default from BYrange
      }
    }
  }
  if (max(MaxAgePO) > 100)  stop("MaxAgePO must be smaller than 100; consider a different time unit")
	names(MaxAgePO) <- c("M", "P")

  MaxT = max(diff(BYrange), MaxAgePO, na.rm=TRUE)
  RR <- c("M", "P", "FS", "MS", "PS")  # relatedness categories considered

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ReturnDefault <- function(MaxP, MaxT,  Return, Disc=Discrete, quietR=quiet) {
    if (is.null(Disc))  Disc <- FALSE
		if (!quietR)  message("Ageprior: Default 0/1, ",
                          ifelse(Disc, "discrete","overlapping"),
                          " generations, MaxAgeParent=", MaxP[1],",",MaxP[2])
    LR.RU.A.default <- MkAPdefault(MaxP, MaxT, Disc)
    if (Plot)  PlotAgePrior(LR.RU.A.default)
  	if (Return == "all") {
	    Specs.AP <- list(Pedigree = call.AP[["Pedigree"]],
                 LifeHistData = call.AP[["LifeHistData"]],
                 Discrete = Discrete,
                 Flatten = Flatten,
                 lambdaNW = lambdaNW,
                 Smooth = Smooth)
      return( list(BirthYearRange = BYrange,
                   MaxAgeParent = MaxP,
                   tblA.R = NA,
                   RelativePairs.AgeKnown = NA,
                   Weights.AgeKnown = NA,
                   LR.RU.A.unweighed = NA,
                   LR.RU.A = LR.RU.A.default,
                   Specs.AP = Specs.AP) )
    } else {
      return( LR.RU.A.default )
    }
  }

	MkAPdefault <- function(MaxP, MaxT, Disc) {
		AP.default <- matrix(1, MaxT+1, 5,
                            dimnames=list(0:MaxT, RR))
		if(!is.null(Disc) && Disc) {  # always: MaxP[1]==MaxP[2]
			AP.default[,] <- 0
			AP.default[MaxP[1]+1, c("M","P")] <- 1
			AP.default[1, c("FS","MS","PS")] <- 1

  	} else {
			AP.default[,] <- 1
			AP.default[1, c("M", "P")] <- 0
			if (MaxP[1] < MaxT) {
				AP.default[((MaxP[1]+1):MaxT)+1, "M"] <- 0
				AP.default[(MaxP[1]:MaxT)+1, c("FS", "MS")] <- 0
			}
			if (MaxP[2] < MaxT) {
				AP.default[((MaxP[2]+1):MaxT)+1, "P"] <- 0
				AP.default[(MaxP[2]:MaxT)+1, c("FS", "PS")] <- 0
			}
  	}
		return( AP.default )
	}
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	# No pedigree or no parents w. known age: return default 0/1 prior
  if (is.null(Pedigree) | all(is.na(LifeHistData$BirthYear))) {
    OUT <- ReturnDefault(MaxAgePO, MaxT, Return)
    return( OUT )
  }
  Ped <- PedPolish(Pedigree, ZeroToNA=TRUE)
  Ped <- merge(Ped, LifeHistData[,c("ID","BirthYear")], by.x="id", by.y="ID", all.x=TRUE)
  if (sum(!is.na(Ped$dam) | !is.na(Ped$sire)) ==0) {
    OUT <- ReturnDefault(MaxAgePO, MaxT, Return)
    return( OUT )
  }

  # age difference per relatedness category
  Ped.R <- Ped[! (is.na(Ped$dam) & is.na(Ped$sire)), ]  # individuals with at least 1 parent (quicker)
  Ped.R <- merge(Ped.R, setNames(Ped[, c("id", "BirthYear")], c("dam", "BirthYear.dam")), all.x=TRUE)
  Ped.R <- merge(Ped.R, setNames(Ped[, c("id", "BirthYear")], c("sire", "BirthYear.sire")), all.x=TRUE)
  Ped.R$Age.dam <- with(Ped.R, BirthYear - BirthYear.dam)
  Ped.R$Age.sire <- with(Ped.R, BirthYear - BirthYear.sire)

  NAK.R <- setNames(rep(0,5), RR)  # Number Age Known given Relationship
  NAK.R["M"] <- sum(!is.na(Ped.R$Age.dam))
  NAK.R["P"] <- sum(!is.na(Ped.R$Age.sire))
	MinPairs.AgeKnown <- 20   # minimum no. mother-offspring or father-offspring pairs w known age diff
	for (p in 1:2) {
	  ParAge <- Ped.R[, paste0("Age.",c("dam","sire")[p])]
		if (NAK.R[p]>0) {
			if (MaxAgePO[p] == diff(BYrange) && any(!is.na(ParAge)) && NAK.R[p] >= MinPairs.AgeKnown) {
				MaxAgePO[p] <- max(ParAge, na.rm=TRUE)
			} else {
				MaxAgePO[p] <- max(ParAge, MaxAgePO[p], na.rm=TRUE)
			}
		}
	}
  if (!is.null(MaxAgeParent) & !quiet) {
		for (p in 1:2) {
			if (!is.na(MaxAgePO[p]) && !is.na(MaxAgeParent[p]) && MaxAgePO[p] > MaxAgeParent[p]) {
				warning("Some ", c("dams", "sires")[p], " older than MaxAgeParent, using new estimate")
			}
		}
  }

  BYrange <- c(min = min(c(Ped$BirthYear, LifeHistData$BirthYear), na.rm=TRUE),
               max = max(c(Ped$BirthYear, LifeHistData$BirthYear), na.rm=TRUE))
  MaxT = max(diff(BYrange), MaxAgePO, na.rm=TRUE)
	if (is.null(Discrete) || !Discrete) {
		if (Smooth)  MaxT <- MaxT +2	# space for a smooth-tail
	}

  if (all(NAK.R[c("M","P")]==0)) {
    OUT <- ReturnDefault(MaxAgePO, MaxT, Return)
    return( OUT )
  }

	# ~~~~~~~~~~~~~~~

  tblA.R <- matrix(NA, MaxT+1, length(RR)+1, dimnames=list(0:MaxT, c(RR, "X")))
  tblA.R[, "M"] <- table(factor(Ped.R$Age.dam, levels=0:MaxT))
  tblA.R[, "P"] <- table(factor(Ped.R$Age.sire, levels=0:MaxT))

  # Siblings
	RCM <- sapply(seq_along(Ped.R$id), GetRelCat, Ped.R[,c("id","dam","sire")], GenBack=1)  # slow.
  AgeDifM <- outer(Ped.R$BirthYear, Ped.R$BirthYear, "-")
  diag(AgeDifM) <- NA
  tblA.R[, "FS"] <- table(factor(abs(AgeDifM[RCM == "FS"]), levels=0:MaxT)) /2
  tblA.R[, "MS"] <- table(factor(abs(AgeDifM[RCM %in% c("FS", "MHS")]), levels=0:MaxT)) /2
  tblA.R[, "PS"] <- table(factor(abs(AgeDifM[RCM %in% c("FS", "PHS")]), levels=0:MaxT)) /2

  # Reference: age difference distribution across all pairs of individuals
  # *NOT* tblA.R[, "X"] <- table(factor(abs(AgeDifM), levels=0:MaxT))/2  : Ped.R is only indivs with at least 1 parent)
	Ped.LH <- merge(Ped[, colnames(Ped)!="BirthYear"], LifeHistData[,c("ID","BirthYear")],
	                by.x="id", by.y="ID",	all=TRUE)
	tbl.AgeDifs <- CountAgeDif(Ped.LH$BirthYear, BYrange)  # quicker than table(outer()), see below
  tblA.R[, "X"] <- tbl.AgeDifs[0:MaxT+1]
	tblA.R[is.na(tblA.R)] <- 0

  # note: some individuals in LifeHistData may be irrelevant to pedigree of SNPd indivs
	NAK.R <- apply(tblA.R, 2, sum, na.rm=TRUE)
	if (!is.null(Flatten) && !Flatten && (is.null(Discrete) || !Discrete)) {
		for (p in 1:2) {
			if (NAK.R[p] < 5) {
				warning("Fewer than 5 ", c("mother","father")[p],
											"-offspring pairs with known age difference, changing to Flatten=TRUE",
								immediate.=TRUE)
	      Flatten <- TRUE
			} else if (NAK.R[p] < MinPairs.AgeKnown) {
				message("Fewer than 20 ", c("mother","father")[p],
											"-offspring pairs with known age difference, please consider Flatten=TRUE")
			}
		}
	}

	#~~~~~~~~~~~~~~~~~
	if (((all(MaxAgePO == 1) & (is.null(MaxAgeParent) || all(MaxAgeParent==1 & !is.na(MaxAgeParent)))) |
	    (!is.null(MaxAgeParent) && !is.na(MaxAgeParent[1]) && all(MaxAgePO == MaxAgeParent[1]))) &
	    all(tblA.R[-1, c("FS", "MS", "PS")]==0)) {
		if (is.null(Discrete) & all(NAK.R[c("M","P", "MS", "PS")]>MinPairs.AgeKnown)) {
			Discrete <- TRUE
			if (Smooth) {
				MaxT <- MaxT -2  # remove space for smooth-tail again
				tblA.R <- tblA.R[1:(MaxT +1), ]
			}
		}
	} else {
		if (!is.null(Discrete) && Discrete) {
			stop("Discrete=TRUE, but some parents have age >1 or some siblings have age difference >0")
		}
	}
	if (is.null(Discrete))  Discrete <- FALSE
  if (Discrete)   Smooth <- FALSE
  if (is.null(Flatten)) {
		if (any(NAK.R[c("M","P", "MS", "PS")] < MinPairs.AgeKnown) & !Discrete) {
			Flatten <- TRUE
		} else {
			Flatten <- FALSE  # more precise AP speeds up computation
		}
  }


	#~~~~~~~~~~~~~~~~~

  PA.R <- sweep(tblA.R, 2, NAK.R, "/")
	PA.R[is.nan(PA.R)] <- 1.0

  FSuseHS <- FALSE
  if (NAK.R["FS"] / min(NAK.R[c("MS", "PS")]) < 0.5 & all(NAK.R[c("MS", "PS")]>MinPairs.AgeKnown)) {
    if (!Smooth | all(NAK.R[c("MS", "PS")] > 5*MinPairs.AgeKnown)) {
      FS.tmp <- PA.R[,"MS"] * PA.R[,"PS"]
      FS.tmp <- FS.tmp/sum(FS.tmp)
    } else {
      FS.tmp <- apply(PA.R[, c("MS", "PS")], 1, mean)
    }
    PA.R[,"FS"] <- (PA.R[,"FS"] + FS.tmp)/2
    FSuseHS <- TRUE
  }

#  PA.R <- sweep(PA.R, 2, colSums(PA.R, na.rm=TRUE), "/")  **superseded**
  for (r in RR) {
    if (any(!PA.R[,r] %in% c(0,1))) {
      PA.R[,r] <- PA.R[,r] / sum(PA.R[,r], na.rm=TRUE)
    }
  }
  PA.R[is.na(PA.R)] <- 0

  middle <- function(V,i) {
    if (i>1 & i<length(V)) {
      (V[i-1] + V[i+1])/2
    } else if (i==1) {
      V[i+1]
    } else if (i==length(V)) {
      V[i-1]
    }
  }
  ungap <- function(V, gap) {
    if(length(gap)==0)  return( V )
    if (gap[1]>1 & gap[length(gap)]<length(V)) {
      dYdX <- (V[gap[length(gap)]+1] - V[gap[1]-1]) / (length(gap)+1)
    } else {
      dYdX <- 0
    }
    for (i in seq_along(gap)) {
      if (gap[1]==1) {
        V[gap[i]] <- V[gap[length(gap)]+1] - i*dYdX
      } else {
        V[gap[i]] <- V[gap[1]-1] + i*dYdX
      }
    }
    V
  }

  # ~~~~~~~~~~~~~~~
  # probability ratio Related/Unrelated given Age difference
  if(any(!PA.R[, RR] %in% c(0,1))) {
    LR.RU.A.par <- sweep(PA.R[, RR], 1, PA.R[,"X"], "/")
  } else {
    LR.RU.A.par <- PA.R[, RR]
  }
  LR.RU.A.par[is.na(LR.RU.A.par)] <- 0   # if PA.R[,"X"]==0
  LR.RU.A <- LR.RU.A.par

  W.R <- 1 - exp(-lambdaNW * NAK.R[RR])
  if (Flatten) {
    # Weight of tblA.R versus flat prior, as function of sample size N
	  # default: lambdaNW = -log(0.5)/100 : <50% weight if N<100, and >50% if N>100
    if (FSuseHS) {
      W.R["FS"] <- 1-exp(-lambdaNW * mean(c(NAK.R["FS"], min(NAK.R[c("MS", "PS")]))))
    }
    W.R <- c(W.R, MGM=W.R[["M"]], PGF=W.R[["P"]], MGF=1/mean(1/W.R[c("P","M")]),
             UA=1/mean(1/W.R))

		LR.RU.A.default <- MkAPdefault(MaxAgePO, MaxT, Disc=Discrete)

    for (r in RR) {
      LR.RU.A[, r] <- W.R[r] * LR.RU.A.par[, r] + (1 - W.R[r]) * LR.RU.A.default[, r]
    }
  }

  MinP = 0.001
  if (Smooth) {
		#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SmoothAP <- function(V, MinP) {
      Front <- max(1, min(which(V > MinP)), na.rm=TRUE)
      End <- min(max(which(V > MinP)), length(V), na.rm=TRUE)
      lowMid <- rep(FALSE, length(V))
      if (End - Front >= 2) {
        for (i in c((Front+1) : (End-1))) {
          lowMid[i] <- ifelse(V[i] <= MinP, TRUE,
                          ifelse(V[i] / ((V[i-1] + V[i+1])/2) < 0.1, TRUE, FALSE))
        }
      }
      W <- V
      if (Front > 1 & W[Front] > 3*MinP)  					W[Front -1] <- MinP
      if (End < length(V))   												W[End +1] <- W[End]/2
			if ((End+1) < length(V) & W[End+1] > 3*MinP)  W[End+2] <- MinP
      for (x in 1:3) {  # in case neighbouring lowMid's
        for (i in seq_along(V)) {
          if (lowMid[i]) {
            W[i] <- (W[i-1] + W[i+1])/2
          }
        }
      }
      return( W )
    }
		#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    for (r in 1:ncol(LR.RU.A)) {
      LR.RU.A[, r] <- SmoothAP(LR.RU.A[, r], MinP = MinP)  # MinP = max(MinP, 1-W.R[r]))
    }
  }
  LR.RU.A[1, c("M", "P")] <- 0
  LR.RU.A[is.na(LR.RU.A)] <- 0
  LR.RU.A <- round(LR.RU.A,3)


  #~~~~~~~~~~~~~~~
  Specs.AP <- list(Pedigree = call.AP[["Pedigree"]],
                   LifeHistData = call.AP[["LifeHistData"]],
                   Discrete = Discrete,
                   Flatten = Flatten,
                   lambdaNW = lambdaNW,
                   Smooth = Smooth)

  OUT <- list(BirthYearRange = BYrange,
              MaxAgeParent = MaxAgePO,
              tblA.R = tblA.R,
              PA.R = PA.R,
              Weights = round(W.R,4),
              LR.RU.A.unweighed = round(LR.RU.A.par,3),
              LR.RU.A = LR.RU.A,
              Specs.AP = Specs.AP)
  if (Plot) {
		PlotAgePrior( AP = OUT[["LR.RU.A"]] )
  }
	if (!quiet)  message("Ageprior: Pedigree-based, ",
	                     ifelse(Discrete, "discrete ","overlapping "), "generations",
	                     ifelse(Flatten, ", flattened",""),
	                     ifelse(Smooth, ", smoothed",""),
	                     ", MaxAgeParent = ", MaxAgePO[1], ",", MaxAgePO[2])
  utils::flush.console()
  if (Return == "all") {
    return( OUT )
  } else {
    MaxA <- min(max(which(apply(LR.RU.A, 1, function(x) any(x>0)))) +1, nrow(LR.RU.A))
    return( LR.RU.A[1:MaxA, ])
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot age priors
#'
#' @description visualise the age-difference based prior probability ratios as a
#'   heatmap
#'
#' @param AP  matrix with age priors (P(A|R)/P(A)) with age differences in rows
#'  and relationships in columns; by default M: maternal parent (mother), P: paternal
#'  parent (father), FS: full siblings, MS: maternal siblings (full + half),
#'  PS: paternal siblings.
#' @param legend if TRUE, a new plotting window is started and
#'   \code{\link{layout}} is used to plot a legend next to the main plot. Set to
#'   FALSE if you want to add it as panel to an existing plot (e.g. with
#'   par(mfcol=c(2,2))).
#'
#' @return a heatmap
#'
#' @seealso \code{\link{MakeAgePrior}}, \code{\link{SummarySeq}}
#'
#' @importFrom graphics layout par image axis mtext abline
#'
#' @examples
#' \dontrun{
#' PlotAgePrior(SeqOUT$AgePriors)
#' }
#'
#' @export

PlotAgePrior <- function(AP = NULL, legend=TRUE)
{
  if (is.data.frame(AP))  AP <- as.matrix(AP)
  if (!is.matrix(AP))  stop("AP must be a matrix")
  if (any(AP < 0 | AP > 1000) | any(!is.double(AP)))
    stop("AP must be a numeric matrix with values >0 & < 1000")
  RR <- colnames(AP)
  if (ncol(AP)==15 && all(c("M", "P", "FS", "MS", "PS", "MGF", "PGM", "MFA", "PPA") %in% RR)) {  # re-arrange for clarity
    RR <- c("M", "P", "FS", "MS", "PS",
        "MGM", "MGF", "PGM", "PGF",
        "MFA", "PFA", "MMA", "MPA", "PMA", "PPA")
  }

  # trim off rows with only zero's. First row: A=0
  MaxA <- min(max(which(apply(AP, 1, function(x) any(x>0)))) +1, nrow(AP))
  if (!is.null(rownames(AP))) {
    AA <- as.numeric(rownames(AP))[1:MaxA]
  } else {
    AA <- 1:MaxA
  }

	stretch <- function(V, x=10) {
  	a <- list()
  	for (i in 1:(length(V)-1)) {
  	  a[[i]] <- seq(V[i], V[i+1], length.out=x+1)[1:x]
  	}
	  return(c(unlist(a), V[length(V)]))
	}
	flip <- function(M, dim=2)  apply(M, dim, function(x) rev(x))

	#~~~  heatmap  ~~~~~~~
	brks <- c(1:5, 10,20,50)
	brks <- c(0, 1/1001, rev(1/brks[-1]), brks) - 1e-5
	brks.f <- stretch(brks, x=10)
	mids <- (brks.f[-length(brks.f)] + brks.f[-1])/2
	cols <- c(1, grDevices::hcl.colors(89, "Light Grays"), grDevices::hcl.colors(70, "Greens", rev=TRUE))

	if (legend) {
	  oldpar <- par(no.readonly = TRUE)
  	ly <- layout(matrix(c(1,2), nrow=1), widths=c(.8, .2))
  	par(mai=c(.9,.9,.2,.1))
	}
	image(x=t(AP[1:MaxA, RR]), y=AA, xaxt="n", yaxt=ifelse(length(AA)<5, "n", "s"), las=1,
	      breaks = brks.f, col = cols,
	      ylab="Age difference (A)", cex.lab=1.1)
  if (length(AA)<5)  axis(side=2, at=AA, labels=AA, las=1, cex.lab=1.1)
	axis(side=1, at=seq(0,1,along.with=RR), line=-0.5, labels=RR,
	     cex.axis=ifelse(length(RR)==5, 1.1, 0.8),
	     las = ifelse(length(RR)==5, 1, 2),
	     tick=FALSE)
	if (any(AA<0)) abline(h=-0.4, col=2)
	mtext("Relationship (R)", side=1, line=2, cex=1.1)

	if (legend) {
  	par(mai=c(1.1,.3,0.5,.7))  # legend
  	E <- try(image(t(as.matrix(mids)), axes=FALSE, frame.plot=TRUE, breaks = brks.f, col = cols))
  	if (is.null(E)) {  # sometimes doesn't fit in plotting window
    	axis(side=4, at=seq(0, 1,length.out=length(brks)), labels=round(brks, 3), las=1, cex.axis=0.8)
    	axis(side=4, at=seq(0, 1,length.out=length(brks)), labels=FALSE, tck=0.2) # tcl=0.5)
      axis(side=2, at=seq(0, 1,length.out=length(brks)), labels=FALSE, tck=0.2)
    	mtext(("P(A|R)/P(A)"), side=3, line=0.5, cex=1)
  	} else {
  	  par(mai=rep(0,4))
      graphics::plot.new()  # so that can continue normally
  	}
  	par(oldpar)  # restore old par settings
	}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CountAgeDif <- function(BirthYear, BYrange = range(BirthYear)) {
	BYf <- factor(BirthYear, levels=c(BYrange[1]:BYrange[2]))
	BYM <- matrix(0, length(BYf), nlevels(BYf),
								dimnames = list(seq_along(BirthYear), levels(BYf)))
	BYM[cbind(seq_along(BYf), BYf)] <- 1
	BYM[rowSums(BYM)==0, ] <- NA

	AA <- outer(as.numeric(levels(BYf)), as.numeric(levels(BYf)), "-")
	BY.cnt <- apply(BYM, 2, sum, na.rm=TRUE)
	tmp <- outer(BY.cnt, BY.cnt, "*")
	A.cnt <- setNames(rep(0, max(AA)+1), 0:max(AA))
	for (a in 0:max(AA)) {
		A.cnt[a+1] <- sum(tmp[AA == a])
	}
	A.cnt["0"] <- A.cnt["0"] -sum(!is.na(BYM[,1]))   # self
	return( A.cnt )
}
