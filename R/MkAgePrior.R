# this is P(age|relationship) / P(Age), to be used in P(relationship|age) =
# P(age|relationship) * P(relationship) / P(age) where P(Age) is the emperical
# age distribution in the sample.

#======================================================================
#======================================================================

#' @title Age priors
#'
#' @description Calculate age-difference based prior probability ratios
#'   \eqn{P(A|R) / P(A)} for various categories of pairwise relatives.
#'
#' @details  Using Bayes' theorem, \deqn{P(Relationship | Age difference) =
#'   P(Age difference | Relationship) / P(Age difference) * P(Relationship)} or
#'   for short \eqn{P(R|A) =P(A|R) / P(A) * P(R)}. During pedigree
#'   reconstruction, the ratios \eqn{P(A|R) / P(A)} calculated here are
#'   multiplied by the age-independent genetic-only \eqn{P(Relationship)} to
#'   obtain a probability that the pair are relatives of type R conditional on
#'   both their age difference and their genotypes. This age-difference prior is
#'   used not only for pairs of genotyped individuals, but also between
#'   genotyped and dummy individuals and between pairs of dummy individuals.
#'   Therefore, \eqn{P(Age difference)} is taken in the broadest sense possible
#'   and calculated across all pairs of individuals in Ped and LifeHistData.
#'   When P(A) is 0 but The resulting distribution is by default flattened and
#'   smoothed to account for finite sample size.
#'
#' Sometimes no age information is available for one or more relative
#' categories, for example when all candidate fathers are unsampled and/or of
#' unknown age. To account for the amount of information available, if
#' Flatten=TRUE a weighed average of \eqn{P(A|R)/P(A)} as estimated from the
#' data and a completely flat distribution is returned. The weighing factor
#' \eqn{W(R)} for each relationship is a function of the number of relative
#' pairs with a known age difference \eqn{N(R)}, following the sigmoid curve
#' \eqn{W(R) = 1 - exp(-lambdaNW * N(R))}. By default,
#' \eqn{lambdaNW=-log(0.5)/100}, corresponding to a weight of <50\% if there are
#' <100 pairs, and >50\% if there are >100 pairs. See example below for a plot
#' of this curve.
#'
#' @section Discrete generations: When generations do not overlap, Flatten and
#'   Smooth should both be set to FALSE. This is done automatically when it is
#'   detected that all parents are 1 year older than their offspring, and
#'   siblings are always born in the same year.
#'
#' @section Single cohort: When no birth year information is given, or all
#'   individuals have the same birthyear, it is assumed that a single cohort has
#'   been analysed and a simple 2 by 9 matrix with 0's and 1's is returned. Note
#'   that by default it is assumed that no avuncular pairs exist within this
#'   single cohort; if they may exist, change LR.RU.A["0","UA"] from 0 to 1.
#'
#' @section Other time units: "Birth year" may be in any arbitrary time unit
#'   relevant to the species (day, month, decade), as long as parents are never
#'   born in the same time unit as their offspring. Time of
#'   birth/hatching/germination must be expressed as whole, positive numbers.
#'
#' @param Ped dataframe with pedigree, with id - dam - sire in columns 1-3, and
#'   optional column with birth years. Other columns are ignored.
#' @param LifeHistData dataframe with columns id - sex (not used) - birth year
#'   (unknown: negative number or NA). Column names are ignored, so the column
#'   order is important. "Birth year" may be in any arbitrary discrete time unit
#'   relevant to the species (day, month, decade), as long as parents are never
#'   born in the same time unit as their offspring.
#' @param MaxAgeParent  maximum age of a parent (max across dams and sires). If
#'   NULL, will be estimated from the data. If there are fewer than 100 parents
#'   of either sex assigned, \code{MaxAgeParent} is set to the maximum age
#'   difference in the birth year column of \code{Ped} or \code{LifeHistData}.
#' @param Flatten Calculate weighed average between the observed age difference
#'   distribution among the relative pairs with known age difference and a
#'   completely flat distribution. The weights are a function of the number of
#'   pairs and \code{lambdaNW}, see \code{Details}. Advisable if the relative
#'   pairs with known age difference non-typical of the pedigree as a whole or
#'   when their number is limited, and therefore automatically set to TRUE when
#'   there are fewer than 20 parents of either sex assigned. Not advisable when
#'   generations do not overlap.
#' @param lambdaNW  When \code{Flatten=TRUE}, weighing factors of the
#'   data-estimated age-difference distribution versus a flat distribution are
#'   calculated as \eqn{W(R) = 1 - exp(-lambdaNW * N(R))}, where \eqn{N(R)} is
#'   the number of pairs with relationship R for which the age difference is
#'   known. See Details below.
#' @param Smooth Smooth the tails of and any dips in the distribution? Sets dips
#'   (<10\% of average of neighbouring ages) to the average of the neigbouring
#'   ages, sets the age after the end (oldest observed age) to LR(end)/2, and
#'   assigns a small value (0.001) to the ages before the front (youngest
#'   observed age) and after the new end. Set to FALSE when generations do not
#'   overlap.
#' @param Plot  plot a 2-panel overview of the results?
#' @param Return  return only a matrix with the likelihood-ratio \eqn{P(A|R) /
#'   P(A)} ("LR") or a list including also various intermediate statistics
#'   ("all") ?
#' @param quiet suppress messages
#'
#' @return A matrix with the probability ratio of the age difference between two
#'   individuals conditional on them being a certain type of relative
#'   (\eqn{P(A|R)}) versus being a random draw from the sample (\eqn{P(A)}). For
#'   siblings and avuncular pairs, this is the absolute age difference.
#'
#'   The matrix has one row per age difference (0 - nAgeClasses) and nine
#'   columns, one for each relationship type, with abbreviations:
#'   \item{M}{Mothers} \item{P}{Fathers} \item{FS}{Full siblings}
#'   \item{MS}{Maternal half-siblings} \item{PS}{Paternal half-siblings}
#'   \item{MGM}{Maternal grandmother} \item{PGF}{Paternal grandfather}
#'   \item{MGF}{Maternal grandfathers and paternal grandmothers}
#'   \item{UA}{Avuncular (aunt/uncle - niece/nephew)}
#'
#'   When \code{Return}='all', a list is returned with in addition to this
#'   matrix ('LR.RU.A') the following elements:
#'   \item{BirthYearRange}{vector length 2}
#'   \item{MaxAgeParent}{single number, as estimated from the data or provided}
#'   \item{tblA.R}{matrix with the counts per age difference (0 - nAgeClasses)
#'   and the five relationship types as for 'LR.RU.A', plus a column 'X' with
#'   age differences across all pairs of individuals, including those in
#'   LifeHistData but not in Ped.}
#'   \item{Weights}{vector length 4, the weights
#'     used to flatten the distributions}
#'   \item{LR.RU.A.unweighed}{matrix with
#'     nAgeClasses+1 rows and 9 columns; LR.RU.A prior to flattening and
#'     smoothing}
#'
#' @seealso \code{\link{sequoia}, \link{PlotAgePrior}}
#'
#' @examples
#' data(LH_HSg5, Ped_HSg5, package="sequoia")
#' \dontrun{
#' MakeAgePrior(Ped_HSg5, LH_HSg5, Flatten=FALSE, Smooth=FALSE)
#' APlist <- MakeAgePrior(Ped_HSg5, LH_HSg5, Flatten=FALSE, Smooth=TRUE,
#'   Return="all") }
#'
#' # effect of lambdaNW on weights when Flatten=TRUE:
#' lambda1 <- -log(0.5)/100
#' lambda2 <- -log(1-.9)/100
#' curve(1-exp(-lambda1*x), lty=2, lwd=2, from=1, to=1000, log="x",
#'       xlab="N known age difference", ylab="W", las=1)
#' curve(1-exp(-lambda2*x), lty=2, lwd=2, col="blue", add=TRUE)
#' abline(h=c(0,.1,.25,.5,.75,.9,1), col="grey", lty=3)
#'
#' N = c(1,10,42,100,200, 432, 1e3, 1e4)
#' data.frame(N = N,
#'            W1 = round(1-exp(-lambda1*N), 3),
#'            W2 = round(1-exp(-lambda2*N), 3))
#'
#' @importFrom stats setNames
#'
#' @export

MakeAgePrior <- function(Ped = NULL,
                         LifeHistData = NULL,
                         MaxAgeParent = NULL,
                         Flatten = TRUE,
                         lambdaNW = -log(0.5)/100,
                         Smooth = TRUE,
                         Plot = TRUE,
												 Return = "LR",
												 quiet = FALSE)

{
  if (is.null(LifeHistData)) {
     LifeHistData <- data.frame(id = Ped[,1],
                                 Sex = 3,
                                 stringsAsFactors = FALSE)
    if (any(c("BY", "BirthYear", "birthyear", "Birthyear") %in% colnames(Ped))) {
      for (yy in c("BY", "BirthYear", "birthyear", "Birthyear")) {
        if (yy %in% colnames(Ped)) {
          LifeHistData$BY <- Ped[, yy]
        }
      }
    } else {
      LifeHistData$BY <- NA
    }
  } else {
    if (any(c("BY", "BirthYear", "birthyear", "Birthyear") %in% colnames(Ped))) {
      stop("Please provide *EITHER* LifeHistData *OR* Birthyear column in Ped")
    }
    LifeHistData <- CheckLH(LifeHistData)[,1:3]
  }
  if (all(is.na(LifeHistData$BY))) {
    if(!quiet)  message("No birth year information, assuming single cohort")
  }

  if (Flatten) {
    if (is.null(lambdaNW))  stop("Please provide lambdaNW")
	  if (is.na(lambdaNW))  lambdaNW <- -log(0.5)/100
  }
	if (!Return %in% c("LR", "all"))  stop("Invalid value of Return")

  RR <- c("M", "P", "FS", "MS", "PS")  # relatedness categories considered
  RRP <- c(RR, c("MGM", "PGF", "MGF", "UA"))

  names(LifeHistData) <- c("id", "Sex", "BY")
  LifeHistData$BY[which(LifeHistData$BY < 0)] <- NA
  LifeHistData$BY <- as.numeric(LifeHistData$BY)
  BYrange <- suppressWarnings(range(LifeHistData$BY, na.rm=TRUE))
  MaxT <- max(1, diff(BYrange), na.rm=TRUE)

  ReturnDefault <- function(MaxT, Return) {
   LR.RU.A.default <- matrix(1, MaxT+1, 9,
                            dimnames=list(0:MaxT, RRP))
  	LR.RU.A.default[1, c("M", "P")] <- 0
  	LR.RU.A.default[1:2, c("MGM", "PGF", "MGF")] <- 0
    if (MaxT==1) {
  		LR.RU.A.default["1", c("FS", "MS", "PS")] <- 0
  		LR.RU.A.default["0", "UA"] <- 0   # !!
    }
  	if (Return == "all") {
      return( list(BirthYearRange = BYrange,
                   MaxAgeParent = MaxT,
                   tblA.R = NA,
                   RelativePairs.AgeKnown = NA,
                   Weights.AgeKnown = NA,
                   LR.RU.A.unweighed = NA,
                   LR.RU.A = LR.RU.A.default) )
    } else {
      return( LR.RU.A.default )
    }
  }

  if (is.null(Ped) | all(is.na(LifeHistData$BY))) {
    OUT <- ReturnDefault(MaxT, Return)
    return( OUT )
  }

  PedIN <- Ped
  Ped <- setNames(Ped[,1:3], c("id", "dam", "sire"))
  for (i in 1:3) Ped[,i] <- as.character(Ped[,i])
  Ped <- AddParPed(Ped)
  Ped <- merge(Ped, LifeHistData[,c("id","BY")], all.x=TRUE)

  if (sum(!is.na(Ped$dam) | !is.na(Ped$sire)) ==0) {
    OUT <- ReturnDefault(MaxT, Return)
    return( OUT )
  }

  #~~~~~~~~~~~~~~~~~
  # age difference per relatedness category
  Ped.R <- Ped[! (is.na(Ped$dam) & is.na(Ped$sire)), ]  # individuals with at least 1 parent (quicker)
  NAK.R <- setNames(rep(0,5), RR)  # Number Age Known given Relationship

  # Parents
  Ped.R <- merge(Ped.R, setNames(Ped[, c("id", "BY")], c("dam", "BY.dam")), all.x=TRUE)
  Ped.R <- merge(Ped.R, setNames(Ped[, c("id", "BY")], c("sire", "BY.sire")), all.x=TRUE)
  Ped.R$Age.dam <- with(Ped.R, BY - BY.dam)
  Ped.R$Age.sire <- with(Ped.R, BY - BY.sire)

  NAK.R["M"] <- sum(!is.na(Ped.R$Age.dam))
  NAK.R["P"] <- sum(!is.na(Ped.R$Age.sire))
  if (any(NAK.R[c("M","P")]>0)) {
    MaxAgePO <- max(unlist(c(Ped.R$Age.dam, Ped.R$Age.sire)), na.rm=TRUE)
  } else {
    MaxAgePO <- max(1, diff(range(LifeHistData$BY, na.rm=TRUE)))
  }

  if (!is.null(MaxAgeParent)) {
    if (MaxAgePO > MaxAgeParent) {
      if(!quiet)  warning("Some parents older than MaxAgeParent, using new estimate")
    }
    MaxAgePO <- max(MaxAgeParent, MaxAgePO, na.rm=TRUE)
  } else if (any(NAK.R[c("M","P")] < 100)) {  # rather arbitrary threshold
    MaxAgePO <- max(1, MaxAgePO, diff(range(LifeHistData$BY, na.rm=TRUE)), na.rm=TRUE)
    for (p in 1:2) {
      if (NAK.R[p] < 100) {
        if(!quiet)  message(paste("Fewer than 100", c("mother","father")[p],
                      "-offspring pairs with known age difference, using MaxAgeParent =", MaxAgePO))
      }
    }
  }

  BYrange <- c(min = min(c(Ped$BY - MaxAgePO, LifeHistData$BY), na.rm=TRUE),
               max = max(c(Ped$BY, LifeHistData$BY), na.rm=TRUE))
  MaxT <- max(1, diff(BYrange), na.rm=TRUE)

  if (all(NAK.R[c("M","P")]==0)) {
    OUT <- ReturnDefault(MaxT, Return)
    return( OUT )
  }

  tblA.R <- matrix(NA, MaxT+1, length(RR)+1, dimnames=list(0:MaxT, c(RR, "X")))
  tblA.R[, "M"] <- table(factor(Ped.R$Age.dam, levels=0:MaxT))
  tblA.R[, "P"] <- table(factor(Ped.R$Age.sire, levels=0:MaxT))

  # Siblings
	RCM <- sapply(1:nrow(Ped.R), GetRelCat, Ped.R[,c("id","dam","sire")], GenBack=1)  # slow.
  AgeDifM <- outer(Ped.R$BY, Ped.R$BY, "-")
  diag(AgeDifM) <- NA
  tblA.R[, "FS"] <- table(factor(abs(AgeDifM[RCM == "FS"]), levels=0:MaxT)) /2
  tblA.R[, "MS"] <- table(factor(abs(AgeDifM[RCM %in% c("FS", "MHS")]), levels=0:MaxT)) /2
  tblA.R[, "PS"] <- table(factor(abs(AgeDifM[RCM %in% c("FS", "PHS")]), levels=0:MaxT)) /2
	tblA.R[is.na(tblA.R)] <- 0

  # Reference: age difference distribution across all pairs of individuals
  Ped.LH <- merge(Ped[, colnames(Ped)!="BY"], LifeHistData[,c("id","BY")], all=TRUE)
  tblA.R[, "X"] <- CountAgeDif(Ped.LH$BY, BYrange)  # quicker than table(outer()), see Utils.R
  # note: doubles! some dummies are same indiv as non-genotyped observed
  # note: some individuals in LifeHistData are irrelevant to pedigree of SNPd indivs


  #~~~~~~~~~~~~~~~~~
  NAK.R <- apply(tblA.R, 2, sum, na.rm=TRUE)
  PA.R <- sweep(tblA.R, 2, NAK.R, "/")
  PA.R[is.nan(PA.R)] <- 0

   FSuseHS <- FALSE
  if (any(NAK.R[c("MS", "PS")] > NAK.R["FS"]) & all(NAK.R[c("MS", "PS")]>20)) {
    if (NAK.R["FS"] / min(NAK.R[c("MS", "PS")]) < 0.05) {
      FS.tmp <- PA.R[,"MS"] * PA.R[,"PS"]
      FS.tmp <- FS.tmp/sum(FS.tmp)
      PA.R[,"FS"] <- (PA.R[,"FS"] + FS.tmp)/2
      FSuseHS <- TRUE
    }
  }

  PA.R <- cbind(PA.R, MGM=NA, PGF=NA, MGF=NA, UA=0)

  AgeSum <- outer(0:MaxT, 0:MaxT, "+")
  PRM <- list()
  PRM[["MGM"]] <- outer(PA.R[,"M"], PA.R[,"M"], "*")
  PRM[["PGF"]] <- outer(PA.R[,"P"], PA.R[,"P"], "*")
  PRM[["MGF"]] <- outer(PA.R[,"M"], PA.R[,"P"], "*")
  for (r in c("MGM", "PGF", "MGF")) {
    for (a in 0:MaxT) {
      PA.R[a+1, r] <- sum(PRM[[r]][AgeSum==a])
    }
  }

  # avuncular: average over FS/MS/PS of mum/dad, sib older/younger
  # ignored: P(FS) vs P(MS) vs P(PS)
  AUtmp <- array(dim=c(2,3,MaxT+1, MaxT+1),
                  dimnames=list(c("M", "P"), c("FS", "MS", "PS"), 0:MaxT, 0:MaxT))
  for (p in c("M", "P")) {
    for (s in c("FS", "MS", "PS")) {
      AUtmp[p,s,,] <- outer(PA.R[,p], PA.R[,s], "*")
    }
  }
  AgeDif.Y <- outer(0:MaxT, 0:MaxT, function(i,j)  abs(i - j))
  for (p in c("M", "P")) {
    for (s in c("FS", "MS", "PS")) {
      for (a in 0:MaxT) {
        PA.R[a+1, "UA"] <- PA.R[a+1, "UA"] + sum(AUtmp[p,s,,][AgeSum==a | AgeDif.Y==a] )
      }
    }
  }

  PA.R <- sweep(PA.R, 2, colSums(PA.R, na.rm=TRUE), "/")
  PA.R[is.na(PA.R)] <- 0


  # ~~~~~~~~~~~~~~~
  # probability ratio Related/Unrelated given Age difference
  if(any(!PA.R[, RRP] %in% c(0,1))) {
    LR.RU.A.par <- sweep(PA.R[, RRP], 1, PA.R[,"X"], "/")
  } else {
    LR.RU.A.par <- PA.R[, RRP]
  }
  LR.RU.A.par[is.na(LR.RU.A.par)] <- 0   # if PA.R[,"X"]==0
  if (any(!is.finite(LR.RU.A.par))) {   # cause: grandparent age, but none sampled
    LR.RU.A.par[!is.finite(LR.RU.A.par)] <- 5    # arbitrary "large value"
  }
  LR.RU.A <- LR.RU.A.par

  for (p in 1:2) {
    if(NAK.R[p]<20 & !Flatten) {
      if (!quiet)  warning(paste("Fewer than 20", c("mother","father")[p],
                    "-offspring pairs with known age difference, changing Flatten to TRUE"))
      Flatten <- TRUE
    }
  }
  if (Flatten) {
    if (all(NAK.R[c("M","P")]>100) & sum(NAK.R[c("MS","FS")])>100 & sum(NAK.R[c("PS","FS")])>100) {
      if (MaxAgePO==1 & all(tblA.R["0", c("FS","MS","PS")] == NAK.R[c("FS","MS","PS")])) {
        Flatten <- FALSE   # discrete generations.
        Smooth <- FALSE
      }
    }
  }

  W.R <- 1 - exp(-lambdaNW * NAK.R[RR])
  if (Flatten) {
    # Weight of tblA.R versus flat prior, as function of sample size N
	  # default: lambdaNW = -log(0.5)/100 : <50% weight if N<100, and >50% if N>100
    if (FSuseHS) {
      W.R["FS"] <- 1-exp(-lambdaNW * mean(c(NAK.R["FS"], min(NAK.R[c("MS", "PS")]))))
    }
    W.R <- c(W.R, MGM=W.R[["M"]], PGF=W.R[["P"]], MGF=1/mean(1/W.R[c("P","M")]),
             UA=1/mean(1/W.R))

		LR.RU.A.default <- matrix(1, MaxT+1, 9,
                              dimnames=list(0:MaxT, RRP))
		LR.RU.A.default[1, c("M", "P")] <- 0
		LR.RU.A.default[1:2, c("MGM", "PGF", "MGF")] <- 0

    for (r in RRP) {
      LR.RU.A[, r] <- W.R[r] * LR.RU.A.par[, r] + (1 - W.R[r]) * LR.RU.A.default[, r]
    }
  }

  MinP = 0.001
  if (Smooth) {
    SmoothAP <- function(V, MinP) {
      Front <- max(1, min(which(V > MinP)), na.rm=TRUE)
      End <- min(max(which(V > MinP)), length(V), na.rm=TRUE)
      lowMid <- rep(FALSE, length(V))
      if (End - Front > 2) {
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
        for (i in 1:length(V)) {
          if (lowMid[i]) {
            W[i] <- (W[i-1] + W[i+1])/2
          }
        }
      }
      return( W )
    }

    for (r in 1:ncol(LR.RU.A)) {
      LR.RU.A[, r] <- SmoothAP(LR.RU.A[, r], MinP = MinP)  # MinP = max(MinP, 1-W.R[r]))
    }
  }

  LR.RU.A[1, c("M", "P")] <- 0
  LR.RU.A[1:2, c("MGM", "PGF", "MGF")] <- 0
  LR.RU.A[is.na(LR.RU.A)] <- 0
  LR.RU.A <- round(LR.RU.A,3)

	OUT <- list(BirthYearRange = BYrange,
                 MaxAgeParent = MaxAgePO,
                 tblA.R = tblA.R,
								 PA.R = PA.R,
                 Weights = round(W.R,4),
                 LR.RU.A.unweighed = round(LR.RU.A.par,3),
                 LR.RU.A = LR.RU.A )

  if (Plot) {
		PlotAgePrior( APlist = OUT )
  }

  if (Return == "all") {
    return( OUT )
  } else {
    return( OUT[["LR.RU.A"]] )
  }

}


#' @title Plot age priors
#'
#' @description visualise the age-difference based prior probability ratios
#'
#' @param AP  matrix with age priors of dimension nAgeclasses by 9, as ' AgePriors' in Sequoia output
#' @param APlist output list from MakeAgePriors when called with Return='all', i.e. including
#'  various intermediate output.
#'
#' @seealso \code{\link{MakeAgePrior}}
#'
#' @importFrom graphics par plot lines points abline legend
#'
#' @export

PlotAgePrior <- function(AP = NULL,
                         APlist = NULL)
{

	if (is.list(AP)) {
		APlist <- AP
		AP <- NULL
	}
	if (!is.null(APlist)) {
		op <- par(mfcol=c(2,1), mai=c(.8, .9, .5, .2))
		LR.RU.A <- APlist[["LR.RU.A"]]
		PA.R <- APlist[["PA.R"]]
	} else {
		# use existing par(), to e.g. add panel to existing plots
		LR.RU.A <- AP
	}
	MaxT <- nrow(LR.RU.A)-1
	MinP <- 0.001
	RR <- c("M", "P", "FS", "MS", "PS")  # relatedness categories considered
  RP <- c("MGM", "PGF", "MGF", "UA")
	RRP <- c(RR, RP)

	COL.R <- c(M="red", P="blue", FS="mediumpurple", MS="lightpink2", PS="skyblue",
               MGM="seagreen2" , PGF="seagreen4" , MGF="seagreen3" , UA="peru")
	jit <- stats::setNames(c(seq(-.07,.07, length.out=5), rep(0,4)), RRP)

	xmax <- ifelse(!is.null(APlist), MaxT, max(which(rowSums(LR.RU.A)>0)))

	if (!is.null(APlist)) {
		plot(0:MaxT, PA.R[,"X"], type="l", lwd=2, col="grey", log="y", cex.lab=1.2,
				 ylim=c(2e-4, 1), ylab="P(A|R)", xlab="Age difference", las=1,
				 main="Age-difference distribution, per relationship")
		abline(v=seq(0,MaxT,10), h=10^(-4:0), col="grey", lty=3)
		points(0:MaxT, PA.R[,"X"], col="grey", pch=20, cex=.9)
		for (r in RR) {
			lines(0:MaxT+jit[r], PA.R[, r], col=COL.R[r], lwd=2)
			points(0:MaxT+jit[r], PA.R[, r], col=COL.R[r], pch=20, cex=.9)
		}
		for (r in RRP) {
			lines(0:MaxT, PA.R[,r], col=COL.R[r], lwd=1)
		}
	}

	plot(0:MaxT, LR.RU.A[,1], type="n", lwd=2, col="grey", las=1, log="y",
			 ylim=c(MinP, 1.1*max(LR.RU.A, na.rm=TRUE)), xlim=c(0,xmax),
			 ylab="P(R|A) / P(A)", xlab="Age difference", cex.lab=1.2,
			 main="Age-difference distribution; scaled, flattened* & smoothed*")
	legend("top", "*: if TRUE", bty="n")
	abline(v=seq(0,MaxT,10), h=10^(-4:0), col="grey", lty=3, xpd=FALSE)
	abline(h=1, col="grey", lwd=2, xpd=FALSE)
	for (r in RRP) {
		lines(0:MaxT+jit[r], LR.RU.A[, r], col=COL.R[r], lwd=ifelse(r %in% RR,2,1), lty=1)
		if (r %in% RR) points(0:MaxT+jit[r], LR.RU.A[, r], col=COL.R[r], pch=20, cex=.9)
	}
	legend("topright", c(RRP, "U"), title="Relationship", cex=0.8,
				 col=c(COL.R, "grey"), lwd=c(rep(2,5),rep(1,4),2),
				 inset=.02, bg=0, ncol=2)

	if (!is.null(APlist))   par(op)
}


