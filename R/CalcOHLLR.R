#=======================================================================
#' @title calculate OH and LLR
#'
#' @description Count opposite homozygous (OH) loci between parent-offspring
#'   pairs and Mendelian errors (ME) between parent-parent-offspring trios, and
#'   calculate the parental log-likelihood ratios (LLR).
#'
#' @details Any individuals in \code{Pedigree} that do not occur in \code{GenoM}
#'   are substituted by dummy individuals; a value of '0' in column
#'   'SNPd.id.dam' in the output means that either the focal individual or the
#'   dam was thus substituted, or both were. Use \code{\link{getAssignCat}} to
#'   distinguish between these cases.
#'
#'   The birth years in \code{LifeHistData} and the \code{AgePrior} are not used
#'   in the calculation and do not affect the value of the likelihoods for the
#'   various relationships, but they _are_ used during some filtering steps, and
#'   may therefore affect the likelihood _ratio_. The default
#'   (\code{AgePrior=FALSE}) assumes all age-relationship combinations are
#'   possible, which may mean that some additional alternatives are considered
#'   compared to the \code{\link{sequoia}} default, resulting in somewhat lower
#'   \code{LLR} values.
#'
#'   A negative LLR for A's parent B indicates either that B is not truely the
#'   parent of A, or that B's parents are incorrect. The latter may cause B's
#'   presumed true, unobserved genotype to greatly divert from its observed
#'   genotype, with downstream consequences for its offspring. In rare cases it
#'   may also be due to 'weird', non-implemented double or triple relationships
#'   between A and B.
#'
#'
#' @param Pedigree  dataframe with columns id-dam-sire. May include
#'   non-genotyped individuals, which will be treated as dummy individuals.
#' @param GenoM the genotype matrix
#' @param CalcLLR  calculate log-likelihood ratios for all assigned parents
#'   (genotyped + dummy/non-genotyped; parent vs. otherwise related). If
#'   \strong{\code{FALSE}}, only number of mismatching SNPs are counted (OH &
#'   ME), and parameters \code{LifeHistData}, \code{AgePrior}, \code{Err},
#'   \code{Tassign}, and \code{Complex} are \strong{ignored}. Note also that
#'   calculating likelihood ratios is much more time consuming than counting OH
#'   & ME.
#' @param LifeHistData   Dataframe with columns ID - Sex - BirthYear,
#'   and optionally columns BY.min and BY.max. If provided, used to delimit
#'   possible alternative relationships.
#' @param AgePrior logical (TRUE/FALSE) to estimate the ageprior from Pedigree
#'   and LifeHistData, or an agepriors matrix (see \code{\link{MakeAgePrior}}).
#'   Affects which alternative relationships are considered (only those where
#'   \eqn{P(A|R) / P(A) > 0}). When TRUE, \code{\link{MakeAgePrior}} is called
#'   using its default values.
#' @param Err estimated genotyping error rate, as a single number or 3x3 matrix.
#'   If a matrix, this should be the probability of observed genotype (columns)
#'   conditional on actual genotype (rows). Each row must therefore sum to 1.
#' @param ErrFlavour function that takes \code{Err} as input, and returns a 3x3
#'   matrix of observed (columns) conditional on actual (rows) genotypes, or
#'   choose from inbuilt ones as used in sequoia 'version2.0', 'version1.3', or
#'   'version1.1'. Ignored if \code{Err} is a matrix. See \code{\link{ErrToM}}.
#' @param Tassign  used to determine whether or not to consider some more exotic
#'   relationships when Complex="full".
#' @param Complex  determines which relationships are considered as
#'   alternatives. Either "full" (default), "simp" (simplified, ignores inbred
#'   relationships), or "mono" (monogamous).
#' @param GDX  call \code{\link{getAssignCat}} to classify individuals as
#'   genotyped (G), substitutable by a dummy (D) or neither (X).
#' @param quiet logical, suppress messages
#'
#' @return the \code{Pedigree} dataframe with additional columns:
#'  \item{LLRdam}{Log10-Likelihood Ratio (LLR) of this female being the mother,
#'  versus the next most likely relationship between the focal individual and
#'  this female (see Details for relationships considered)}
#'  \item{LLRsire}{idem, for male parent}
#'  \item{LLRpair}{LLR for the parental pair, versus the next most likely
#'   configuration between the three individuals (with one or neither parent
#'   assigned)}
#'  \item{OHdam}{Number of loci at which the offspring and mother are
#'    opposite homozygotes}
#'  \item{OHsire}{idem, for father}
#'  \item{MEpair}{Number of Mendelian errors between the offspring and the
#'    parent pair, includes OH as well as e.g. parents being opposing
#'    homozygotes, but the offspring not being a heterozygote. The offspring
#'    being OH with both parents is counted as 2 errors.}
#'  \item{SNPd.id.dam}{Number of SNPs scored (non-missing) for both individual
#'    and dam}
#'  \item{SNPd.id.sire}{Number of SNPs scored for both individual and sire}
#'  \item{id.cat}{Character denoting whether the focal individual is genotyped
#'  (G), substitutable by a dummy (D), or neither (X).}
#'  \item{dam.cat}{as id.cat, for dams. If id.cat and/or dam.cat is 'X', the
#'  dam cannot be assigned.}
#'  \item{sire.cat}{as dam.cat, for sires}
#'  \item{Sexx}{Sex in LifeHistData, or inferred Sex when assigned as part of
#'    parent-pair}
#'  \item{BY.est}{mode of birth year probability distribution}
#'  \item{BY.lo}{lower limit of 95\% highest density region of birth year
#'  probability distribution}
#'  \item{BY.hi}{higher limit}
#'
#' The columns 'LLRdam', 'LLRsire' and 'LLRpair' are only included when
#' \code{CalcLLR=TRUE}. The columns 'dam.cat' and 'sire.cat' are only included
#' when \code{GDX=TRUE}. The columns 'Sexx', 'BY.est', 'BY.lo' and 'BY.hi' are
#' only included when \code{LifeHistData} is provided, and at least one
#' genotyped individual has an unknown birthyear or unknown sex.
#'
#'
#' @seealso \code{\link{SummarySeq}} for visualisation of OH & LLR
#'   distributions; \code{\link{GenoConvert}} to read in various genotype data
#'   formats, \code{\link{CheckGeno}}; \code{\link{PedPolish}} to check and
#'   'polish' the pedigree; \code{\link{getAssignCat}} to find which id-parent
#'   pairs are both genotyped or can be substituted by dummy individuals;
#'   \code{\link{sequoia}} for pedigree reconstruction
#'
#' @examples
#' \dontrun{
#' # have a quick look for errors in an existing pedigree,
#' # without running pedigree reconstruction
#' PedA <- CalcOHLLR(Pedigree = MyOldPedigree, GenoM = MyNewGenotypes,
#'   CalcLLR=FALSE)
#'
#' # or run sequoia with CalcLLR=FALSE, and add OH + LLR later
#' SeqOUT <- sequoia(Genotypes, LifeHist, CalcLLR=FALSE)
#' PedA <- CalcOHLLR(Pedigree = SeqoUT$Pedigree[, 1:3], GenoM = Genotypes,
#'   LifeHistData = LIfeHist, AgePrior = TRUE, Complex = "full")
#'
#' # visualise
#' SummarySeq(PedA, Panels=c("LLR", "OH"))
#' }
#'
#' @useDynLib sequoia, .registration = TRUE
#
#' @export

CalcOHLLR <- function(Pedigree = NULL,
                      GenoM = NULL,
                      CalcLLR = TRUE,
                      LifeHistData = NULL,
                      AgePrior = FALSE,
                      Err = 1e-4,
                      ErrFlavour = "version2.0",
                      Tassign = 0.5,
                      Complex = "full",
                      GDX = TRUE,
                      quiet = FALSE) {

  on.exit(.Fortran(deallocall), add=TRUE)

  # input check & prep
  if (Err<0 || Err>1 || !is.double(Err)) stop("'Err' must be a number between 0 and 1")
  if (!is.double(Tassign) || Tassign<0)  stop("'Tassign' must be a positive number")
  if (!CalcLLR %in% c(TRUE, FALSE))  stop("'CalcLLR' must be TRUE or FALSE")

  if (Complex %in% c("full", "simp", "mono")) {
    Complx <- switch(Complex, full = 2, simp = 1, mono = 0)
  } else {
    stop("'Complex' must be 'full', 'simp', or 'mono'")
  }

  Excl <- CheckGeno(GenoM, quiet=TRUE)
  if (any(c("ExcludedSnps", "ExcludedSnps-mono") %in% names(Excl))) {
    stop("Please run 'CheckGeno' and filter SNPs")
  }
  if ("ExcludedIndiv" %in% names(Excl)) {
    stop("Please run 'CheckGeno' and remove individuals scored for <5% of SNPs")
  }
  gID <- rownames(GenoM)
  Ng <- nrow(GenoM)

  if (!is.null(LifeHistData)) {
    LhIN <- CheckLH(LifeHistData)
  } else {
    LhIN <- CheckLH(data.frame(ID = rownames(GenoM),
                               Sex = 3,
                               BirthYear = NA))
  }
  LHF <- orderLH(LhIN[LhIN$Sex %in% c(1:3), ], gID)

	if (is.null(Pedigree)) {
    stop("Must provide 'Pedigree'")
	}
  PedX <- PedPolish(Pedigree, gID, DropNonSNPd=FALSE)

  #~~~~~~~~~
  if (!AgePrior %in% c(TRUE, FALSE) && !is.matrix(AgePrior))
    stop("'AgePrior' must be TRUE or FALSE, or a matrix")
  if (is.matrix(AgePrior)) {
    if (any(AgePrior < 0 | AgePrior > 1000) | any(!is.double(AgePrior)))
      stop("AP must be a numeric matrix with values >0 & < 1000")
    if (!all(c("M", "P", "FS", "MS", "PS") %in% colnames(AgePrior)))
      stop("AgePrior must have columns M, P, FS, MS, PS")
    AP <- AgePrior[, c("M", "P", "FS", "MS", "PS")]
  } else if (AgePrior) {
    AP <- MakeAgePrior(PedX, LifeHistData, Plot = FALSE, quiet=quiet)
  } else {
    AP <- matrix(1, nrow=100, ncol=5,
                 dimnames=list(0:99, c("M", "P", "FS", "MS", "PS")))
    AP["0", c("M", "P")] <- 0
  }

  #~~~~~~~~~
  # turn non-genotyped parents (dummies + observed) into temporary dummies
  Dummifiable <- GetDummifiable(PedX, gID)
  n.FakeDummy <- sapply(Dummifiable, length)
  par <- c("dam", "sire")
  FakeDumPrfx <- c("FF0", "MM0")
  DPnc <- nchar(FakeDumPrfx)[1]
  Renamed <- list()
  PedR <- PedX
  for (k in 1:2) {
    if (n.FakeDummy[k] == 0) next
    Renamed <- c(Renamed,
                 list(data.frame(orig = Dummifiable[[k]],
                                 fkdum = paste0(FakeDumPrfx[k],
                                                formatC(1:min(n.FakeDummy[k],9999),
                                                        width=4, flag=0)),
                                 stringsAsFactors = FALSE)))
    PedR[, par[k]] <- Replace(PedR[, par[k]], Renamed[[k]][,"orig"], Renamed[[k]][, "fkdum"])
    PedR[, "id"] <- Replace(PedR[, "id"], Renamed[[k]][,"orig"], Renamed[[k]][, "fkdum"])
  }
  names(Renamed) <- par[n.FakeDummy > 0]

  PedNum <- IDToNum(PedR, gID, DumPrefix=FakeDumPrfx)
  PedPar <- PedNum[match(gID, PedR$id), 2:3]

  # Dummies  (some complicated folding, adapted from GetMaybeRel)
  DumParRF <- rep(0, 4*as.integer(Ng/2))
  if (any(n.FakeDummy > 0)) {
    SibshipGPs <- array(0, dim=c(2,max(n.FakeDummy),2),
                        dimnames=list(c("grandma", "granddad"),
                                      1:max(n.FakeDummy), c("mat", "pat")))
    for (k in 1:2) {
      if (n.FakeDummy[k] == 0) next
      PedDum.k <- PedNum[substr(PedR$id,1,DPnc)==FakeDumPrfx[k], ]
      SibshipGPs[,1:n.FakeDummy[k],k] <- t(PedDum.k[order(-PedDum.k[,"id"]), 2:3])
      for (s in 1:n.FakeDummy[k]) {
        for (g in 1:2) {
          x <- (k-1)*2*as.integer(Ng/2) + (s-1)*2 + g
          DumParRF[x] <- SibshipGPs[g,s,k]
        }
      }
    }
  }

	ErrM <- ErrToM(Err, flavour = ErrFlavour, Return = "matrix")  # performs checks if Err is already matrix
	sts <- SnpStats(GenoM, Plot=FALSE)
	MaxMismatchV <- CalcMaxMismatch(Err=ErrM, MAF=sts[,"AF"],
	                                qntl=0.999^(1/Ng))

  #~~~~~~~~~
  # call fortran
	SpecsInt <- c(nSnp = as.integer(ncol(GenoM)),     # 1
	              MaxMisDUP = as.integer(MaxMismatchV["DUP"]),  # 2
	              MaxMisOH = as.integer(ncol(GenoM)),    # 3
	              MaxMisME = as.integer(ncol(GenoM)),    # 4
	              SMax = as.integer(Ng/10),           # 5
	              Complx = as.integer(Complx),        # 6
	              quiet = as.integer(quiet),          # 7
	              nAgeCl = as.integer(nrow(AP)))      # 8
	SpecsDbl <- c(TF = -2.0,  # not used
	              TA = Tassign)  # minimal use, to select alternatives


  TMP <- .Fortran(getpedllr,
                  ng = as.integer(Ng),
                  specsint = as.integer(SpecsInt),
                  specsdbl = as.double(SpecsDbl),
									errv = as.double(ErrM),
                  calcllr = as.integer(CalcLLR),
                  genofr = as.integer(GenoM),
                  sexrf = as.integer(LHF$Sex),
                  byrf = as.integer(c(LHF$BirthYear, LHF$BY.min, LHF$BY.max)),
                  aprf = as.double(AP),
                  parentsrf = as.integer(PedPar),

                  ohrf = integer(3*Ng),
                  lrrf = double(3*Ng),
                  snpdboth = integer(2*Ng),
                  dumparrf = as.integer(DumParRF),
                  dumlrrf = double(3*Ng),
									dumbyrf = integer(3*Ng))

  # polish output
  TMP$lrrf[abs(TMP$lrrf - 999) < 0.1] <- NA
  TMP$lrrf <- round(TMP$lrrf, 2)
	TMP$dumlrrf[abs(TMP$dumlrrf - 999) < 0.1] <- NA
	TMP$dumlrrf <- round(TMP$dumlrrf, 2)
  TMP$ohrf[TMP$ohrf < 0] <- NA
  TMP$snpdboth[TMP$SnpdBoth < 0] <- NA

  NewCols <- data.frame(id = gID,
                        VtoM(TMP$lrrf, nc=3),
                        VtoM(TMP$ohrf, nc=3),
                        VtoM(TMP$snpdboth, nc=2))
  names(NewCols) <- c("id", "LLRdam", "LLRsire", "LLRpair",
                      "OHdam", "OHsire", "MEpair",
                      "SNPd.id.dam", "SNPd.id.sire")
	if (!is.null(LifeHistData) & any(LhIN$BirthYear < 0 | LhIN$Sex == 3)) {
		LhOUT <- data.frame(Sexx = TMP$sexrf,
												BY.est = TMP$byrf[1:Ng],
												BY.lo = TMP$byrf[1:Ng + Ng],
												BY.hi = TMP$byrf[1:Ng + 2*Ng],
												stringsAsFactors = FALSE)
		LhOUT$BY.est[LhOUT$BY.est < 0] <- NA
		NewCols <- cbind(NewCols, LhOUT)
	}

  if (any(n.FakeDummy > 0)) {
    NgOdd <- Ng%%2==1
    DumCols <- data.frame(id=c(Renamed[[1]][,"orig"],
                               Renamed[[2]][,"orig"]),
                          VtoM(TMP$dumlrrf, sum(n.FakeDummy), 3, NgOdd),
                          OHdam=NA, OHsire=NA, MEpair=NA)
    names(DumCols)[2:4] <- c("LLRdam", "LLRsire", "LLRpair")
    DumCols$SNPd.id.dam <- ifelse(!is.na(DumCols$LLRdam), 0, NA)
    DumCols$SNPd.id.sire <- ifelse(!is.na(DumCols$LLRsire), 0, NA)
    if (!is.null(LifeHistData) & any(LhIN$BirthYear < 0 | LhIN$Sex == 3)) {
      DumCols <- cbind(DumCols,
                       Sexx = rep(1:2, n.FakeDummy),
                       VtoM(TMP$dumbyrf, sum(n.FakeDummy),3, NgOdd))
      names(DumCols)[11:13] <- c("BY.est", "BY.lo", "BY.hi")
    }
    NewCols <- rbind(NewCols, DumCols)
  }
  if (!CalcLLR) {
    NewCols <- NewCols[, -c(2:4)]
  }
  PedX$rownum <- seq_along(PedX$id)
  PedX <- merge(PedX, NewCols, all.x=TRUE)  # merge changes order
	if (GDX) {
		PedA <- getAssignCat(PedX[, c("id", "dam", "sire")], gID)
		PedX <- merge(PedX,
		              PedA[, c("id", "id.cat", "dam.cat", "sire.cat")],
		              all.x=TRUE)
	}
  PedX <- PedX[order(PedX$rownum), ]
  PedX <- PedX[, names(PedX)!="rownum"]
  rownames(PedX) <- seq_along(PedX$id)

  ColumnOrder <- c("id", "dam", "sire",
                   "OHdam", "OHsire", "MEpair",
                   "LLRdam", "LLRsire", "LLRpair",
                   "SNPd.id.dam", "SNPd.id.sire",
                   "id.cat", "dam.cat", "sire.cat",
                   "BY.est", "BY.lo", "BY.hi")
  PedX <- PedX[, intersect(ColumnOrder, names(PedX))]

  return( PedX )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

