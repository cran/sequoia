#' @title Pedigree Reconstruction
#'
#' @description Perform pedigree reconstruction based on SNP data, including
#'   parentage assignment and sibship clustering.
#'
#' @details Dummy parents of sibships are denoted by F0001, F0002, ... (mothers)
#'   and M0001, M0002, ... (fathers), are appended to the bottom of the
#'   pedigree, and may have been assigned real or dummy parents themselves (i.e.
#'   sibship-grandparents). A dummy parent is not assigned to singletons.
#'
#'   For each pair of candidate relatives, the likelihoods are calculated
#'   of them being parent-offspring (PO), full siblings (FS), half siblings
#'   (HS), grandparent-grandoffspring (GG), full avuncular (niece/nephew -
#'   aunt/uncle; FA), half avuncular/great-grandparental/cousins (HA), or
#'   unrelated (U). Assignments are made if the likelihood ratio (LLR) between
#'   the focal relationship and the most likely alternative exceed the threshold
#'   Tassign.
#'
#'   Further explanation of the various options and interpretation of the output
#'   is provided in the vignette.
#'
#' @param GenoM  numeric matrix with genotype data: One row per individual, and
#'   one column per SNP, coded as 0, 1, 2 or -9 (missing). Use
#'   \code{\link{GenoConvert}} to convert genotype files created in PLINK using
#'   --recodeA or in Colony's 2-column format to this format.
#' @param LifeHistData Dataframe with 3 columns (optionally 5):
#'  \describe{
#'  \item{ID}{max. 30 characters long,}
#'  \item{Sex}{1 = females, 2 = males, other = unknown, except 4 =
#'    hermaphrodite,}
#' \item{BirthYear }{birth or hatching year, integer, with missing values as NA
#'   or any negative value.}
#' \item{BY.min}{minimum birth year, only used if BirthYear is missing}
#' \item{BY.max}{maximum birth year, only used if BirthYear is missing} }
#' If the species has multiple generations per year, use an integer coding such
#' that the candidate parents' `Birth year' is at least one smaller than their
#' putative offspring's. Column names are ignored, so ensure column order is ID
#' - sex - birth year (- BY.min - BY.max).
#' @param SeqList list with output from a previous run, containing the elements
#'   `Specs', `AgePriors' and/or `PedigreePar', as described below, to be used
#'   in the current run. If \code{SeqList$Specs} is provided, all other input
#'   parameter values except \code{MaxSibIter} are ignored.
#' @param MaxSibIter number of iterations of sibship clustering, including
#'   assignment of grandparents to sibships and avuncular relationships between
#'   sibships. Set to 0 to not (yet) perform this step, which is by far the most
#'   time consuming and may take several hours for large datasets. Clustering
#'   continues until convergence or until MaxSibIter is reached.
#' @param Err estimated genotyping error rate, as a single number or 3x3 matrix.
#'   If a matrix, this should be the probability of observed genotype (columns)
#'   conditional on actual genotype (rows). Each row must therefore sum to 1.
#' @param ErrFlavour function that takes \code{Err} as input, and returns a 3x3
#'   matrix of observed (columns) conditional on actual (rows) genotypes, or
#'   choose from inbuilt ones as used in sequoia 'version2.0', 'version1.3', or
#'   'version1.1'. Ignored if \code{Err} is a matrix. See \code{\link{ErrToM}}.
#' @param MaxMismatch DEPRECATED AND IGNORED. Now calculated using
#'   \code{\link{CalcMaxMismatch}}.
#' @param Tfilter threshold log10-likelihood ratio (LLR) between a proposed
#'   relationship versus unrelated, to select candidate relatives. Typically a
#'   negative value, related to the fact that unconditional likelihoods are
#'   calculated during the filtering steps. More negative values may decrease
#'   non-assignment, but will increase computational time.
#' @param Tassign minimum LLR required for acceptance of proposed relationship,
#'   relative to next most likely relationship. Higher values result in more
#'   conservative assignments. Must be zero or positive.
#' @param MaxSibshipSize  maximum number of offspring for a single individual (a
#'   generous safety margin is advised).
#' @param DummyPrefix character vector of length 2 with prefixes for dummy dams
#'   (mothers) and sires (fathers); maximum 20 characters each.
#' @param Complex  either "full" (default), "simp" (simplified, no explicit
#'   consideration of inbred relationships), "mono" (monogamous) or "herm"
#'   (hermaphrodites, otherwise like "full").
#' @param UseAge  either "yes" (default), "no", or "extra" (additional rounds
#'   with extra reliance on ageprior, may boost assignments but increased risk
#'   of erroneous assignments); used during full reconstruction only.
#' @param args.AP list with arguments to be passed on to
#'   \code{\link{MakeAgePrior}}.
#' @param FindMaybeRel DEPRECATED, advised to run \code{\link{GetMaybeRel}}
#'   separately. TRUE/FALSE to identify pairs of non-assigned likely relatives
#'   after pedigree reconstruction. Can be time-consuming in large datasets.
#' @param CalcLLR  calculate log-likelihood ratios for all assigned parents
#'   (genotyped + dummy; parent vs. otherwise related). Time-consuming in large
#'   datasets. Can be done separately with \code{\link{CalcOHLLR}}.
#' @param quiet suppress messages: TRUE/FALSE/"verbose".
#' @param Plot display plots from \code{\link{SnpStats}, \link{MakeAgePrior}},
#'   and \code{\link{SummarySeq}}. Defaults (NULL) to TRUE when quiet=FALSE or
#'   "verbose", and FALSE when quiet=TRUE. If you get errors an error 'figure
#'   margins too large', enlarge the plotting area (drag with mouse). 'invalid
#'   graphics state' error can be dealt with by clearing the plotting area with
#'   dev.off().
#'
#' @return A list with some or all of the following components:
#' \item{AgePriors}{Matrix with age-difference based probability ratios for
#'   each relationship, used for full pedigree reconstruction; see
#'   \code{\link{MakeAgePrior}} for details. When running only parentage
#'   assignment (MaxSibIter=0) the returned AgePriors has been updated to
#'   incorporate the information of the assigned parents, and is ready for use
#'   during full pedigree reconstruction.}
#' \item{DummyIDs}{Dataframe with pedigree for dummy individuals, as well as
#' their sex, estimated birth year (point estimate, upper and lower bound of
#' 95\% confidence interval), number of offspring, and offspring IDs (genotyped
#' offspring only).}
#' \item{DupGenotype}{Dataframe, duplicated genotypes (with different IDs,
#'  duplicate IDs are not allowed). The specified number of maximum mismatches
#'   is used here too. Note that this dataframe may include pairs of closely
#'   related individuals, and monozygotic twins.}
#' \item{DupLifeHistID}{Dataframe, row numbers of duplicated IDs in life
#'   history dataframe. For convenience only, but may signal a problem. The
#'   first entry is used.}
#' \item{ErrM}{Error matrix; probability of observed genotype (columns)
#'   conditional on actual genotype (rows)}
#' \item{ExcludedInd}{Individuals in GenoM which were excluded because of a
#'   too low genotyping success rate (<50\%).}
#' \item{ExcludedSNPs}{Column numbers of SNPs in GenoM which were excluded
#'   because of a too low genotyping success rate (<10\%).}
#' \item{LifeHist}{Provided dataframe with sex and birth year data.}
#' \item{LifeHistPar}{LifeHist with additional columns 'Sexx' (inferred Sex when
#' assigned as part of parent-pair), 'BY.est' (mode of birth year probability
#' distribution), 'BY.lo' (lower limit of 95\% highest density region), 'BY.hi'
#' (higher limit), inferred after parentage assignment. 'BY.est' is NA when the
#' probability distribution is flat between 'BY.lo' and 'BY.hi'.}
#' \item{LifeHistSib}{as LifeHistPar, but estimated after full pedigree
#' reconstruction}
#' \item{MaybeParent}{Dataframe with pairs of individuals who are more likely
#'   parent-offspring than unrelated, but which could not be phased due to
#'   unknown age difference or sex, or for whom LLR did not pass Tassign.}
#' \item{MaybeRel}{Dataframe with pairs of individuals who are more likely
#'   to be first or second degree relatives than unrelated, but which could not
#'   be assigned.}
#' \item{MaybeTrio}{Dataframe with non-assigned parent-parent-offspring trios
#' (both parents are of unknown sex), with similar columns as the pedigree}
#' \item{NoLH}{Vector, IDs in genotype data for which no life history data is
#'  provided.}
#' \item{Pedigree}{Dataframe with assigned genotyped and dummy parents from
#'   Sibship step; entries for dummy individuals are added at the bottom.}
#' \item{PedigreePar}{Dataframe with assigned parents from Parentage step.}
#' \item{Specs}{Named vector with parameter values.}
#' \item{TotLikParents}{Numeric vector, Total likelihood of the genotype data
#'   at initiation and after each iteration during Parentage.}
#' \item{TotLikSib}{Numeric vector, Total likelihood of the genotype data
#'   at initiation and after each iteration during Sibship clustering.}
#' \item{AgePriorExtra}{As AgePriors, but including columns for grandparents
#'  and avuncular pairs. NOT updated after parentage assignment, but returned
#'  as used during the run.}
#'
#' List elements PedigreePar and Pedigree both have the following columns:
#'  \item{id}{Individual ID}
#'  \item{dam}{Assigned mother, or NA}
#'  \item{sire}{Assigned father, or NA}
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
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. (2017) Pedigree reconstruction from SNP data:
#'   Parentage assignment, sibship clustering, and beyond. Molecular Ecology
#'   Resources 17:1009--1024.
#'
#' @section Disclaimer:
#' While every effort has been made to ensure that sequoia provides what it
#' claims to do, there is absolutely no guarantee that the results provided are
#' correct. Use of sequoia is entirely at your own risk.
#'
#' @seealso \code{\link{GenoConvert}} to read in various data formats,
#'   \code{\link{CheckGeno}}, \code{\link{SnpStats}} to calculate missingness
#'   and allele frequencies, \code{\link{MakeAgePrior}} to estimate effect of
#'   age on relationships, \code{\link{GetMaybeRel}} to find pairs of potential
#'   relatives, \code{\link{SummarySeq}} and \code{\link{PlotAgePrior}} to
#'   visualise results, \code{\link{GetRelCat}} to turn a pedigree into pairwise
#'   relationships, \code{\link{CalcOHLLR}} to calculate OH and LLR,
#'   \code{\link{PedCompare}} and \code{\link{ComparePairs}} to compare to a
#'   previous pedigree, \code{\link{EstConf}} and \code{\link{SimGeno}} to
#'   estimate assignment errors, \code{\link{writeSeq}} to save results,
#'   vignette("sequoia") for further details & FAQ.
#'
#' @examples
#' # ===  EXAMPLE 1: simulate data  ===
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#' head(SimGeno_example[,1:10])
#' head(LH_HSg5)
#' SeqOUT <- sequoia(GenoM = SimGeno_example, Err = 0.005,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 0)
#' names(SeqOUT)
#' SeqOUT$PedigreePar[34:42, ]
#'
#' # compare to true (or old) pedigree:
#' PC <- PedCompare(Ped_HSg5, SeqOUT$PedigreePar)
#' PC$Counts["GG",,]
#'
#' \dontrun{
#' SeqOUT2 <- sequoia(GenoM = SimGeno_example, Err = 0.005,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 10)
#' SeqOUT2$Pedigree[34:42, ]
#'
#' PC2 <- PedCompare(Ped_HSg5, SeqOUT2$Pedigree)
#' PC2$Counts["GT",,]
#'
#' # important to run with (approx.) correct genotyping error rate:
#' SeqOUT2.b <- sequoia(GenoM = SimGeno_example, #  Err = 1e-4 by default
#'                   LifeHistData = LH_HSg5, MaxSibIter = 10)
#' PC2.b <- PedCompare(Ped_HSg5, SeqOUT2.b$Pedigree)
#' PC2.b$Counts["GT",,]
#'
#'
#' # ===  EXAMPLE 2: real data  ===
#' # ideally, select 400-700 SNPs: high MAF & low LD
#' # save in 0/1/2/NA format (PLINK's --recodeA)
#' GenoM <- GenoConvert(InFile = "inputfile_for_sequoia.raw",
#'                      InFormat = "raw")  # can also do Colony format
#' SNPSTATS <- SnpStats(GenoM)
#' # perhaps after some data-cleaning:
#' write.table(GenoM, file="MyGenoData.txt", row.names=T, col.names=F)
#'
#' # later:
#' GenoM <- as.matrix(read.table("MyGenoData.txt", row.names=1, header=F))
#' LHdata <- read.table("LifeHistoryData.txt", header=T) # ID-Sex-birthyear
#' SeqOUT <- sequoia(GenoM, LHdata, Err=0.005)
#' SummarySeq(SeqOUT)
#'
#' writeSeq(SeqOUT, folder="sequoia_output")  # several text files
#'
#' # runtime:
#' SeqOUT$Specs$TimeEnd - SeqOUT$Specs$TimeStart
#' }
#'
#' @export

sequoia <- function(GenoM = NULL,
                    LifeHistData = NULL,
                    SeqList = NULL,
                    MaxSibIter = 10,
                    Err = 0.0001,
                    ErrFlavour = "version2.0",
                    MaxMismatch = NA,  # DEPRECATED
                    Tfilter = -2.0,
                    Tassign = 0.5,
                    MaxSibshipSize = 100,
                    DummyPrefix = c("F", "M"),
                    Complex = "full",
                    UseAge = "yes",
                    args.AP = list(Flatten=NULL, Smooth=TRUE),
                    FindMaybeRel = FALSE,
                    CalcLLR = TRUE,
                    quiet = FALSE,
                    Plot = NULL)
{
  if (!is.logical(quiet)) {
    if (quiet == "verbose") {
      quietF = -1
      quiet = FALSE
    } else {
      stop("`quiet' must be TRUE/FALSE or 'verbose'")
    }
  } else {
    quietF <- as.numeric(quiet)
  }
  if (!is.null(Plot) && !is.logical(Plot))  stop("'Plot' must be TRUE/FALSE")
  if (is.null(Plot)) {
    Plot <- ifelse(quiet, FALSE, TRUE)
  }

  Excl <- CheckGeno(GenoM, quietF, Plot)
  if ("ExcludedSnps" %in% names(Excl))  GenoM <- GenoM[, -Excl[["ExcludedSnps"]]]
  if ("ExcludedSnps-mono" %in% names(Excl))  GenoM <- GenoM[, -Excl[["ExcludedSnps-mono"]]]
  if ("ExcludedIndiv" %in% names(Excl))  GenoM <- GenoM[!rownames(GenoM) %in% Excl[["ExcludedIndiv"]], ]

  if (!"Specs" %in% names(SeqList)) {
    if (is.null(LifeHistData)) {
      warning("no LifeHistData provided, expect lower assignment rate",
              immediate.=TRUE)
    } else if (all(is.na(LifeHistData))) {
      stop("invalid value for LifeHistData, provide NULL or dataframe")
    } else if (length(intersect(LifeHistData[,1], rownames(GenoM)))==0) {
      warning("none of the genotyped individuals included in lifehistory data",
              immediate.=TRUE)
    }
  }
  if (!is.na(MaxMismatch)) {
    warning("NOTE: 'MaxMismatch' is deprecated & ignored; now calculated internally by 'CalcMaxMismatch'",
            immediate.=TRUE)
  }
  if (FindMaybeRel)  message("NOTE: 'FindMaybeRel' will be deprecated, please run GetMaybeRel() separately")

  if ("LifeHist" %in% names(SeqList)) {
    if(!quiet)  message("using LifeHistData in SeqList")
    LifeHistData <- SeqList$LifeHist
  }
  if (!is.null(LifeHistData)) {
    LifeHistData <- CheckLH(LifeHistData)
    if (!quiet & (any(LifeHistData$BY.min>0 | LifeHistData$BY.max > 0))) {
      message("reading minimum/maximum birth years in LifeHistData")
    }
  } else {
    LifeHistData <- data.frame(ID = rownames(GenoM),
                               Sex = 3,
                               BirthYear = -999,
                               BY.min = -999,
                               BY.max = -999,
                               stringsAsFactors = FALSE)
  }
  utils::flush.console()

  if ("Specs" %in% names(SeqList)) {
    if(!quiet)  message("settings in SeqList will overrule other settings")

    if (!"MaxMismatchOH" %in% names(SeqList$Specs)) {  # backwards compatability, version 1.x
      ErrM <- ErrToM(FacToNum(SeqList$Specs["GenotypingErrorRate"]),
                     flavour = ErrFlavour, Return = "matrix")
			sts <- SnpStats(GenoM, Plot=FALSE)
      MaxMismatchV <- CalcMaxMismatch(Err=ErrM, MAF=sts[,"AF"], ErrFlavour=ErrFlavour,
                                      qntl=0.999^(1/nrow(GenoM)))  # length 3 vector
    } else {
      MaxMismatchV <- setNames(FacToNum(SeqList$Specs[c("MaxMismatchDUP", "MaxMismatchOH",
                                                        "MaxMismatchME")]),
                               c("DUP", "OH", "ME"))
    }
    Specs <- SeqPrep(GenoM = GenoM,
                     LifeHistData = LifeHistData,
                     nAgeClasses = FacToNum(SeqList$Specs["nAgeClasses"]),
                     MaxSibIter = MaxSibIter,   # take new value
                     Err = FacToNum(SeqList$Specs["GenotypingErrorRate"]),
                     MaxMismatchV = MaxMismatchV,
                     Tfilter = FacToNum(SeqList$Specs["Tfilter"]),
                     Tassign = FacToNum(SeqList$Specs["Tassign"]),
                     MaxSibshipSize = FacToNum(SeqList$Specs["MaxSibshipSize"]),
                     DummyPrefix = as.character(SeqList$Specs[c("DummyPrefixFemale",
                                                   "DummyPrefixMale")]),
                     Complexity = as.character(SeqList$Specs["Complexity"]),
                     UseAge = as.character(SeqList$Specs["UseAge"]),
                     FindMaybeRel = as.logical(SeqList$Specs["FindMaybeRel"]),
                     CalcLLR = as.logical(SeqList$Specs["CalcLLR"]))
  } else {
    ErrM <- ErrToM(Err, flavour = ErrFlavour, Return = "matrix")  # performs checks if Err is already matrix
		sts <- SnpStats(GenoM, Plot=FALSE)
		MaxMismatchV <- CalcMaxMismatch(Err=ErrM, MAF=sts[,"AF"], ErrFlavour=ErrFlavour,
		                                qntl=0.999^(1/nrow(GenoM)))  # length 3 vector

    Specs <- SeqPrep(GenoM = GenoM,
                     LifeHistData = LifeHistData,
                     nAgeClasses = 1,
                     MaxSibIter = as.numeric(MaxSibIter),
                     Err = ifelse(length(Err)==1,
                                  Err,
                                  signif(1 - mean(diag(ErrM)), 3)),  # 1 - mean prob. correct (unweighed)
                     MaxMismatchV = MaxMismatchV,
                     Tfilter = Tfilter,
                     Tassign = Tassign,
            				 MaxSibshipSize = MaxSibshipSize,
            				 DummyPrefix = DummyPrefix,
            				 Complexity = Complex,
            				 UseAge = UseAge,
            				 FindMaybeRel = FindMaybeRel,
            				 CalcLLR = CalcLLR)
  }
  Specs$TimeStart <- Sys.time()
  utils::flush.console()

  if ("AgePriors" %in% names(SeqList)) {
    if(!quiet)  message("using ageprior in SeqList")
    AgePriors <- SeqList$AgePriors
  } else {
    AgePriors <- do.call(MakeAgePrior, c(list(Pedigree = NULL,
                                              LifeHistData = LifeHistData,
                                              Plot = Plot,
                                              quiet = quiet),
                                         args.AP))
  }
  Specs[,"nAgeClasses"] <- nrow(AgePriors)
	if (!all(apply(AgePriors, 2, function(x) any(x > 0)))) {
		stop("AgePriors error: some relationships are impossible for all age differences")
	}

  if (MaxSibIter>-9) {
		DupList <- DuplicateCheck(GenoM = GenoM, Specs = Specs, ErrM = ErrM,
		                          LhIN = LifeHistData[,1:5], quiet=quietF)
    utils::flush.console()
    if ("DupGenoID" %in% names(DupList)) {
      return(DupList)
    }
  } else DupList <- NULL
  utils::flush.console()

  if ("PedigreePar" %in% names(SeqList)) {
    PedParents <- SeqList$PedigreePar
    n.shared.ids <- length(intersect(PedParents[,1], rownames(GenoM)))
    if (n.shared.ids==0) {
      stop("GenoM and SeqList$PedigreePar do not share any common individuals")
    } else if (n.shared.ids < nrow(GenoM)/10 && n.shared.ids < nrow(PedParents)) {
      warning("GenoM and SeqList$PedigreePar share few common individuals")
    }
    PedParents <- PedPolish(PedParents, GenoNames = rownames(GenoM), ZeroToNA = TRUE)
  }

  if (any(LifeHistData$Sex==4)) {  # hermaphrodites - pretend 2 clones of opposite sex
    if (!quietF) message("detected hermaphrodites (sex=4), changing Complex to 'herm'")
    GenoM <- herm_clone_Geno(GenoM, LifeHistData, herm.suf=c("f", "m"))
    LifeHistData <- herm_clone_LH(LifeHistData, herm.suf=c("f", "m"))
    Specs[1, "NumberIndivGenotyped"] <- nrow(GenoM)
    Specs[1, "Complexity"] <- "herm"

    if ("PedigreePar" %in% names(SeqList)) {
      PedParents <- herm_clone_Ped(Ped = PedParents, LH = LifeHistData[,1:3], herm.suf=c("f", "m"))
      PedParents <- PedParents[match(rownames(GenoM),PedParents[,1]), ]
    }
  }

  if ("PedigreePar" %in% names(SeqList)) {
    if (MaxSibIter>0) {
      if(!quiet) message("using parents in SeqList")
      if (!"AgePriors" %in% names(SeqList) && !is.null(LifeHistData)) {
        AgePriors <- do.call(MakeAgePrior, c(list(Pedigree = PedParents[, 1:3],
                                                  LifeHistData = LifeHistData,
                                                  Plot = Plot,
                                                  quiet = quiet),
                                             args.AP))
        if (nrow(AgePriors) != Specs[,"nAgeClasses"]) {
          Specs[,"nAgeClasses"] <- nrow(AgePriors)
        }
      }
      ParList <- list(PedigreePar = PedParents)
      if ("TotLikParents" %in% names(SeqList)) {
        ParList <- c(ParList, SeqList["TotLikParents"])
      }
      if ("MaybeParent" %in% names(SeqList)) {
        ParList <- c(ParList, SeqList["MaybeParent"])
      }
    } else if (MaxSibIter==0) {  # re-run parentage with pedigree-prior
      ParList <- SeqParSib(ParSib = "par", Specs = Specs, ErrM = ErrM, GenoM = GenoM,
                           LhIN = LifeHistData[,1:5], AgePriors = AgePriors,
                           Parents=PedParents, quiet=quietF)
			if (Plot)  SummarySeq(ParList$PedigreePar, Panels="G.parents")
    } # else keep old parents

  } else if (MaxSibIter>=0){
    ParList <- SeqParSib(ParSib="par", Specs=Specs,  ErrM = ErrM, GenoM=GenoM,
                         LhIN=LifeHistData[,1:5],
                         AgePriors=AgePriors, Parents=NULL, quiet=quietF)
    if (Plot)  SummarySeq(ParList$PedigreePar, Panels="G.parents")
    if (!"AgePriors" %in% names(SeqList) && !is.null(LifeHistData)) {
      if (Plot)  Sys.sleep(1)  # else no time to notice previous plot
      AgePriors <- do.call(MakeAgePrior, c(list(Pedigree = ParList$PedigreePar[, 1:3],
                                                LifeHistData = LifeHistData,
                                                Plot = Plot,
                                                quiet = quiet),
                                           args.AP))
      if (nrow(AgePriors) != Specs[,"nAgeClasses"]) {
        Specs[,"nAgeClasses"] <- nrow(AgePriors)
      }
    } else if ("AgePriors" %in% names(SeqList)) {
      if(!quiet)  message("using ageprior in SeqList")
    }
  } else if (MaxSibIter <0) ParList <- NULL
  utils::flush.console()
	if (!all(apply(AgePriors, 2, function(x) any(x > 0)))) {
		stop("AgePriors error: some relationships are impossible for all age differences")
	}

  if (MaxSibIter>0) {
    SibList <- SeqParSib(ParSib = "sib", Specs = Specs, ErrM = ErrM, GenoM = GenoM,
                         LhIN = LifeHistData[,1:5], AgePriors = AgePriors,
                         Parents = ParList$PedigreePar, quiet = quietF)
    ParList <- ParList[names(ParList) != "AgePriorExtra"]  # else included 2x w same name
    if (Plot) {
      PlotAgePrior(SibList$AgePriorExtra)
      SummarySeq(SibList$Pedigree, Panels="G.parents")
    }
  } else SibList <- NULL

	# non-assigned putative relatives
  if (Specs[,"FindMaybeRel"]) {
		if (MaxSibIter > 0) {
			Ped <- SibList$Pedigree
			ParSib = "sib"
		} else if (MaxSibIter == 0) {
			Ped <- ParList$PedigreePar
			ParSib = "par"
		} else {
			Ped <- NULL
			ParSib = "par"
		}
    MaybeRelPairTrio <- GetMaybeRel(Pedigree = Ped,
                                    GenoM = GenoM,
                                    SeqList = list(Specs = Specs,
                                                   LifeHist = LifeHistData[,1:5],
                                                   AgePriors = AgePriors),
                                    ParSib = ParSib,
                                    quiet = quiet)
  } else   MaybeRelPairTrio <- NULL

  if (MaxSibIter>=0 & any(LifeHistData$Sex==4)) {
    ParList$PedigreePar <- herm_unclone_Ped(Ped = ParList$PedigreePar, LH = LifeHistData, herm.suf=c("f", "m"))
    if (MaxSibIter>0) {
      SibList$Pedigree <- herm_unclone_Ped(Ped = SibList$Pedigree, LH = LifeHistData, herm.suf=c("f", "m"))
    }
  }

  #=====================
  Specs$TimeEnd <- Sys.time()
  OUT <- list()
  OUT[["Specs"]] <- Specs
  OUT[["ErrM"]] <- ErrM
  if (length(Excl)>0)  OUT <- c(OUT, Excl)
  OUT[["AgePriors"]] <- AgePriors
  OUT[["LifeHist"]] <- LifeHistData
  OUT <- c(OUT, DupList, ParList)
  if (MaxSibIter>0)  OUT <- c(OUT, SibList)
  if (Specs[,"FindMaybeRel"])  OUT <- c(OUT, MaybeRelPairTrio)
  return(OUT)
}
