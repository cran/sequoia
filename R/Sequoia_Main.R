#' @title Pedigree Reconstruction
#'
#' @description Perform pedigree reconstruction based on SNP data, including
#'   parentage assignment and sibship clustering.
#'
#' @details For each pair of candidate relatives, the likelihoods are calculated
#'   of them being parent-offspring (PO), full siblings (FS), half siblings
#'   (HS), grandparent-grandoffspring (GG), full avuncular (niece/nephew -
#'   aunt/uncle; FA), half avuncular/great-grandparental/cousins (HA), or
#'   unrelated (U). Assignments are made if the likelihood ratio (LLR) between
#'   the focal relationship and the most likely alternative exceed the threshold
#'   Tassign.
#'
#'   Dummy parents of sibships are denoted by F0001, F0002, ... (mothers)
#'   and M0001, M0002, ... (fathers), are appended to the bottom of the
#'   pedigree, and may have been assigned real or dummy parents themselves (i.e.
#'   sibship-grandparents). A dummy parent is not assigned to singletons.
#'
#'   The genotyping error rate `Err` is by default at locus level, not allele
#'   level: the probability to observe true homozygote \emph{aa} as heterozygote
#'   \emph{Aa} is \eqn{\approx E}, and as alternate homozygote \emph{AA}
#'   \eqn{(E/2)^2}; the probability to observe a true heterozygote as \emph{aa}
#'   = the probability to observe it as \emph{AA} \eqn{= E/2}. This error
#'   structure can be fully customised by providing a 3x3 matrix of observed
#'   genotype (columns) conditional on actual genotype (rows) instead.
#'
#'   Full explanation of the various options and interpretation of the output is
#'   provided in the vignette.
#'
#' @param GenoM  numeric matrix with genotype data: One row per individual, and
#'   one column per SNP, coded as 0, 1, 2 or -9 (missing). See also
#'   \code{\link{GenoConvert}}.
#' @param LifeHistData dataframe with 3 columns (optionally 5):
#'  \describe{
#'  \item{ID}{max. 30 characters long}
#'  \item{Sex}{1 = female, 2 = male, 3 = unknown, 4 = hermaphrodite,
#'            other numbers or NA = unknown}
#' \item{BirthYear }{birth or hatching year, integer, with missing values as NA
#'   or any negative value.}
#' \item{BY.min}{minimum birth year, only used if BirthYear is missing}
#' \item{BY.max}{maximum birth year, only used if BirthYear is missing} }
#' If the species has multiple generations per year, use an integer coding such
#' that the candidate parents' `Birth year' is at least one smaller than their
#' putative offspring's. Column names are ignored, so ensure column order is ID
#' - sex - birth year (- BY.min - BY.max). Individuals do not need to be in the
#' same order as in `GenoM', nor do all genotyped individuals need to be
#' included.
#' @param SeqList list with output from a previous run, to be re-used in the
#'   current run. Used are elements `PedigreePar', `LifeHist', `AgePriors',
#'   `Specs', and `ErrM', and these override the corresponding input parameters.
#'   Not all of these elements need to be present, and all other elements are
#'   ignored. If \code{SeqList$Specs} is provided, all  input parameters with
#'   the same name as its items are ignored, except
#'   \code{Module}/\code{MaxSibIter}.
#' @param Module one of
#'   \describe{
#'     \item{pre}{Only input check, return \code{SeqList$Specs}}
#'     \item{dup}{Also check for duplicate genotypes}
#'     \item{par}{Also perform parentage assignment (genotyped parents to
#'       genotyped offspring)}
#'     \item{ped}{(Also) perform full pedigree reconstruction, including
#'       sibship clustering and grandparent assignment. By far the most time
#'       consuming, and may take several hours for large datasets.}
#'   }
#'   NOTE: \emph{Until `MaxSibIter` is fully deprecated: if `MaxSibIter` differs
#'   from the default (\code{42}), and `Module` equals the default
#'   (\code{'ped'}), MaxSibIter overrides `Module`.}
#' @param MaxSibIter \strong{[will be deprecated]} number of iterations of
#'   sibship clustering, including assignment of grandparents to sibships and
#'   avuncular relationships between sibships. Clustering continues until
#'   convergence or until MaxSibIter is reached. Set to 0 for parentage
#'   assignment only.
#' @param Err estimated genotyping error rate, as a single number or 3x3 matrix.
#'   Details below. The error rate is presumed constant across SNPs, and
#'   missingness is presumed random with respect to actual genotype.
#' @param ErrFlavour function that takes \code{Err} (single number) as input,
#'   and returns a 3x3 matrix of observed (columns) conditional on actual (rows)
#'   genotypes, or choose from inbuilt options 'version2.0', 'version1.3', or
#'   'version1.1', referring to the sequoia version in which they were the
#'   default. Ignored if \code{Err} is a matrix. See \code{\link{ErrToM}}.
#' @param MaxMismatch \strong{DEPRECATED AND IGNORED}. Now calculated
#'   automatically using \code{\link{CalcMaxMismatch}}.
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
#'   (mothers) and sires (fathers); maximum 20 characters each. Length 3 vector
#'   in case of hermaphrodites (or default prefix 'H').
#' @param Complex  Breeding system complexity. Either "full" (default), "simp"
#'   (simplified, no explicit consideration of inbred relationships), "mono"
#'   (monogamous).
#' @param Herm  Hermaphrodites, either "no", "A" (distinguish between dam and
#'   sire role, default if at least 1 individual with sex=4), or "B" (no
#'   distinction between dam and sire role). Both of the latter deal with
#'   selfing.
#' @param UseAge  either "yes" (default), "no", or "extra" (additional rounds
#'   with extra reliance on ageprior, may boost assignments but increased risk
#'   of erroneous assignments); used during full reconstruction only.
#' @param args.AP list with arguments to be passed on to
#'   \code{\link{MakeAgePrior}}.
#' @param FindMaybeRel \strong{DEPRECATED AND IGNORED}, advised to run
#'   \code{\link{GetMaybeRel}} separately. TRUE/FALSE to identify pairs of
#'   non-assigned likely relatives after pedigree reconstruction. Can be
#'   time-consuming in large datasets.
#' @param CalcLLR  TRUE/FALSE; calculate log-likelihood ratios for all assigned
#'   parents (genotyped + dummy; parent vs. otherwise related). Time-consuming
#'   in large datasets. Can be done separately with \code{\link{CalcOHLLR}}.
#' @param quiet suppress messages: TRUE/FALSE/"verbose".
#' @param Plot display plots from \code{\link{SnpStats}, \link{MakeAgePrior}},
#'   and \code{\link{SummarySeq}}. Defaults (NULL) to TRUE when quiet=FALSE or
#'   "verbose", and FALSE when quiet=TRUE. If you get error 'figure margins too
#'   large', enlarge the plotting area (drag with mouse). Error 'invalid
#'   graphics state' can be dealt with by clearing the plotting area with
#'   dev.off().
#'
#' @return A list with some or all of the following components:
#' \item{AgePriors}{Matrix with age-difference based probability ratios for
#'   each relationship, used for full pedigree reconstruction; see
#'   \code{\link{MakeAgePrior}} for details. When running only parentage
#'   assignment (\code{Module="par"}) the returned AgePriors has been updated to
#'   incorporate the information of the assigned parents, and is ready for use
#'   during full pedigree reconstruction.}
#' \item{DummyIDs}{Dataframe with pedigree for dummy individuals, as well as
#' their sex, estimated birth year (point estimate, upper and lower bound of
#' 95\% confidence interval; see also \code{\link{CalcBYprobs}}), number of
#' offspring, and offspring IDs. From version 2.1 onwards, this includes dummy
#' offspring.}
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
#' \item{DummyClones}{Hermaphrodites only: female-male dummy ID pairs that refer
#'   to the same non-genotyped individual}
#'
#' List elements PedigreePar and Pedigree both have the following columns:
#'  \item{id}{Individual ID}
#'  \item{dam}{Assigned mother, or NA}
#'  \item{sire}{Assigned father, or NA}
#'  \item{LLRdam}{Log10-Likelihood Ratio (LLR) of this female being the mother,
#'  versus the next most likely relationship between the focal individual and
#'  this female. See Details below for relationships considered, and see
#'  \code{\link{CalcPairLL}} for underlying likelihood values and further
#'  details)}
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
#' @seealso
#' \itemize{
#'   \item \code{\link{GenoConvert}} to read in various data formats,
#'   \item \code{\link{CheckGeno}}, \code{\link{SnpStats}} to calculate
#'     missingness and allele frequencies,
#'   \item \code{\link{SimGeno}}  to simulate SNP data from a pedigree
#'  \item \code{\link{MakeAgePrior}} to estimate effect of age on relationships,
#'   \item \code{\link{GetMaybeRel}} to find pairs of potential relatives,
#'   \item \code{\link{SummarySeq}} and \code{\link{PlotAgePrior}} to visualise
#'   results,
#'   \item \code{\link{GetRelM}} to turn a pedigree into pairwise relationships,
#'   \item \code{\link{CalcOHLLR}} to calculate Mendelian errors and LLR for any
#'    pedigree,
#'   \item \code{\link{CalcPairLL}} for likelihoods of various relationships
#'    between specific pairs,
#'   \item \code{\link{CalcBYprobs}} to estimate birth years,
#'   \item \code{\link{PedCompare}} and \code{\link{ComparePairs}} to compare to
#'   two pedigrees,
#'   \item \code{\link{EstConf}} to estimate assignment errors,
#'   \item \code{\link{writeSeq}} to save results,
#'   \item vignette("sequoia") for detailed manual & FAQ.
#' }
#'
#' @examples
#' # ===  EXAMPLE 1: simulated data  ===
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#' head(SimGeno_example[,1:10])
#' head(LH_HSg5)
#' # parentage assignment:
#' SeqOUT <- sequoia(GenoM = SimGeno_example, Err = 0.005,
#'                   LifeHistData = LH_HSg5, Module="par", Plot=TRUE)
#' names(SeqOUT)
#' SeqOUT$PedigreePar[34:42, ]
#'
#' # compare to true (or old) pedigree:
#' PC <- PedCompare(Ped_HSg5, SeqOUT$PedigreePar)
#' PC$Counts["GG",,]
#'
#' \donttest{
#' # parentage assignment + full pedigree reconstruction:
#' SeqOUT2 <- sequoia(GenoM = SimGeno_example, Err = 0.005,
#'                   LifeHistData = LH_HSg5, Module="ped", quiet="verbose")
#' SeqOUT2$Pedigree[34:42, ]
#'
#' PC2 <- PedCompare(Ped_HSg5, SeqOUT2$Pedigree)
#' PC2$Counts["GT",,]
#' PC2$Counts[,,"dam"]
#'
#' # different kind of pedigree comparison:
#' ComparePairs(Ped1=Ped_HSg5, Ped2=SeqOUT$PedigreePar, patmat=TRUE)
#'
#' # results overview:
#' SummarySeq(SeqOUT2)
#'
#' # important to run with approx. correct genotyping error rate:
#' SeqOUT2.b <- sequoia(GenoM = SimGeno_example, #  Err = 1e-4 by default
#'                   LifeHistData = LH_HSg5, Module="ped", Plot=FALSE)
#' PC2.b <- PedCompare(Ped_HSg5, SeqOUT2.b$Pedigree)
#' PC2.b$Counts["GT",,]
#' }
#' \dontrun{
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
                    Module = "ped",
                    MaxSibIter = 42,
                    Err = 0.0001,
                    ErrFlavour = "version2.0",
                    MaxMismatch = NA,  # DEPRECATED
                    Tfilter = -2.0,
                    Tassign = 0.5,
                    MaxSibshipSize = 100,
                    DummyPrefix = c("F", "M"),
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP = list(Flatten=NULL, Smooth=TRUE),
                    FindMaybeRel = FALSE,  # DEPRECATED
                    CalcLLR = TRUE,
                    quiet = FALSE,
                    Plot = NULL)
{

  TimeStart <- Sys.time()

  # set quiet & Plot ----
  if (!quiet %in% c(TRUE, FALSE, "verbose"))
    stop("'quiet' must be TRUE or FALSE or 'verbose'")
  quietR <- ifelse(quiet == "verbose", FALSE, quiet)  # 'verbose' determines chattiness of Fortran only
  if (is.null(Plot))   # default
    Plot <- ifelse(quietR, FALSE, TRUE)


  # Backwards compatibility: Module vs MaxSibIter ----
  # if MaxSibIter is not default (42), and Module is default (ped), MaxSibIter overrides Module
  if (MaxSibIter != 42 && Module == "ped") {
    Module <- cut(MaxSibIter,
                  breaks= c(-Inf, -9, -1, 0, Inf),
                  labels = c("pre", "dup", "par", "ped"))
    if (!quietR)  message("NOTE: 'MaxSibIter' will be deprecated in the future, ",
      "please consider using 'Module' instead")
  } else {
    Module <- factor(Module, levels = c("pre", "dup", "par", "ped"))
    if (is.na(Module))  stop("'Module'  must be 'pre', 'dup', 'par', or 'ped'")
  }

  if (FindMaybeRel)
    warning("'FindMaybeRel' has been deprecated and is ignored,",
            "instead run GetMaybeRel() afterwards", immediate.=TRUE)
  if (!is.na(MaxMismatch))
    warning("'MaxMismatch' has been deprecated and is ignored,",
            "now calculated automatically via CalcMaxMismatch()", immediate.=TRUE)

  if (!is.null(SeqList)) {
    if (any(is.na(names(SeqList))))
      stop("You seem to have misspelled one or more names of elements of SeqList")
  }


  # Check genotype matrix ----
  Excl <- CheckGeno(GenoM, quietR, Plot, Return = "excl", DumPrefix=DummyPrefix)
  if ("ExcludedSnps" %in% names(Excl))  GenoM <- GenoM[, -Excl[["ExcludedSnps"]]]
  if ("ExcludedSnps-mono" %in% names(Excl))  GenoM <- GenoM[, -Excl[["ExcludedSnps-mono"]]]
  if ("ExcludedIndiv" %in% names(Excl))  GenoM <- GenoM[!rownames(GenoM) %in% Excl[["ExcludedIndiv"]], ]

  if (!is.null(SeqList) && (!is.list(SeqList) | is.data.frame(SeqList)))  # dataframes are lists
    stop("SeqList must be a list or NULL")


  # Check life history data ----
  if ("LifeHist" %in% names(SeqList)) {
    if (!quietR)  message("using LifeHistData in SeqList")
    LifeHistData <- SeqList$LifeHist
  } else if (is.null(LifeHistData) & Module != "dup") {
    warning("no LifeHistData provided, expect lower assignment rate\n",
            immediate.=TRUE)
  }
  ChkLH.L <- CheckLH(LifeHistData, gID = rownames(GenoM),  sorted=FALSE,
                     returnDups = TRUE)
  # keep non-genotyped IDs (for future reference); orderLH() called by SeqParSib()
  DupList <- ChkLH.L[c("DupLifeHistID", "NoLH")]
  LifeHistData <- ChkLH.L$LifeHistData   # duplicates removed, if any

  utils::flush.console()    # print all warnings thus far

  # Check/make ageprior ----
  if ("AgePriors" %in% names(SeqList)) {
    if(!quietR)  message("using AgePriors in SeqList")
    AgePriors <- CheckAP( SeqList[["AgePriors"]] )
  } else {
    AgePriors <- do.call(MakeAgePrior, c(list(Pedigree = SeqList[["PedigreePar"]],
                    LifeHistData = LifeHistData[LifeHistData$ID %in% rownames(GenoM),],
                                              Plot = Plot,
                                              quiet = ifelse(Module=="dup", TRUE, quietR)),
                                         args.AP))
  }


  # Check Pedigree ----
  if ("PedigreePar" %in% names(SeqList) & Module != "dup") {
    if (!quietR)  message("using PedigreePar in SeqList")
    PedParents <- PedPolish(SeqList[["PedigreePar"]], gID = rownames(GenoM),
                            ZeroToNA = TRUE, DropNonSNPd = TRUE)
  } else {
    PedParents <- NULL
  }

  utils::flush.console()


  # make list with parameter values (loose input or SeqList$Specs) ----
  if ("Specs" %in% names(SeqList)) {
    PARAM <- SpecsToParam(SeqList$Specs, SeqList$ErrM, ErrFlavour,  # re-wrapping + ErrToM check
                          dimGeno = dim(GenoM), Module, quiet) # overrule old values
    if (Module=="ped")  PARAM$MaxSibIter <- 42
  } else {
    if (is.logical(UseAge))   UseAge <- ifelse(UseAge, "yes", "no")
    if ((Herm != "no" | any(LifeHistData$Sex==4)) & length(DummyPrefix)==2)
      DummyPrefix <- c(DummyPrefix, "H")

    PARAM <- namedlist(dimGeno = dim(GenoM),
                       Err,
                       ErrFlavour,
                       Tfilter,
                       Tassign,
                       nAgeClasses = nrow(AgePriors),
                       MaxSibshipSize,
                       Module = as.character(Module),
                       MaxSibIter,
                       DummyPrefix,
                       Complex,
                       Herm,
                       UseAge,
                       CalcLLR,
                       quiet)
    PARAM$ErrM <- ErrToM(Err, flavour = ErrFlavour, Return = "matrix")
  }

  # MaxMismatch ----
  # vector with max. mismatches for duplicates, PO pairs, PPO trios
  if (!"MaxMismatchV" %in% names(PARAM)) {  # DUP/OH/ME from version 2.0 onwards
    sts <- SnpStats(GenoM, Plot=FALSE)
    PARAM$MaxMismatchV <- setNames(CalcMaxMismatch(Err=PARAM$ErrM,
                                                   MAF=sts[,"AF"],
                                                   ErrFlavour=PARAM$ErrFlavour,
                                                   qntl=0.999^(1/nrow(GenoM))),
                                   c("DUP", "OH", "ME"))
  }


  # hermaprhodites ----
  if (any(LifeHistData$Sex==4) && PARAM$Herm == "no") {  #!grepl("herm", PARAM$Complex)) {
    if (!quietR) message("\nDetected hermaphrodites (sex=4), changing Herm to 'A'\n")
#    PARAM$Complex <- "herm"
    PARAM$Herm <- "A"
  }
  if (PARAM$Herm == "B" && any(LifeHistData$Sex %in% c(1,2)))
    warning("Results may be inconsistent when combining Herm='B' with known-sex individuals")


  # check parameter values ----
  CheckParams(c(PARAM, list(Plot=Plot)))
  utils::flush.console()

  # turn into fortran-friendly parameters ----
  FortPARAM <- MkFortParams(PARAM, fun="main")  # SpecsInt, SpecsDbl, character choice to numeric


  # @@ 2 @@ Duplicate check ----
  if (Module != "pre") {
    DupList <- c(DuplicateCheck(GenoM, FortPARAM, quiet=quietR),
                 DupList)  # from CheckLH()
    utils::flush.console()
    if ("DupGenoID" %in% names(DupList)) {  # fatal error
      return(DupList)
    } else if (length(DupList)==0 & !quietR) {
      message("No duplicates found")
    }
  } else DupList <- NULL
  utils::flush.console()


  # @@ 3 @@ Parentage assignment ----
  if (Module == "par"  | (Module  == "ped" & is.null(PedParents))) {
    if(!quietR & Module == "par" & !is.null(PedParents)) {
      message("\n~~~ Parentage assignment with pedigree-prior ~~~")  # only sensible with Herm=A or B (?)
    } else if(!quietR) {
      message("\n~~~ Parentage assignment ~~~")
    }

    ParList <- SeqParSib(ParSib = "par", FortPARAM, GenoM,
                         LhIN=LifeHistData, AgePriors=AgePriors,
                         Parents=PedParents, DumPfx = PARAM$DummyPrefix, quiet=quietR)
    if (Plot) {
      SummarySeq(ParList$PedigreePar, Panels="G.parents")
    }

  } else if (Module != "dup" & "PedigreePar" %in% names(SeqList)) {  # don't include for 'dup'; confusing.
    ParList <- list(PedigreePar = PedParents)   # re-use assigned parents given as SeqList input (after polishing)
  } else {
    ParList <- NULL
  }
  utils::flush.console()

  W <- tryCatch.W.E(getGenerations(ParList$PedigreePar, StopIfInvalid=FALSE))$warning
  if (!is.null(W)) {
    print(W)
    if (Module=="ped")  warning("Cancelling full pedigree reconstruction.")
  }


  # Update ageprior ----
  if (Module %in% c("par", "ped") & !"AgePriors" %in% names(SeqList) & is.null(W) &
      !"PedigreePar" %in% names(SeqList) && !is.null(LifeHistData)) {
    if (Plot)  Sys.sleep(1)  # else no time to notice previous plot
    # in case MakeAgePrior throws error, do return parentage results:
    AgePriors <- tryCatch( do.call(MakeAgePrior,
                                   c(list(Pedigree = ParList$PedigreePar[, 1:3],
                                          LifeHistData = LifeHistData[LifeHistData$ID %in% rownames(GenoM),],
                                          Plot = Plot & Module=="ped",
                                          quiet = !(!quietR & Module=="ped")),
                                     args.AP)),
                           error = function(e) {
                             message("AgePrior: error!")
                             return(ParList)
                           })
  } else if ("AgePriors" %in% names(SeqList) & !"PedigreePar" %in% names(SeqList)) {
    if(!quietR)  message("using AgePriors in SeqList again")
  }

  if (nrow(AgePriors) != PARAM$nAgeClasses) {
    PARAM$nAgeClasses <- nrow(AgePriors)
    FortPARAM$SpecsInt[["nAgeCl"]] <- nrow(AgePriors)
  }


  # @@ 4 @@ Full pedigree reconstruction ----
  if (Module == "ped" & is.null(W)) {
    if (!all(apply(AgePriors, 2, function(x) any(x > 0))))
      stop("AgePriors error: some relationships are impossible for all age differences")
    if(!quietR)  message("\n~~~ Full pedigree reconstruction ~~~")
    SibList <- SeqParSib(ParSib = "sib", FortPARAM, GenoM,
                         LhIN = LifeHistData, AgePriors = AgePriors,
                         Parents = ParList$PedigreePar, DumPfx = PARAM$DummyPrefix,
                         quiet = quietR)
    ParList <- ParList[names(ParList) != "AgePriorExtra"]  # else included 2x w same name
    if (Plot) {
      PlotAgePrior(SibList$AgePriorExtra)
      SummarySeq(SibList$Pedigree, Panels="G.parents")
    }
  } else SibList <- NULL


  #=====================
  # Output ----
  OUT <- list()
  OUT[["Specs"]] <- ParamToSpecs(PARAM, TimeStart, ErrFlavour)
  OUT[["ErrM"]] <- PARAM$ErrM
  if (is.function(ErrFlavour))  OUT[["ErrFlavour"]] <- ErrFlavour
  if (length(Excl)>0)  OUT <- c(OUT, Excl)
  OUT[["AgePriors"]] <- AgePriors
  OUT[["LifeHist"]] <- LifeHistData
  if (as.numeric(Module) > 1)  OUT <- c(OUT, DupList)
  if (as.numeric(Module) > 2)  OUT <- c(OUT, ParList)
  if (as.numeric(Module) > 3)  OUT <- c(OUT, SibList)
  return(OUT)
}
