#' Pedigree Reconstruction
#'
#' Perform pedigree reconstruction based on SNP data, including parentage
#' assignment and sibship clustering.
#'
#' For each pair of candidate relatives, the likelihoods are calculated of them
#'  being parent-offspring (PO), full siblings (FS), half siblings (HS),
#'  grandparent-grandoffspring (GG), full avuncular (niece/nephew - aunt/uncle;
#'  FA), half avuncular/great-grandparental/cousins (HA), or unrelated (U).
#'  Assignments are made if the likelihood ratio (LLR) between the focal
#'  relationship and the most likely alternative exceed the threshold Tassign.
#'
#' Further explanation of the various options and interpretation of the output
#'  is provided in the vignette.
#'
#' @param GenoM  numeric matrix with genotype data: One row per individual, and
#'   one column per SNP, coded as 0, 1, 2 or -9 (missing). Use
#'   \code{\link{GenoConvert}} to convert genotype files created in PLINK using
#'   --recodeA or in Colony's 2-column format to this format.
#' @param LifeHistData Dataframe with 3 columns:
#'  \itemize{
#'  \item{ID: }{max. 30 characters long,}
#'  \item{Sex: }{1 = females, 2 = males, other = unkown, except 4 = hermaphrodite,}
#'  \item{BY: }{(birth or hatching year) Negative numbers are
#'    interpreted as missing values.}}
#'  If the species has multiple generations per year, use an integer coding
#'   such that the candidate parents' `Birth year' is at least one larger than
#'   their putative offspring.
#' @param SeqList list with output from a previous run, containing the elements
#'   `Specs', `AgePriors' and/or `PedigreePar', as described below, to be used
#'    in the current run. If SeqList$Specs is provided, all other input
#'    parameter values except MaxSibIter are ignored.
#' @param MaxSibIter Number of iterations of sibship clustering, including
#'   assignment of grandparents to sibships and avuncular relationships between
#'   sibships. Set to 0 to not (yet) perform this step, which is by far the
#'   most time consuming and may take several hours for large datasets.
#'   clustering is iterated until convergence or until MaxSibIter is reached.
#' @param Err Estimated genotyping error rate. The error model aims to deal
#'   with scoring errors typical for SNP arrays.
#' @param MaxMismatch Maximum number of loci at which candidate parent and
#'  offspring are allowed to be opposite homozygotes.
#' @param Tfilter Threshold log10-likelihood ratio (LLR) between a proposed
#'   relationship versus unrelated, to select candidate relatives. Typically a
#'   negative value, related to the fact that unconditional likelihoods are
#'   calculated during the filtering steps. More negative values may decrease
#'   non-assignment, but will increase computational time.
#' @param Tassign Minimum LLR required for acceptance of
#'  proposed relationship, relative to next most likely relationship. Higher
#'  values result in more conservative assignments. Must be zero or positive.
#' @param MaxSibshipSize  Maximum number of offspring for a single individual
#'  (a generous safety margin is advised).
#' @param DummyPrefix character vector of length 2 with prefixes for dummy
#'   dams (mothers) and sires (fathers); maximum 20 characters each.
#' @param Complex  Either "full" (default), "simp" (simplified, no explicit
#'  consideration of inbred relationships; not fully implemented yet) or
#'  "mono" (monogamous).
#' @param FindMaybeRel  Identify pairs of non-assigned likely relatives after
#'   pedigree reconstruction. Can be time-consuming in large datasets.
#' @param CalcLLR  Calculate log-likelihood ratios for all assignments. Can be
#'  time-consuming in large datasets.
#' @param quiet suppress messages.
#'
#' @return A list with some or all of the following components:
#' \item{AgePriors}{Matrix with age-difference based prior probabilities.}
#' \item{DummyIDs}{Dataframe with pedigree for dummy individuals, as well as
#'   their sex, estimated birth year (point estimate, upper and lower bound of
#'   95\% confidence interval), number of offspring, and offspring IDs.}
#' \item{DupGenoID}{Dataframe, rownumbers of duplicated IDs in genotype file.
#'   Please do remove or relabel these to avoid downstream confusion.}
#' \item{DupGenotype}{Dataframe, duplicated genotypes (with or without
#'   identical IDs). The specified number of maximum mismatches is allowed,
#'   and this dataframe may include pairs of closely related individuals.}
#' \item{DupLifeHistID}{Dataframe, rownumbers of duplicated IDs in life
#'   history dataframe.}
#' \item{LifeHist}{Provided dataframe with sex and birth year data.}
#' \item{MaybeParent}{Dataframe with pairs of individuals who are more likely
#'   parent-offspring than unrelated, but which could not be phased due to
#'   unknown age difference (coded as 999) or sex, or for whom LLR did not pass
#'   Tassign.}
#' \item{MaybeRel}{Dataframe with pairs of individuals who are more likely
#'   to be first or second degree relatives than unrelated, but which could not
#'   be assigned.}
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
#'
#'
#' List elements PedigreePar and Pedigree both have the following columns:
#'  \item{id}{Individual ID}
#'  \item{dam}{Assigned mother, or NA}
#'  \item{sire}{Assigned father, or NA}
#'  \item{LLRdam}{Log10-Likelihood Ratio (LLR) of this female being the mother,
#'    versus the next most likely relationship between the focal individual and
#'    this female (see Details for relationships considered)}
#'  \item{LLRsire}{idem, for male parent}
#'  \item{LLRpair}{LLR for the parental pair, versus the next most likely
#'   configuration between the three individuals (with one or neither parent
#'   assigned)}
#' In addition, PedigreePar has the columns
#'  \item{OHdam}{Number of loci at which the offspring and mother are
#'    opposite homozygotes}
#'  \item{OHsire}{idem, for father}
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. Pedigree reconstruction from SNP data: Parentage
#'   assignment, sibship clustering, and beyond. (accepted manuscript) Molecular
#'   Ecology Resources
#'
#' @seealso \code{\link{GenoConvert}, \link{SimGeno}, \link{PedCompare}}, vignette("sequoia")
#'
#' @examples
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#' head(SimGeno_example[,1:10])
#' head(LH_HSg5)
#' SeqOUT <- sequoia(GenoM = SimGeno_example,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 0)
#' names(SeqOUT)
#' SeqOUT$PedigreePar[34:42, ]
#' \dontrun{
#' SeqOUT2 <- sequoia(GenoM = SimGeno_example,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 10)
#' SeqOUT2$Pedigree[34:42, ]
#'
#' # reading in data from text files:
#' GenoM <- as.matrix(read.table("MyGenoData.txt", row.names=1, header=FALSE))
#' LH <- read.table("MyLifeHistData.txt", header=TRUE)
#' MySeqOUT <- sequoia(GenoM = GenoM, LifeHistData = LH)
#' }
#' @export

sequoia <- function(GenoM = NULL,
                    LifeHistData = NULL,
                    SeqList = NULL,
                    MaxSibIter = 10,
                    Err = 0.0001,
                    MaxMismatch = 3,
                    Tfilter = -2.0,
                    Tassign = 0.5,
                    MaxSibshipSize = 100,
                    DummyPrefix = c("F", "M"),
                    Complex = "full",
                    FindMaybeRel = TRUE,
                    CalcLLR = TRUE,
                    quiet = FALSE)
{
  if (is.null(GenoM)) stop("please provide 'GenoM'")
  if (!is.matrix(GenoM)) stop("'GenoM' should be a numeric matrix")
  if (!"Specs" %in% names(SeqList)) {
    if (is.null(LifeHistData)) warning("no lifehistory data provided; assuming single cohort")
  }
  if (!is.logical(quiet)) stop("`quiet' must be TRUE/FALSE")

  if ("LifeHist" %in% names(SeqList)) {
    if(!quiet)  message("using LifeHistData in SeqList")
    LifeHistData <- SeqList$LifeHist
  }
  if (!is.null(LifeHistData)) {
    names(LifeHistData) <- c("ID", "Sex", "BY")
    LifeHistData$ID <- as.character(LifeHistData$ID)
    for (x in c("Sex", "BY")) LifeHistData[, x] <- as.integer(LifeHistData[, x])
    LifeHistData$BY[is.na(LifeHistData$BY)] <- -999
    LifeHistData$Sex[is.na(LifeHistData$Sex)] <- 3
    LifeHistData$Sex[!LifeHistData$Sex %in% 1:4] <- 3
  }

  if (any(LifeHistData$Sex==4)) {  # hermaphrodites - pretend 2 clones of opposite sex
    GenoM <- herm_clone_Geno(GenoM, LifeHistData, herm.suf=c("f", "m"))
    LifeHistData <- herm_clone_LH(LifeHistData, herm.suf=c("f", "m"))
  }

  if ("Specs" %in% names(SeqList)) {
    if(!quiet)  message("settings in SeqList will overrule other settings")
    Specs <- SeqPrep(GenoM = GenoM,
                     LifeHistData = LifeHistData,
                     nAgeClasses = FacToNum(SeqList$Specs["nAgeClasses"]),
                     MaxSibIter = MaxSibIter,   # take novel value
                     Err = FacToNum(SeqList$Specs["GenotypingErrorRate"]),
                     MaxMismatch = FacToNum(SeqList$Specs["MaxMismatch"]),
                     Tfilter = FacToNum(SeqList$Specs["Tfilter"]),
                     Tassign = FacToNum(SeqList$Specs["Tassign"]),
                     MaxSibshipSize = FacToNum(SeqList$Specs["MaxSibshipSize"]),
                     DummyPrefix = as.character(SeqList$Specs[c("DummyPrefixFemale",
                                                   "DummyPrefixMale")]),
                     Complexity = as.character(SeqList$Specs["Complexity"]),
                     FindMaybeRel = as.logical(SeqList$Specs["FindMaybeRel"]),
                     CalcLLR = as.logical(SeqList$Specs["CalcLLR"]))
  } else {
    Specs <- SeqPrep(GenoM = GenoM,
                     LifeHistData = LifeHistData,
                     nAgeClasses = 1,
                     MaxSibIter = as.numeric(MaxSibIter),
                     Err =Err,
                     MaxMismatch = MaxMismatch,
                     Tfilter = Tfilter,
                     Tassign = Tassign,
            				 MaxSibshipSize = MaxSibshipSize,
            				 DummyPrefix = DummyPrefix,
            				 Complexity = Complex,
            				 FindMaybeRel = FindMaybeRel,
            				 CalcLLR = CalcLLR)
  }



  if ("AgePriors" %in% names(SeqList)) {
    if(!quiet)  message("using ageprior in SeqList")
    AgePriors <- SeqList$AgePriors
  } else {
    AgePriors <- MakeAgeprior(UseParents = FALSE,
                              FacToNum(Specs[,"nAgeClasses"]))
  }

  DupList <- SeqDup(Specs, GenoM, LifeHistData, quiet)
  if ("DupGenoID" %in% names(DupList)) {
    return(DupList)
  }

  if ("PedigreePar" %in% names(SeqList)) {
    PedParents <- SeqList$PedigreePar  #[, c("id", "dam", "sire")]
    for (x in 1:3)  PedParents[, x] <- as.character(PedParents[, x])
    if (any(LifeHistData$Sex==4)) {  # hermaphrodites
      PedParents <- herm_clone_Ped(PedParents, LifeHistData, herm.suf=c("f", "m"))
    }
    if (MaxSibIter>0) {
      if(!quiet) message("using parents in SeqList")
      if (!"AgePriors" %in% names(SeqList) && !is.null(LifeHistData)) {
        AgePriors <- MakeAgeprior(UseParents = TRUE,
                                  nAgeClasses = FacToNum(Specs[,"nAgeClasses"]),
                                  Parents = PedParents[,1:3],
                                  LifeHistData = LifeHistData)
      }
      ParList <- list(PedigreePar = PedParents)
      if ("TotLikParents" %in% names(SeqList)) {
        ParList <- c(ParList, SeqList[c("MaybeParent", "TotLikParents")])
      }
    } else {  # 'hidden' for now: re-run parentage with pedigree-prior
      ParList <- SeqParSib("par", Specs, GenoM, LifeHistData,
                         AgePriors, Parents=PedParents[,1:3], quiet)

    }
  } else {
    ParList <- SeqParSib("par", Specs, GenoM, LifeHistData,
                         AgePriors, Parents=NULL, quiet)
    if (!is.null(LifeHistData)) {
      AgePriors <- MakeAgeprior(UseParents = TRUE,
                                nAgeClasses = FacToNum(Specs[,"nAgeClasses"]),
                                Parents = ParList$PedigreePar[,1:3],
                                LifeHistData = LifeHistData)  # update agepriors
    }
  }

  if (MaxSibIter>0) {
    SibList <- SeqParSib(ParSib = "sib", Specs = Specs, GenoM = GenoM,
                         LhIN = LifeHistData, AgePriors = AgePriors,
                         Parents = ParList$PedigreePar[, 1:3], quiet = quiet)
  } else SibList <- NULL


  #=====================
  OUT <- list()
  OUT[["Specs"]] <- Specs
  OUT[["AgePriors"]] <- AgePriors
  OUT[["LifeHist"]] <- LifeHistData
  OUT <- c(OUT, DupList, ParList)
  if (MaxSibIter>0)  OUT <- c(OUT, SibList)
  return(OUT)
}
