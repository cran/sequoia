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
#' When calling with GenoFormat="raw" or GenoFormat="col", the genotype file is
#' converted and written to a temporary directory. For large files, this may be
#'  time consuming; please consider calling \code{\link{GenoConvert}} directly
#'  before pedigree reconstruction.
#'
#' Further explanation of the various options and interpretation of the output
#'  is provided in the vignette.
#'
#' @param GenoM  numeric matrix with genotype data: One row per individual, and
#'   one column per SNP, coded as 0, 1, 2 or -9 (missing). Use
#'   \code{\link{GenoConvert}} to convert genotype files created in PLINK using
#'   --recodeA or in Colony's 2-column format to this format.
#' @param GenoFile character string with name of genotype file, with no header.
#'   Assumed to be in current working directory unless complete path is
#'   provided, which may not contain any spaces and must be less than 255
#'   characters long.
#' @param LifeHistData Dataframe with 3 columns:
#'  \itemize{
#'  \item{ID: }{max. 30 characters long,}
#'  \item{Sex: }{1 = females, 2 = males, other numbers = unkown,}
#'  \item{Birth Year: }{(or hatching year) Negative numbers are
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
#' @param quiet suppress messages
#'
#' @return A list with some or all of the following components:
#' \item{AgePriors}{Matrix with age-difference based prior probabilities}
#' \item{DummyIDs}{Dataframe with pedigree for dummy individuals, as well as
#'   their sex, estimated birth year (point estimate, upper and lower bound of
#'   95\% confidence interval), number of offspring, and offspring IDs.}
#' \item{DupGenoID}{Dataframe, rownumbers of duplicated IDs in genotype file.
#'   Please do remove or relabel these to avoid downstream confusion.}
#' \item{DupGenotype}{Dataframe, duplicated genotypes (with or without
#'   identical IDs). The specified number of maximum mismatches is allowed,
#'   and this dataframe may include pairs of closely related individuals.}
#' \item{DupLifeHistID}{Dataframe, rownumbers of duplicated IDs in life
#'   history file}
#' \item{LifeHist}{Provided dataframe with sex and birth year data}
#' \item{MaybeParent}{Dataframe with pairs of individuals who are more likely
#'   parent-offspring than unrelated, but which could not be phased due to
#'   unknown age difference (coded as 999) or sex, or for whom LLR did not pass
#'   Tassign}
#' \item{MaybeRel}{Dataframe with pairs of individuals who are more likely
#'   to be first or second degree relatives than unrelated, but which could not
#'   be assigned.}
#' \item{NoLH}{Vector, IDs (in genotype file) for which no life history data is
#'  provided}
#' \item{Pedigree}{Dataframe with assigned genotyped and dummy parents from
#'   Sibship step; entries for dummy individuals are added at the bottom.}
#' \item{PedigreePar}{Dataframe with assigned parents from Parentage step}
#' \item{Specs}{Named vector with parameter values}
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
#'  \item{LLR_dam}{Log10-Likelihood Ratio (LLR) of this female being the mother,
#'    versus the next most likely relationship between the focal individual and
#'    this female (see Details for relationships considered)}
#'  \item{LLR_sire}{idem, for male parent}
#'  \item{LLR_pair}{LLR for the parental pair, versus the next most likely
#'   configuration between the three individuals (with one or neither parent
#'   assigned)}
#' In addition, PedigreePar has the columns
#'  \item{OH_dam}{Number of loci at which the offspring and mother are
#'    opposite homozygotes}
#'  \item{OH_sire}{idem, for male parent}
#'  \item{rowID}{Row number in genotype file for id}
#'  \item{rowDam}{Row number in genotype file for dam}
#'  \item{rowSire}{Row number in genotype file for sire}
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. Pedigree reconstruction from SNP data: Parentage
#'   assignment, sibship clustering, and beyond. (under review) Molecular
#'   Ecology Resources
#'
#' @seealso \code{\link{GenoConvert}, \link{SimGeno}, \link{PedCompare}}
#'
#' @examples
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#'
#' SeqOUT <- sequoia(GenoM = SimGeno_example,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 0)
#' names(SeqOUT)
#' SeqOUT$PedigreePar[34:42, ]
#' \dontrun{
#' SeqOUT2 <- sequoia(GenoM = SimGeno_example,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 10)
#' SeqOUT2$Pedigree[34:42, ]
#' }
#' @export

sequoia <- function(GenoM = NULL,
                    GenoFile = NULL,
                    LifeHistData = NULL,
                    SeqList = NULL,
                    MaxSibIter = 10,
                    Err = 0.0001,
                    MaxMismatch = 3,
                    Tfilter = -2.0,
                    Tassign = 1.0,
                    MaxSibshipSize = 100,
                    DummyPrefix = c("F", "M"),
                    Complex = "full",
                    quiet = FALSE)
{
  if (!"Specs" %in% names(SeqList)) {
    if (is.null(GenoM) & is.null(GenoFile)) stop("please provide 'GenoM' or 'GenoFile'")
    if (!is.null(GenoM) & !is.matrix(GenoM)) stop("'GenoM' should be a numeric matrix")
    if (length(GenoFile)>1) stop("Provide single name for 'GenoFile'")
    if (is.null(LifeHistData)) warning("no lifehistory data provided; assuming single cohort")
    if (Err<0 | Err>1 | !is.double(Err)) stop("invalid value for 'Err'")
    if (MaxMismatch<0 | !is.wholenumber(MaxMismatch))  stop("invalid value for 'MaxMismatch'")
    if (!is.double(Tfilter))  stop("invalid value for 'Tfilter'")
    if (!is.double(Tassign) | Tassign<0)  stop("invalid value for 'Tassign'")
    if (MaxSibshipSize<0 | !is.wholenumber(MaxSibshipSize)) {
      stop("invalid value for MaxSibshipSize")
    }
    if (!Complex %in% c("full", "simp", "mono")) stop("invalid value for 'Complex'")
  }
  if (MaxSibIter<0 | !is.wholenumber(MaxSibIter))  stop("invalid value for 'MaxSibIter'")
  if (!is.logical(quiet)) stop("`quiet' must be TRUE/FALSE")


  if ("Specs" %in% names(SeqList)) {
    if(!quiet)  message("settings in SeqList will overrule other settings")
    # TODO: check if values are valid
    Specs <- SeqList$Specs
    Specs["MaxSibIter"] <- MaxSibIter
    if (Specs["GenotypeFilename"] != "USEMATRIX" | !is.null(GenoM)) {
      GenoFile <- Specs["GenotypeFilename"]
      if (nchar(GenoFile)>=255)  stop("'GenoFile' longer than 255 characters")
      if (length(grep(" ", GenoFile))>0)  stop("'GenoFile' may not contain spaces")
    } else {
      stop("please provide 'GenoM' or 'GenoFile'")
    }
    if (GenoFile != "USEMATRIX") {
      if (!file.exists(GenoFile)) stop("cannot find GenoFile")
    }
  } else {
    Specs <- SeqPrep(GenoFile, GenoM, LifeHistData,
                     MaxSibIter, Err, MaxMismatch, Tfilter, Tassign,
                     MaxSibshipSize, DummyPrefix, Complex)
  }

  if ("LifeHist" %in% names(SeqList)) {
    if(!quiet)  message("using LifeHistData in SeqList")
    LifeHistData <- SeqList$LifeHist
  }
  if (!is.null(LifeHistData)) {
    LifeHistData[, 1] <- as.character(LifeHistData[, 1])
    for (x in 2:3) LifeHistData[, x] <- as.numeric(LifeHistData[, x])
    names(LifeHistData) <- c("ID", "Sex", "BY")
  }

  if ("AgePriors" %in% names(SeqList)) {
    if(!quiet)  message("using ageprior in SeqList")
    AgePriors <- SeqList$AgePriors
  } else {
    AgePriors <- MakeAgeprior(UseParents = FALSE,
                              FacToNum(Specs["nAgeClasses"]))
  }

  DupList <- SeqDup(Specs, GenoM, LifeHistData, quiet)
  if ("DupGenoID" %in% names(DupList)) {
    return(DupList)
  }


  if ("PedigreePar" %in% names(SeqList)) {
    if (MaxSibIter>0) {
      if(!quiet) message("using parents in SeqList")

      if (!"AgePriors" %in% names(SeqList) & !is.null(LifeHistData)) {
        AgePriors <- MakeAgeprior(UseParents = TRUE,
                                  FacToNum(Specs["nAgeClasses"]),
                                  Parents = SeqList$PedigreePar[,1:3],
                                  LifeHistData)
      }
    }
    ParList <- NULL

  } else {
    ParList <- SeqParSib("par", Specs, GenoM, LifeHistData,
                         AgePriors, Parents=NULL, quiet)
    if (!is.null(LifeHistData)) {
      AgePriors <- MakeAgeprior(UseParents = TRUE,
                                FacToNum(Specs["nAgeClasses"]),
                                Parents = ParList$PedigreePar[,1:3],
                                LifeHistData)  # update agepriors
    }
  }


  if (MaxSibIter>0) {
    if ("PedigreePar" %in% names(SeqList)) {
      if (any(SeqList$PedigreePar$rowID != seq(1,nrow(SeqList$PedigreePar)))) {
        stop("rows in SeqList$PedigreePar changed, cannot use for sibship inference")
      }
      PedPar <- SeqList$PedigreePar[, c("rowDam", "rowSire")]
    } else {
      PedPar <- ParList$PedigreePar[, c("rowDam", "rowSire")]
    }
    SibList <- SeqParSib("sib", Specs, GenoM, LifeHistData,
                         AgePriors, PedPar, quiet)
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
