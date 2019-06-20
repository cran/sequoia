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
#'   Further explanation of the various options and interpretation of the output
#'   is provided in the vignette.
#'
#' @param GenoM  numeric matrix with genotype data: One row per individual, and
#'   one column per SNP, coded as 0, 1, 2 or -9 (missing). Use
#'   \code{\link{GenoConvert}} to convert genotype files created in PLINK using
#'   --recodeA or in Colony's 2-column format to this format.
#' @param LifeHistData Dataframe with 3 columns:
#'  \itemize{
#'  \item{ID: }{max. 30 characters long,}
#'  \item{Sex: }{1 = females, 2 = males, other = unkown, except 4 =
#'    hermaphrodite,}
#' \item{BirthYear: }{(birth or hatching year) Integer, negative numbers are
#'    interpreted as missing values.}}
#' If the species has multiple generations per
#' year, use an integer coding such that the candidate parents' `Birth year' is
#' at least one smaller than their putative offspring's. Column names are
#' ignored, so ensure column order is ID - sex - birth year.
#' @param SeqList list with output from a previous run, containing the elements
#'   `Specs', `AgePriors' and/or `PedigreePar', as described below, to be used
#'   in the current run. If \code{SeqList$Specs} is provided, all other input
#'   parameter values except \code{MaxSibIter} are ignored.
#' @param MaxSibIter number of iterations of sibship clustering, including
#'   assignment of grandparents to sibships and avuncular relationships between
#'   sibships. Set to 0 to not (yet) perform this step, which is by far the most
#'   time consuming and may take several hours for large datasets. Clustering
#'   continues until convergence or until MaxSibIter is reached.
#' @param Err estimated genotyping error rate. The error model aims to deal with
#'   scoring errors typical for SNP arrays.
#' @param MaxMismatch maximum number of loci at which candidate parent and
#'   offspring are allowed to be opposite homozygotes. Setting a more liberal
#'   threshold can improve performance if the error rate is high, at the cost of
#'   decreased speed.
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
#' of erroneous assignments); used during full reconstruction only.
#' @param args.AP list with arguments to be passed on to
#'   \code{\link{MakeAgePrior}}.
#' @param FindMaybeRel  identify pairs of non-assigned likely relatives after
#'   pedigree reconstruction. Can be time-consuming in large datasets. NOTE:
#'   from v1.2 default changed from TRUE to FALSE; \code{\link{GetMaybeRel}} can
#'   now be called separately.
#' @param CalcLLR  calculate log-likelihood ratios for all assigned parents (is
#'   parent vs. is otherwise related). Time-consuming in large datasets.
#' @param quiet suppress messages: TRUE/FALSE/"verbose".
#'
#' @return A list with some or all of the following components:
#'   \item{AgePriors}{Matrix with age-difference based prior probability ratios,
#' used for full pedigree reconstruction. See \code{\link{MakeAgePrior}} for
#' details.}
#' \item{DummyIDs}{Dataframe with pedigree for dummy individuals, as well as
#' their sex, estimated birth year (point estimate, upper and lower bound of
#' 95\% confidence interval), number of offspring, and offspring IDs (genotyped
#' offspring only).}
#' \item{DupGenotype}{Dataframe, duplicated genotypes (with different IDs,
#'  duplicate IDs are not allowed). The specified number of maximum mismatches
#'   is used here too. Note that this dataframe may include pairs of closely
#'   related individuals, and monozygotic twins.}
#' \item{DupLifeHistID}{Dataframe, rownumbers of duplicated IDs in life
#'   history dataframe. For convenience only, but may signal a problem. The
#'   first entry is used.}
#' \item{ExcludedInd}{Individuals in GenoM which were excluded because of a
#'   too low genotyping success rate (<50\%).}
#' \item{ExcludedSNPs}{Column numbers of SNPs in GenoM which were excluded
#'   because of a too low genotyping success rate (<10\%).}
#' \item{LifeHist}{Provided dataframe with sex and birth year data.}
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
#'    homzygotes, but the offspring not being a heterozygote. The offspring
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
#' @seealso \code{\link{GenoConvert}, \link{SnpStats}, \link{GetMaybeRel},
#'   \link{EstConf}, \link{SummarySeq}, \link{writeSeq}}, vignette("sequoia")
#'
#' @examples
#' # == EXAMPLE 1 ==
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#' head(SimGeno_example[,1:10])
#' head(LH_HSg5)
#' SeqOUT <- sequoia(GenoM = SimGeno_example,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 0)
#' names(SeqOUT)
#' SeqOUT$PedigreePar[34:42, ]
#'
#' \dontrun{
#' SeqOUT2 <- sequoia(GenoM = SimGeno_example,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 10)
#' SeqOUT2$Pedigree[34:42, ]
#'
#' # == EXAMPLE 2 ==
#' # ideally, select 400-700 SNPs: high MAF & low LD
#' # save in 0/1/2/NA format (PLINK's --recodeA)
#' GenoM <- GenoConvert(InFile = "inputfile_for_sequoia.raw")
#' SNPSTATS <- SnpStats(GenoM)
#' # perhaps after some data-cleaning:
#' write.table(GenoM, file="MyGenoData.txt", row.names=T, col.names=F)
#' # later:
#' GenoM <- as.matrix(read.table("MyGenoData.txt", row.names=1, header=F))
#' LHdata <- read.table("LifeHistoryData.txt", header=T) # ID-Sex-birthyear
#' SeqOUT <- sequoia(GenoM, LHdata, Err=0.005, MaxMismatch=10)
#' SummarySeq(SeqOUT)
#' writeSeq(SeqOUT, OutFormat = "xls")  # needs library xlsx
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
                    UseAge = "yes",
                    args.AP = list(Flatten=TRUE, Smooth=TRUE),
                    FindMaybeRel = FALSE,
                    CalcLLR = TRUE,
                    quiet = FALSE)
{
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

  Excl <- CheckGeno(GenoM, quietF)

  if ("LifeHist" %in% names(SeqList)) {
    if(!quiet)  message("using LifeHistData in SeqList")
    LifeHistData <- SeqList$LifeHist
  }
  if (!is.null(LifeHistData)) {
    LifeHistData <- CheckLH(LifeHistData)
  } else {
    LifeHistData <- data.frame(ID = rownames(GenoM),
                              Sex = 3,
                              BY = -999,
                              stringsAsFactors = FALSE)
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
                     UseAge = as.character(SeqList$Specs["UseAge"]),
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
            				 UseAge = UseAge,
            				 FindMaybeRel = FindMaybeRel,
            				 CalcLLR = CalcLLR)
  }
  utils::flush.console()

  if ("AgePriors" %in% names(SeqList)) {
    if(!quiet)  message("using ageprior in SeqList")
    AgePriors <- SeqList$AgePriors
  } else {
    AgePriors <- do.call(MakeAgePrior, c(list(Ped = NULL,
                                              LifeHistData = LifeHistData,
                                              Plot = FALSE,
                                              quiet = TRUE),
                                         args.AP))
  }
  Specs[,"nAgeClasses"] <- nrow(AgePriors)

  if (MaxSibIter>-9) {
    DupList <- SeqParSib(ParSib = "dup", Specs = Specs, GenoM = GenoM,
		LhIN = LifeHistData[,1:3], quiet=quietF)

    if ("DupGenoID" %in% names(DupList)) {
      return(DupList)
    }
  } else DupList <- NULL
  utils::flush.console()

  if (any(LifeHistData$Sex==4)) {  # hermaphrodites - pretend 2 clones of opposite sex
    if (!quietF) message("detected hermaphrodites (sex=4), changing Complex to 'herm'")
    GenoM <- herm_clone_Geno(GenoM, LifeHistData, herm.suf=c("f", "m"))
    LifeHistData <- herm_clone_LH(LifeHistData, herm.suf=c("f", "m"))
    Specs[1, "NumberIndivGenotyped"] <- nrow(GenoM)
    Specs[1, "Complexity"] <- "herm"
  }

  if ("PedigreePar" %in% names(SeqList)) {
    PedParents <- SeqList$PedigreePar  #[, c("id", "dam", "sire")]
    for (x in 1:3)  PedParents[, x] <- as.character(PedParents[, x])
    for (x in 2:3)  PedParents[which(PedParents[, x]=="0"), x] <- NA
    if (!all(all(PedParents$id == rownames(GenoM))))  stop("SeqList$PedigreePar should match GenoM")
    DumPfx = paste0(DummyPrefix, "0")
#    if (any(substr(unlist(PedParent),1,nchar(DumPfx)) %in% DumPfx)) {
#      stop("(candidate) parent list cannot contain dummy IDs")
#    }
    if (any(LifeHistData$Sex==4)) {  # hermaphrodites
      PedParents <- herm_clone_Ped(Ped = PedParents, LH = LifeHistData[,1:3], herm.suf=c("f", "m"))
      rownames(PedParents) <- PedParents$id
      PedParents <- PedParents[rownames(GenoM), ]
    }

    if (MaxSibIter>0) {
      if(!quiet) message("using parents in SeqList")
      if (!"AgePriors" %in% names(SeqList) && !is.null(LifeHistData)) {
				AgePriors <- do.call(MakeAgePrior, c(list(Ped = PedParents[, 1:3],
                                              LifeHistData = LifeHistData,
                                              Plot = FALSE,
                                              quiet = TRUE),
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
      ParList <- SeqParSib(ParSib = "par", Specs = Specs, GenoM = GenoM,
                           LhIN = LifeHistData[,1:3], AgePriors = AgePriors,
                           Parents=PedParents, quiet=quietF)

    } # else keep old parents

  } else if (MaxSibIter>=0){
    ParList <- SeqParSib(ParSib="par", Specs=Specs, GenoM=GenoM, LhIN=LifeHistData[,1:3],
                         AgePriors=AgePriors, Parents=NULL, quiet=quietF)
    if (!"AgePriors" %in% names(SeqList) && !is.null(LifeHistData)) {
      AgePriors <- do.call(MakeAgePrior, c(list(Ped = ParList$PedigreePar[, 1:3],
                                                LifeHistData = LifeHistData,
                                                Plot = FALSE,
                                                quiet = TRUE),
                                         args.AP))
       if (nrow(AgePriors) != Specs[,"nAgeClasses"]) {
          Specs[,"nAgeClasses"] <- nrow(AgePriors)
        }
    }
  } else if (MaxSibIter <0) ParList <- NULL
  utils::flush.console()

  if (MaxSibIter>0) {
    SibList <- SeqParSib(ParSib = "sib", Specs = Specs, GenoM = GenoM,
                         LhIN = LifeHistData[,1:3], AgePriors = AgePriors,
                         Parents = ParList$PedigreePar, quiet = quietF)
    if (any(LifeHistData$Sex==4)) {
      ParList$PedigreePar <- herm_unclone_Ped(Ped = ParList$PedigreePar, LH = LifeHistData, herm.suf=c("f", "m"))
    }
  } else SibList <- NULL


  #=====================
  OUT <- list()
  OUT[["Specs"]] <- Specs
  if (length(Excl)>0)  OUT <- c(OUT, Excl)
  OUT[["AgePriors"]] <- AgePriors
  OUT[["LifeHist"]] <- LifeHistData
  OUT <- c(OUT, DupList, ParList)
  if (MaxSibIter>0)  OUT <- c(OUT, SibList)
  return(OUT)
}
