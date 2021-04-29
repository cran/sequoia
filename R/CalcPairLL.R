#' @title Calculate Likelihoods for Alternative Relationships
#'
#' @description For each specified pair of individuals, calculate the
#'   log10-likelihoods of being PO, FS, HS, GP, FA, HA, U (see Details).
#'   Individuals must be genotyped or have at least one genotyped offspring.
#'
#'   \strong{NOTE} values \eqn{>0} are various \code{NA} types, see 'Likelihood
#'   special codes' in 'Value' section below.
#'
#' @details The same pair may be included multiple times, e.g. with different
#'   sex, age difference, or focal relationship, to explore their effect on the
#'   likelihoods. Likelihoods are only calculated for relationships that are
#'   possible given the age difference, e.g. PO (parent-offspring) is not
#'   calculated for pairs with an age difference of 0.
#'
#'   Non-genotyped individuals can be included if they have at least one
#'   genotyped offspring and can be turned into a dummy (see
#'   \code{\link{getAssignCat}}); to establish this a pedigree must be provided.
#'
#'   \strong{Warning 1}: There is no check whether the input pedigree is genetically
#'   sensible, it is simply conditioned upon. Checking whether a pedigree is
#'   compatible with the SNP data can be done with \code{\link{CalcOHLLR}}.
#'
#'   \strong{Warning 2}: Conditioning on a \code{Pedigree} can make computation
#'   orders of magnitude slower.
#'
#'
#' @param Pairs  dataframe with columns \code{ID1} and \code{ID2}, and
#'   optionally
#'   \describe{
#'  \item{Sex1}{Sex of ID1, 1=female, 2=male, 3=unknown, or \code{NA} to take
#'    from \code{LifeHistData}. The sex of individuals occurring as parent in
#'    \code{Pedigree} cannot be altered.}
#'  \item{Sex2}{Sex of ID2}
#'  \item{AgeDif}{Age difference in whole time units, BirthYear1 - BirthYear2
#'    (i.e. positive if ID2 is born before ID1). If \code{NA}, calculated from
#'    \code{LifeHistData}. Use '999' to explicitly specify 'unknown'.}
#'  \item{focal}{relationship character abbreviation; PO, FS, HS, GP or U. See
#'  Details for its effect and explanation of abbreviations. Default: U}
#'  \item{patmat}{1=maternal relatives, 2=paternal relatives. Only relevant for
#'    HS & GP, for which it defaults to Sex1, or 1 if Sex1=3, but is currently
#'    only predictably implemented for pairs of two genotyped individuals.
#'    Always equal to Sex2 for PO pairs when Sex2 is known.}
#'  \item{dropPar1}{Drop the parents of \code{ID1} before calculating the pair
#'  likelihood, rather than conditioning on them; choose from 'none', 'dam',
#'  'sire', or 'both'. See example. If e.g. the pair shares a common mother,
#'  'none' and 'sire' will condition on this shared mother and not calculate the
#'  likelihood that they are maternal siblings, while dropPar1='dam' or 'both'
#'  will calculate that likelihood, and the other likelihoods as if the mother
#'  of ID1 were unknown.}
#'  \item{dropPar2}{as \code{dropPar1}, for \code{ID2}}
#'  }
#' @param Pedigree  dataframe with columns id-dam-sire; likelihoods will be
#'   calculated conditional on the pedigree. May include non-genotyped
#'   individuals, which will be treated as dummy individuals.
#' @param SeqList list with output from \code{\link{sequoia}}. If input
#'   parameter \code{Pedigree=NULL}, \code{SeqList$Pedigree} will be used if
#'   present, and \code{SeqList$PedigreePar} otherwise. If \code{SeqList$Specs}
#'   is present, input parameters with the same name as its items are ignored.
#'   The list elements 'LifeHist', 'AgePriors', and 'ErrM' are also used if
#'   present, and override the corresponding input parameters.
#' @param Plot logical, display scatter plots by \code{\link{PlotPairLL}}.
#' @inheritParams CalcOHLLR
#'
#'
#' @return The \code{Pairs} dataframe including all optional columns listed
#'   above, plus the additional columns:
#'  \item{LL_xx}{Log10-Likelihood of this pair having relationship xx, with xx
#'    being one of PO, FS, etc. as detailed below.}
#'  \item{TopRel}{Abbreviation of most likely relationship}
#'  \item{LLR}{Likelihood ratio between most-likely and second most likely
#'    relationships}
#'
#' \strong{Relationship abbreviations:}
#'   \item{PO}{Parent - offspring}
#'   \item{FS}{Full siblings}
#'   \item{HS}{Half siblings}
#'   \item{GP}{Grandparent}
#'   \item{FA}{Full avuncular}
#'   \item{HA}{Half avuncular and other 3rd degree relationships}
#'   \item{U}{Unrelated}
#'   \item{2nd}{Unclear which type of 2nd degree relatives
#'     (HS, GP, or FA)}
#'   \item{??}{Unclear which type of 1st, 2nd or 3rd degree
#'     relatives}
#'
#' \strong{Likelihood special codes:}
#'   \item{222}{Maybe (via) other parent (e.g. focal="GP", but as likely to be
#'     maternal as paternal grandparent, and therefore not assignable)}
#'   \item{333}{Excluded from comparison (shouldn't occur)}
#'   \item{444}{Not implemented (e.g. would create an odd double/triple
#'     relationship in combination with the provided pedigree)}
#'   \item{777}{Impossible (e.g. cannot be both full sibling and grandparent)}
#'   \item{888}{Already assigned in the provided pedigree (see \code{dropPar}
#'     arguments)}
#'   \item{999}{NA}
#'
#'
#' @section Double relationships & focal relationship:
#'   Especially when Complex='full', not only the seven relationship
#'   alternatives listed above are considered, but a whole range of possible
#'   double and even triple relationships. For example, mother A and offspring B
#'   (PO) may also be paternal half-siblings (HS, A and A's mother mated with
#'   same male), grandmother and grand-offspring (GP, B's father is A's son), or
#'   paternal aunt (B's father is a full or half sib of A).
#'
#'   The likelihood reported as 'LL_PO' is the most-likely one of the possible
#'   alternatives, among those that are not impossible due to age differences or
#'   due to the pedigree (as reconstructed up to that point). Whether e.g. the
#'   likelihood to be both PO & HS is counted as PO or as HS, depends on the
#'   situation and is determined by the variable 'focal': During parentage
#'   assignment, it is counted as PO but not HS, while during sibship
#'   clustering, it is counted as HS but not PO -- not omitting from the
#'   alternative relationship would result in a deadlock.
#'
#'
#' @seealso \code{\link{PlotPairLL}} to plot alternative relationship pairs from
#'   the output; \code{\link{CalcOHLLR}} to calculate LLR for parents &
#'   parent-pairs in a pedigree; \code{\link{GetRelM}} to find all pairwise
#'   relatives according to the pedigree; \code{\link{GetMaybeRel}} to get
#'   likely relative pairs not in the pedigree.
#'
#' @examples
#' # likelihoods underlying parent LLR in pedigree:
#' data(LH_HSg5, SimGeno_example)
#' Seq.HSg5 <- sequoia(SimGeno_example, LH_HSg5, Module="par")
#' tail(Seq.HSg5$PedigreePar)
#' # take bottom 3 individuals:
#' Pairs <- data.frame(ID1 = c("a01190", "b01191", "a01192"),
#'                     ID2 = rep(c("a00012", "b00007"), each=3),
#'                     AgeDif = 1, focal = "PO")
#' # LLRdam & LLRsire:
#' CalcPairLL(Pairs, SimGeno_example)
#' # LLRpair is min. of dam & sire LLR, conditional on co-parent:
#' CalcPairLL(cbind(Pairs, dropPar1=rep(c("dam", "sire"), each=3)),
#'            SimGeno_example, SeqList = Seq.HSg5)
#'
#' \donttest{
#' # likelihoods underlying LLR in getMaybeRel output:
#' data(Ped_griffin, SeqOUT_griffin, package="sequoia")
#' Geno.griffin <- SimGeno(Ped_griffin, nSnp=200, SnpError = 0.01, ParMis=0.4)
#' MR <- GetMaybeRel(GenoM = Geno.griffin,
#'                   LifeHistData = SeqOUT_griffin$LifeHist,
#'                   Module = "par", Err = 0.001)
#' FivePairs <- MR$MaybePar[1:5, c("ID1", "ID2", "Sex1", "Sex2")]
#' FivePairs$AgeDif <- NA   # pretend unknown age differences
#' PairLL <- CalcPairLL(Pairs = rbind( cbind(FivePairs, focal = "PO"),
#'                                     cbind(FivePairs, focal = "HS"),
#'                                     cbind(FivePairs, focal = "GP")),
#'                      GenoM = Geno.griffin,
#'                      Err = 0.005, Plot=FALSE)
#' PairLL[c(1, 6, 11), ]
#' # LL(FS)==222 : HSHA, HSGP, FAHA more likely than FS
#' # LL(GP) higher when focal=HS: GP via 'other' parent also considered
#' # LL(FA) higher when focal=PO: FAHA, or FS of 'other' parent
#' }
#'
#' @useDynLib sequoia, .registration = TRUE
#
#' @export

CalcPairLL <- function(Pairs = NULL,
                       GenoM = NULL,
                       Pedigree = NULL,
                       LifeHistData = NULL,
                       AgePrior = TRUE,
                       SeqList = NULL,
                       Complex = "full",
                       Herm = "no",
                       Err = 1e-4,
                       ErrFlavour = "version2.0",
                       Tassign = 0.5,
                       Tfilter = -2.0,
                       quiet = FALSE,
                       Plot = TRUE)
{
  on.exit(.Fortran(deallocall), add=TRUE)

  # check genotype data ----
  GenoM <- CheckGeno(GenoM, quiet=TRUE, Plot=FALSE)
  gID <- rownames(GenoM)

  # unpack SeqList ----
  if (!is.null(SeqList)) {
    NewName = c(Pedigree = "Pedigree",
                PedigreePar = "Pedigree",
                LifeHist = "LifeHistData",
                AgePriors = "AgePrior")
    for (x in names(NewName)) {
      if (x %in% c("Pedigree", "PedigreePar") & !is.null(Pedigree))  next
      if (x %in% names(SeqList)) {
        if (x == "PedigreePar" & "Pedigree" %in% names(SeqList))  next
        if (x == "AgePriors" & AgePrior == FALSE)  next
        if (!quiet)  message("using ", x, " in SeqList")
        assign(NewName[x], SeqList[[x]])
      }
    }
  }

  Ped <- PedPolish(Pedigree, gID, DropNonSNPd=FALSE, NullOK = TRUE)
  # Ped=NULL if Pedigree=NULL & gID=NULL

  LHF <- CheckLH(LifeHistData, rownames(GenoM), sorted=TRUE)
  LHF$Sex[LHF$Sex==4] <- 3



  # turn IDs into numbers, & turn non-genotyped parents into temporary dummies ----
  PedN <- PedToNum(Ped, gID, DoDummies = "new")  # list: PedPar - DumPar - Renamed - Nd
  PairL <- FortifyPairs(Pairs, gID, PedN$Renamed, LHF)  # also checks input


  # check/make ageprior ----
  if (is.logical(AgePrior) && AgePrior==TRUE) {
    if (is.null(Pedigree) | is.null(LifeHistData))
      AgePrior <- FALSE
  }

  MaxAgeDif <- max(c(0,PairL$AgeDif[PairL$AgeDif<999]), na.rm=TRUE)
  if (is.logical(AgePrior) && AgePrior %in% c(TRUE, FALSE)) {
    if (AgePrior) {
      AP <- MakeAgePrior(Ped, LifeHistData,
                         MaxAgeParent = max(c(99, MaxAgeDif, na.rm=TRUE)),
                         quiet = quiet, Plot = Plot)
    } else {
      AP <- MakeAgePrior(quiet = TRUE, Plot = FALSE)    # non-informative ageprior
    }
  } else if (is.matrix(AgePrior) | is.data.frame(AgePrior)) {
    AP <- CheckAP(AgePrior)
    if (nrow(AP) < MaxAgeDif) {
      APz <- matrix(0, MaxAgeDif - nrow(AP), 5,
                    dimnames=list(
                      seq.int(from=901, length.out = MaxAgeDif - nrow(AP)),
                      colnames(AP)))
      AP <- rbind(AP, APz)
    }
  } else {
    stop("'AgePrior' must be TRUE or FALSE, or a matrix")
  }

  # Specs / param ----
  if ("Specs" %in% names(SeqList)) {
    if(!quiet)  message("settings in SeqList$Specs will overrule input parameters")
    PARAM <- SpecsToParam(SeqList$Specs, SeqList$ErrM, ErrFlavour,
                          dimGeno = dim(GenoM), quiet, Plot)  # overrule values in SeqList
    PARAM$nAgeClasses <- nrow(AP)
  } else {
    PARAM <- namedlist(dimGeno = dim(GenoM),
#                       dropPar,
                       Err,
                       ErrFlavour,
                       Tfilter,
                       Tassign,
                       nAgeClasses = nrow(AP),
                       MaxSibshipSize = max(table(Ped$dam), table(Ped$sire), 90,
                                            na.rm=TRUE) +10,
                       Complex,
                       Herm,
                       quiet,
                       Plot)
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


  # check parameter values ----
  CheckParams(PARAM)
  utils::flush.console()


  # call Fortran ----
  Np <- nrow(Pairs)
  FortPARAM <- MkFortParams(PARAM, fun="CalcPairs")

  TMP <- .Fortran(getpairll,
                  ng = as.integer(nrow(GenoM)),   # no. genotyped indiv
                  np = as.integer(Np),   # no. pairs
                  specsint = as.integer(FortPARAM$SpecsInt),
                  specsdbl = as.double(FortPARAM$SpecsDbl),
                  errv = as.double(FortPARAM$ErrM),
                  genofr = as.integer(GenoM),

                  byrf = as.integer(c(LHF$BirthYear, LHF$BY.min, LHF$BY.max)),
                  aprf = as.double(AP),

                  pairids = as.integer(PairL$ID),
                  pairsex = as.integer(PairL$Sex),
                  pairagediff = as.integer(PairL$AgeDif),
                  pairfocal = as.integer(PairL$focal),
                  pairk = as.integer(PairL$patmat),
                  dropp = as.integer(PairL$dropPar),

                  parentsrf = as.integer(PedN$PedPar),
                  dumparrf = as.integer(PedN$DumPar),
                  llrf = double(7*Np),
                  toprf = integer(Np),
                  dlrf = double(Np)
  )

  # wrap output ----
  RelNames <- c("PO", "FS", "HS", "GP", "FA", "HA", "U", "??", "2nd")
  Pairs.OUT <- with(PairL, data.frame(Pairs[,1:2],
                                      Sex1 = Sex[1:Np],
                                      Sex2 = Sex[(Np+1):(2*Np)],
                                      AgeDif,
                                      focal = RelNames[PairL$focal],
                                      patmat,
                                      dropPar1 = c("none", "dam", "sire", "both")[dropPar[1:Np] +1],
                                      dropPar2 = c("none", "dam", "sire", "both")[dropPar[(Np+1):(2*Np)] +1],
                                      stringsAsFactors = FALSE))
  Pairs.OUT$AgeDif[Pairs.OUT$AgeDif == 999] <- NA

  Pairs.OUT <- cbind(Pairs.OUT,
                     setNames(as.data.frame(round(VtoM(TMP$llrf, nc=7), 2)),
                              RelNames[1:7]),
                     TopRel = RelNames[TMP$toprf],
                     LLR = round(TMP$dlrf, 2))
  Pairs.OUT$LLR[Pairs.OUT$LLR == -777] <- NA

  # plot ----
  if (Plot) {
    if (nrow(Pairs) > 1e4) {
      message("Plot not shown when >10.000 pairs")
    } else {
      PlotPairLL(Pairs.OUT)
    }
  }

  return( Pairs.OUT )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Make Pairs Fortran Compatible
#'
#' @description Convert dataframe \code{Pairs} into a list of integer vectors.
#'   Called only by \code{\link{CalcPairLL}}.
#'
#' @param Pairs  dataframe with columns ID1 - ID2 - Sex1 - Sex2 - AgeDif - focal
#'   - k.
#' @param gID  character vector with IDs of genotyped individuals.
#' @param Renamed  length-2 list (dams, sires) each with a 2-column dataframe.
#'   matching character IDs to negative numbers, for dummified individuals.
#'   Element of the list returned by \code{\link{PedToNum}}.
#' @param LH  lifehistory dataframe, ID - Sex - BirthYear.
#'
#' @return A named list, with elements ID - Sex - AgeDif - focal. The first two
#'   are per individual and thus each have length 2*nrow(Pairs), while the last
#'   two have length 1*nrow(Pairs).
#'
#' @keywords internal

FortifyPairs <- function(Pairs,   # pairs with character IDs etc
                         gID,
                         Renamed,
                         LH)
{
  Pairs <- as.data.frame(Pairs)   # in case it's a matrix
  if (ncol(Pairs)==2) {
    colnames(Pairs) <- c("ID1", "ID2")
  } else if (!all(c("ID1", "ID2") %in% names(Pairs))) {
    stop("Pairs must  have at least columns 'ID1' and 'ID2'")
  }
  for (x in c("ID1", "ID2"))   Pairs[,x] <- as.character(Pairs[,x])

  names(Pairs)[-c(1:2)] <- tolower(names(Pairs)[-c(1:2)])
  if ("agediff" %in% names(Pairs)) {
    names(Pairs)[names(Pairs) == "agediff"] <- "agedif"
    #    warning("")
  }
  PairsColNames <- c(sex1 = "Sex1", sex2 = "Sex2", agedif = "AgeDif",
                     focal = "focal", patmat = "patmat",
                     droppar1 = "dropPar1", droppar2="dropPar2")
  for (x in names(PairsColNames)) {
    if (!x %in% names(Pairs))  Pairs[,x] <- NA   # add column if not already present
    names(Pairs)[names(Pairs) == x] <- PairsColNames[x]
  }


  # add LH data ----
  Pairs$zz <- seq_along(Pairs$ID1)   # to fix merge() row order mess
  if (any(is.na(Pairs$AgeDif)) | any(is.na(Pairs$Sex1)) | any(is.na(Pairs$Sex2))) {
    Pairs <- merge(Pairs, setNames(LH[, 1:3], c("ID1", "Sex1.LH", "BY1")), all.x=TRUE)
    Pairs <- merge(Pairs, setNames(LH[, 1:3], c("ID2", "Sex2.LH", "BY2")), all.x=TRUE)
  }
  Pairs <- Pairs[order(Pairs$zz), ]

  # age difference ----
  if (any(!is.wholenumber(Pairs$AgeDif) & !is.na(Pairs$AgeDif)))
    stop("'AgeDif' in 'Pairs' must be whole numbers", call.=FALSE)
  if (any(is.na(Pairs$AgeDif))) {
    Pairs$AgeDif <- with(Pairs, ifelse(!is.na(AgeDif),
                                       AgeDif,
                                       ifelse(BY1 >=0 & BY2 >= 0,
                                              BY1 - BY2,
                                              999)))
    Pairs$AgeDif[is.na(Pairs$AgeDif)] <- 999
  }

  # sex ----
  for (x in c("Sex1", "Sex2")) {
    if (!all(Pairs[,x] %in% c(1:3, NA)))
      stop("'Sex1' and 'Sex2' in 'Pairs' must be 1=female, 2=male, 3/NA= unknown", call.=FALSE)
    if (any(is.na(Pairs[,x]))) {
      Pairs[,x] <- ifelse(!is.na(Pairs[,x]),
                          Pairs[,x],
                          Pairs[,paste0(x, ".LH")])
    }
    Pairs[is.na(Pairs[,x]), x] <- 3
  }

  # focal ----
  Pairs$focal[is.na(Pairs$focal)] <- "U"
  RelNames <- c("PO", "FS", "HS", "GP", "FA", "HA", "U")
  if (!all(Pairs$focal %in% c(RelNames, NA))) {
    stop("Some 'focal' values in 'Pairs' are invalid. ",
         "Valid values are: PO, FS, HS, GP, FA, HA, U", call.=FALSE)
  }
  Pairs$focal[is.na(Pairs$focal)] <- "U"


  # dropPar ----
  for (x in c("dropPar1", "dropPar2")) {
    Pairs[,x] <- as.character(Pairs[,x])
    if (!all(Pairs[,x] %in% c("none", "dam", "sire", "both", NA)))
      stop("'dropPar1' and 'dropPar2' in 'Pairs' must be 'none', 'dam', 'sire', or 'both'", call.=FALSE)
    Pairs[is.na(Pairs[,x]), x] <- "none"
  }


  # character IDs to numbers ----
  ID <- c("ID1", "ID2")
  SEX <- c("Sex1", "Sex2")
  GenoNums <- setNames(seq_along(gID), gID)

  for (x in 1:2) {
    if (!all(Pairs[,ID[x]] %in% c(gID, unlist(Renamed))))
      stop("All individuals must be genotyped or dummifiable", call.=FALSE)
    Pairs[, paste0(ID[x], ".num")] <- ifelse(Pairs[, ID[x]] %in% gID,
                                             GenoNums[Pairs[, ID[x]]],
                                             0)
    if (all(Pairs[, paste0(ID[x], ".num")] > 0))  next   # no dummies
    for (k in 1:2) {
      dum.k <- Pairs[,ID[x]] %in% Renamed[[k]]$name
      if (any(Pairs[dum.k, SEX[x]] == 3-k))
        stop("Some ", c("dams", "sires")[k], " in Pedigree have Sex = ",
             3-k, " in Pairs or LifeHistData", call.=FALSE)
      Pairs[dum.k, SEX[x]] <- k   # fixes any with unknown sex
      Pairs[dum.k, paste0(ID[x], ".num")] <- Renamed[[k]][match(Pairs[dum.k, ID[x]],
                                                                Renamed[[k]]$name),
                                                          "num"]
    }
  }
  if (any(Pairs$ID1.num == 0 | Pairs$ID2.num == 0))  stop("Something went wrong")


  # paternal/maternal ----
  if (!all(Pairs$patmat %in% c(1,2,NA)))
    stop("'patmat' in 'Pairs' must be 1, 2, or NA")
  Pairs$patmat <- with(Pairs, ifelse(!is.na(patmat), patmat,
                                     ifelse(focal == "PO" & Sex2 != 3, Sex2,
                                            ifelse(Sex1 != 3, Sex1, 1))))


  # output list, for Fortran ----
  PairL <- with(Pairs, list(ID = c(ID1.num, ID2.num),
                            Sex = c(Sex1, Sex2),
                            AgeDif = AgeDif,
                            focal = as.numeric(factor(as.character(focal),
                                                      levels = RelNames)),
                            patmat = patmat,
                            dropPar = FacToNum(factor(c(dropPar1, dropPar2),
                                                      levels=c("none", "dam", "sire", "both"),
                                                      labels=c(0,1,2,3)))
  ))

  return( PairL )
}


