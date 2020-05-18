#' @title Find putative relatives
#'
#' @description Identify pairs of individuals likely to be related, but not
#'   assigned as such in the provided pedigree.
#'
#' @param GenoM matrix with genotype data, size nInd x nSnp
#' @param SeqList list with output from \code{\link{sequoia}}. If provided, the
#' elements `Specs', `AgePriors' and 'LifeHist' are used, and all other input
#' parameters except 'GenoM', 'ParSib' and 'quiet' are ignored.
#' @param Pedigree dataframe with id - dam - sire in columns 1-3
#' @param LifeHistData dataframe with columns id - sex (1=female, 2=male,
#'   3=unknown) - birth year
#' @param ParSib either 'par' to check for putative parent-offspring pairs only,
#'   or 'sib' to check for all types of first and second degree relatives. When
#'   'par', all pairs are returned that are more likely parent-offspring than
#'   unrelated, including pairs that are even more likely to be otherwise
#'   related.
#' @param Complex either "full" (default), "simp" (simplified, no explicit
#'   consideration of inbred relationships), "mono" (monogamous) or "herm"
#'   (hermaphrodites, otherwise like "full").
#' @param Err estimated genotyping error rate, as a single number or 3x3 matrix.
#'   If a matrix, this should be the probability of observed genotype (columns)
#'   conditional on actual genotype (rows). Each row must therefore sum to 1.
#' @param ErrFlavour function that takes \code{Err} as input, and returns a 3x3
#'   matrix of observed (columns) conditional on actual (rows) genotypes, or
#'   choose from inbuilt ones as used in sequoia 'version2.0', 'version1.3', or
#'   'version1.1'. Ignored if \code{Err} is a matrix. See \code{\link{ErrToM}}.
#' @param MaxMismatch DEPRECATED AND IGNORED. Now calculated using
#'   \code{\link{CalcMaxMismatch}}.
#' @param Tassign minimum LLR required for acceptance of proposed relationship,
#'   relative to next most likely relationship. Higher values result in more
#'   conservative assignments. Must be zero or positive.
#' @param MaxPairs  The maximum number of putative pairs to return.
#' @param DumPrefix  character vector of length 2 with prefixes for dummy dams
#'   (mothers) and sires (fathers) used in \code{Pedigree}.
#' @param quiet suppress messages
#'
#' @return A list with
#'   \item{MaybeParent or MaybeRel}{A dataframe with non-assigned likely
#'     relatives, with columns ID1 - ID2 - TopRel - LLR - OH - BirthYear1 -
#'     BirthYear2 - AgeDif - Sex1 - Sex2 - SNPdBoth}
#'   \item{MaybeTrio}{A dataframe with non-assigned parent-parent-offspring
#'     trios, with columns id - parent1 - parent2 - LLRparent1 - LLRparent2 -
#'     LLRpair - OHparent1 - OHparent2 - MEpair - SNPd.id.parent1 -
#'     SNPd.id.parent2}
#' The following categories are used in column 'TopRel', indicating the most
#' likely relationship category:
#' \item{PO}{Parent-Offspring}
#' \item{FS}{Full Siblings}
#' \item{HS}{Half Siblings}
#' \item{GP}{GrandParent - grand-offspring}
#' \item{FA}{Full Avuncular (aunt/uncle)}
#' \item{2nd}{2nd degree relatives, not enough information to distinguish between
#'   HS,GP and FA}
#' \item{Q}{Unclear, but probably 1st, 2nd or 3rd degree relatives}
#'
#' @examples
#' \dontrun{
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#' SeqOUT <- sequoia(GenoM = SimGeno_example,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 0)
#' MaybePO <- GetMaybeRel(GenoM = SimGeno_example,
#'                       SeqList = SeqOUT)
#'
#' Maybe <- GetMaybeRel(GenoM = SimGeno_example,
#'                       Pedigree = SeqOUT$PedigreePar, ParSib="sib")
#' }
#'
#' @useDynLib sequoia, .registration = TRUE
#'
#' @export

GetMaybeRel <- function(GenoM = NULL,
                        SeqList = NULL,
                        Pedigree = NULL,
                        LifeHistData = NULL,
                        ParSib = "par",
                        Complex = "full",
                        Err = 0.0001,
                        ErrFlavour = "version2.0",
                        MaxMismatch = NA,   # DEPRECATED
                        Tassign = 0.5,
                        MaxPairs =  7*nrow(GenoM),
                        DumPrefix = c("F0", "M0"),
                        quiet = FALSE)
{
  on.exit(.Fortran(deallocall), add=TRUE)

  # input check
	Excl <- CheckGeno(GenoM, quiet=quiet)
	if (ParSib == "full") ParSib <- "sib"
  if (!ParSib %in% c("par", "sib"))  stop("Invalid value for 'ParSib', choose 'par' or 'sib'")
	if (nchar(DumPrefix[1]) != nchar(DumPrefix[2])) stop("DumPrefix must have same number of characters")

  if (!is.null(SeqList) && (!is.list(SeqList) | is.data.frame(SeqList))) {
      stop("SeqList must be a list or NULL")
  }
  if (!is.na(MaxMismatch)) {
    warning("NOTE: 'MaxMismatch' is deprecated & ignored; now calculated internally by 'CalcMaxMismatch'",
            immediate.=TRUE)
  }
  if ("Pedigree" %in% names(SeqList)) {
    Pedigree <- SeqList$Pedigree
    if (!all(rownames(GenoM) %in% Pedigree$id))   stop("SeqList Pedigree does not match GenoM")
  } else if ("PedigreePar" %in% names(SeqList)) {
    Pedigree <- SeqList$PedigreePar
    if (!all(Pedigree$id == rownames(GenoM)))  stop("SeqList PedigreePar does not match GenoM")
  } else if (!is.null(Pedigree)) {
    if (!all(rownames(GenoM) %in% Pedigree[,1]))  stop("Pedigree must include all individuals in GenoM, or be NULL")
  }

  PrSb <- switch(ParSib, dup=0, par = 1, sib = 2)
  if (is.null(MaxPairs) || MaxPairs<0 || !is.wholenumber(MaxPairs)) {
    stop("'MaxPairs' but be a positive whole number")
  }

  if ("LifeHist" %in% names(SeqList)) {
    LhIN <- SeqList$LifeHist
  } else if (!is.null(LifeHistData)) {
    LhIN <- CheckLH(LifeHistData)
  } else {
    LhIN <- orderLH(LH=NULL, gID = rownames(GenoM))
  }

	if (any(LhIN$Sex==4)) {  # hermaphrodites - pretend 2 clones of opposite sex
    if (!quiet) message("detected hermaphrodites (sex=4), changing Complex to 'herm'")
    GenoM <- herm_clone_Geno(GenoM, LhIN, herm.suf=c("f", "m"))
    if (!"LifeHist" %in% names(SeqList)) {
      LhIN <- herm_clone_LH(LhIN, herm.suf=c("f", "m"))
    }  # else already cloned
    if (!is.null(Pedigree)) {
      Pedigree <- herm_clone_Ped(Ped = Pedigree, LH = LhIN[,1:3], herm.suf=c("f", "m"))
      Pedigree <- Pedigree[match(rownames(GenoM), Pedigree[,1]), ]
    }
	}

	gID <- rownames(GenoM)
  GenoV <- as.integer(GenoM)

  if (!is.null(Pedigree)) {
    Ped <- PedPolish(Pedigree, GenoNames = rownames(GenoM), NAToZero = TRUE)
    DPnc <- nchar(DumPrefix)[1]
    PedNum <- IDToNum(Ped[,1:3], gID, DumPrefix)
    PedPar <- as.matrix(PedNum[gID, 2:3])   # not dummy indivs
    PedPar[is.na(PedPar)] <- 0   # in case not all genotyped indivs in Ped
  } else {
    PedPar <- rep(0, 2*nrow(GenoM))
  }

  if ("Specs" %in% names(SeqList)) {
    Specs <- SeqList$Specs
    Ng <- FacToNum(Specs[,"NumberIndivGenotyped"])
    SMax <- FacToNum(Specs[,"MaxSibshipSize"])
		ErrM <- ErrToM(FacToNum(SeqList$Specs["GenotypingErrorRate"]),
                     flavour = ErrFlavour, Return = "matrix")
    if (!"MaxMismatchOH" %in% names(SeqList$Specs)) {  # backwards compatability, version 1.x
			sts <- SnpStats(GenoM, Plot=FALSE)
			MaxMismatchV <- CalcMaxMismatch(Err=ErrM, MAF=sts[,"AF"], ErrFlavour=ErrFlavour,
			                                qntl=0.999^(1/Ng))  # length 3 vector
    } else {
      MaxMismatchV <- setNames(FacToNum(SeqList$Specs[c("MaxMismatchDUP", "MaxMismatchOH",
                                                        "MaxMismatchME")]),
                               c("DUP", "OH", "ME"))
    }

    Cmplx <- switch(Specs[,"Complexity"], full = 2, simp = 1, mono = 0, herm = 4)
    AP <- SeqList$AgePriors[, c("M", "P", "FS", "MS", "PS")]

    SpecsInt <- c(nSnp = FacToNum(Specs[,"NumberSnps"]),           # 1
                  MaxMisDUP = FacToNum(Specs[,"MaxMismatchDUP"]),  # 2
                  MaxMisOH = FacToNum(Specs[,"MaxMismatchOH"]),    # 3
                  MaxMisME = FacToNum(Specs[,"MaxMismatchME"]),    # 4
                  SMax = as.integer(FacToNum(Specs[,"MaxSibshipSize"])),   # 5
                  Complx = as.integer(Cmplx),                 # 6
                  quiet = as.integer(quiet),                   # 7
                  nAgeCl = FacToNum(Specs[,"nAgeClasses"]))    # 8
    SpecsDbl <- c(TF = FacToNum(Specs[,"Tfilter"]),
                  TA = FacToNum(Specs[,"Tassign"]))

  } else {
    Ng <- nrow(GenoM)
    if (is.null(Pedigree)) {
      SMax <- 100
    } else {
      SMax <- max(table(Pedigree$dam), table(Pedigree$sire)) +1
    }
    Complx <- switch(Complex, full = 2, simp = 1, mono = 0, herm = 4)
    AP <- MakeAgePrior(Pedigree, LifeHistData, Plot=FALSE, quiet=TRUE)

    ErrM <- ErrToM(Err, flavour = ErrFlavour, Return = "matrix")  # performs checks if Err is already matrix

    sts <- SnpStats(GenoM, Plot=FALSE)
    MaxMismatchV <- CalcMaxMismatch(Err=ErrM, MAF=sts[,"AF"], ErrFlavour=ErrFlavour,
                                    qntl=0.999^(1/Ng))

    SpecsInt <- c(nSnp = as.integer(ncol(GenoM)),     # 1
                  MaxMisDUP = as.integer(MaxMismatchV["DUP"]),  # 2
                  MaxMisOH = as.integer(MaxMismatchV["OH"]),    # 3
                  MaxMisME = as.integer(MaxMismatchV["ME"]),    # 4
                  SMax = as.integer(SMax),            # 5
                  Complx = as.integer(Complx),        # 6
                  quiet = as.integer(quiet),          # 7
                  nAgeCl = as.integer(nrow(AP)))      # 8
    SpecsDbl <- c(TF = as.double(-2.0),  # not used
                  TA = as.double(Tassign))
  }
  SpecsIntAmb <- c(ParSib = as.integer(PrSb),
                   nAmbMax = as.integer(MaxPairs))




	# Dummies
  Nd <- 0
  DumParRF <- rep(0, 4*as.integer(Ng/2))
  dID <- NULL
  if (!is.null(Pedigree)) {
    Nd <- c(sum(substr(Pedigree$id,1,DPnc)==DumPrefix[1]),
            sum(substr(Pedigree$id,1,DPnc)==DumPrefix[2]))
    if (max(Nd)>0) {
      SibshipGPs <- array(0, dim=c(2,max(Nd),2),
                          dimnames=list(c("grandma", "granddad"), 1:max(Nd), c("mat", "pat")))
      for (k in 1:2) {
        if (Nd[k]>0) {
          SibshipGPs[,1:Nd[k],k] <- t(as.matrix(PedNum[substr(Pedigree$id,1,DPnc)==DumPrefix[k], 2:3]))
          for (s in 1:Nd[k]) {
            for (g in 1:2) {
              x <- (k-1)*2*as.integer(Ng/2) + (s-1)*2 + g
              DumParRF[x] <- SibshipGPs[g,s,k]
            }
          }
        }
      }
      dID <- c(Pedigree$id[substr(Pedigree$id,1,DPnc)==DumPrefix[1]],
               Pedigree$id[substr(Pedigree$id,1,DPnc)==DumPrefix[2]])
    }
  }

  LHF <- orderLH(LhIN[LhIN$Sex %in% c(1:3), ], gID)
  #=========================
	# call fortran
  TMP <- .Fortran(findambig,
                  ng = as.integer(Ng),
                  specsint = as.integer(SpecsInt),
                  specsintamb = as.integer(SpecsIntAmb),
                  specsdbl = as.double(SpecsDbl),
									errv = as.double(ErrM),
                  genofr = as.integer(GenoV),
									sexrf = as.integer(LHF$Sex),
                  byrf = as.integer(c(LHF$BirthYear, LHF$BY.min, LHF$BY.max)),
									aprf = as.double(AP),
									parentsrf = as.integer(PedPar),
                  dumparrf = as.integer(DumParRF),

									namb = as.integer(0),
                  ambigid = integer(2*MaxPairs),
                  ambigrel = integer(2*MaxPairs),
                  ambiglr = double(2*MaxPairs),
                  ambigoh = integer(MaxPairs),

									ntrio = as.integer(0),
                  trioids = integer(3*Ng),
                  triolr = double(3*Ng),
                  triooh = integer(3*Ng))

	TMP$ambiglr[abs(TMP$ambiglr - 999) < 0.1] <- NA
  TMP$ambiglr <- round(TMP$ambiglr, 2)
  TMP$ambigoh[TMP$ambigoh < 0] <- NA
  TMP$triolr[abs(TMP$triolr - 999) < 0.1] <- NA
  TMP$triolr <- round(TMP$triolr, 2)
  TMP$triooh[TMP$triooh < 0] <- NA

  #=========================
  # non-assigned probable relative pairs

  if (TMP$namb > 0) {
    RelName <- c("PO", "FS", "HS", "GP", "FA", "HA", "U ", "Q", "2nd")
    Na <- TMP$namb
    TMP$ambigid <- NumToID(TMP$ambigid, 0, gID, dID)
    AmbigRel <- factor(TMP$ambigrel, levels=1:9, labels=RelName)

    MaybeRel <- data.frame(VtoM(TMP$ambigid, Na),
                           VtoM(AmbigRel, Na),
                           VtoM(TMP$ambiglr, Na),
                           stringsAsFactors=FALSE)
    names(MaybeRel) <- c("ID1", "ID2", "Relx", "TopRel", "LLR_Rx_U", "LLR")
    MaybeRel <- MaybeRel[,-which(names(MaybeRel) %in% c("Relx", "LLR_Rx_U"))]  # drop; confusing.
    MaybeRel$OH <-  TMP$ambigoh[s(Na)]

    if (!is.null(LhIN) & nrow(MaybeRel)>0) {
      LhIN$BirthYear[LhIN$BirthYear<0] <- NA
      MaybeRel <- merge(MaybeRel, setNames(LhIN[,1:3], c("ID1","Sex1","BirthYear1")), all.x=TRUE)
      MaybeRel <- merge(MaybeRel, setNames(LhIN[,1:3], c("ID2","Sex2","BirthYear2")), all.x=TRUE)
      MaybeRel$AgeDif <- with(MaybeRel, BirthYear1 - BirthYear2)
      MaybeRel <- MaybeRel[, c("ID1", "ID2", "TopRel", "LLR", "OH",
                               "BirthYear1", "BirthYear2", "AgeDif", "Sex1", "Sex2")]
      for (i in 1:Na) {
        if (is.na(MaybeRel$AgeDif[i]))  next
        if (MaybeRel$AgeDif[i] < 0) {
          tmpRel <- MaybeRel[i,]
          tmpRel$AgeDif <- abs(tmpRel$AgeDif)
          MaybeRel[i,] <- tmpRel[,c("ID2","ID1","TopRel", "LLR", "OH",
                                 "BirthYear2", "BirthYear1","AgeDif", "Sex2", "Sex1")]
        }
      }
    }

    if (grepl("sib", ParSib)) {
      MaybeRel <- with(MaybeRel, MaybeRel[TopRel %in% c("PO", "FS", "HS","GG", "FA", "2nd", "Q"), ])  # drop HA,XX,U
    }
    MaybeRel <- MaybeRel[order(ordered(MaybeRel$TopRel, levels=RelName),
                                     -MaybeRel$LLR),]

    if (nrow(MaybeRel)==0) {
      MaybeRel <- NULL
    } else {
      rownames(MaybeRel) <- 1:nrow(MaybeRel)
      MaybeRel$SNPdBoth <- CalcSnpdBoth(MaybeRel[, c("ID1", "ID2")], GenoM)

      if (any(LhIN$Sex==4)) {  # hermaphrodites
        MaybeRel <- herm_unclone_MaybeRel(MaybeRel, Ped = NULL, LH=LhIN, herm.suf=c("f", "m"))
      }
    }
  } else  MaybeRel <- NULL

  if (quiet<1) {   #  && nrow(MaybeRel)>0
      if (!is.null(MaybeRel)) {
        nRel <- nrow(MaybeRel)
        nPO <- sum(MaybeRel$TopRel=="PO")
      } else {
        nRel <- 0
        nPO <- 0
      }
      message("Found ", nPO, " likely parent-offspring pairs, and ",
              nRel - nPO, " other non-assigned pairs of possible relatives")
  }



  #=========================
  # non-assigned parent-parent-offspring trios
  if (TMP$ntrio>0) {
    trios <- data.frame(VtoM(TMP$trioids, nr=TMP$ntrio, nc=3),
                         VtoM(TMP$triolr, nr=TMP$ntrio, nc=3),
                        VtoM(TMP$triooh, nr=TMP$ntrio, nc=3),
                         stringsAsFactors=FALSE)
    names(trios) <- c("id", "parent1", "parent2", "LLRparent1", "LLRparent2", "LLRpair",
                      "OHparent1", "OHparent2", "MEpair")
    for (k in 1:3) trios[, k] <- NumToID(trios[, k], k-1, gID, dID)
    trios$SNPd.id.parent1 <- CalcSnpdBoth(trios[, c("id", "parent1")], GenoM)
    trios$SNPd.id.parent2 <- CalcSnpdBoth(trios[, c("id", "parent2")], GenoM)

    if (any(LhIN$Sex==4)) {  # hermaphrodites
      trios <- herm_unclone_Trios(trios, LH=LhIN, herm.suf=c("f", "m"))
    }

    if (quiet<1) {
      message("Found ", nrow(trios), " parent-parent-offspring trios")
    }
  } else  trios <- NULL

  if (grepl("par", ParSib)) {
    return( list(MaybePar = MaybeRel,
               MaybeTrio = trios) )
  } else {
    return( list(MaybeRel = MaybeRel,
               MaybeTrio = trios) )
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title  SNPs scored in both individuals
#'
#' @description count no. SNPs for which both individuals of pair are
#'   successfully SNPd
#'
#' @param Pairs  dataframe or matrix with 2 columns with IDs
#' @param GenoM matrix with genotype data, individuals in rows and SNPs in columns,
#'   missingness as -9
#'
#' @return an unnamed vector with counts
#'
#' @keywords internal

CalcSnpdBoth <- function(Pairs, GenoM) {
  sapply(seq_along(Pairs[,1]), function(i, G=GenoM) {
    sum(G[Pairs[i,1], ] != -9 & G[Pairs[i,2], ] != -9)
  })
}
