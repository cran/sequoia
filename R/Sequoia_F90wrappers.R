#============================================================================
#============================================================================
# preparations

#' Check parameter values and save as named vector.
#'
#'  Please do not increasing the number of SNPs or individuals beyond the
#'  numbers present in the datasets, as this may cause R to crash.
#'
#' @param GenoM matrix with genotype data, size nInd x nSnp
#' @param LifeHistData Dataframe with 3 columns:
#'  \itemize{
#'  \item{ID: }{max. 30 characters long,}
#'  \item{Sex: }{1 = females, 2 = males, other numbers = unkown,}
#'  \item{Birth Year: }{(or hatching year) Zero and negative numbers are
#'    interpreted as missing values.}}
#' @param nAgeClasses Number of age classes (= no. rows in AgePriors)
#' @param MaxSibIter Maximum number of iterations of sibship clustering
#'   (up to 42).
#' @param Err Estimated genotyping error rate.
#' @param MaxMismatch Maximum number of loci at which candidate parent and
#'  offspring are allowed to be opposite homozygotes, or be excluded.
#' @param Tfilter Threshold log-likelihood ratio between a proposed
#'   relationship versus unrelated, to select candidate relatives. Typically a
#'   negative value, related to the fact that unconditional likelihoods are
#'   calculated during the filtering steps. More negative values may decrease
#'   non-assignment, but will increase computational time.
#' @param Tassign Minimum log-likelihood ratio required for acceptance of
#'  proposed relationship, relative to next most likely relationship. Higher
#'  values result in more conservative assignments.
#' @param MaxSibshipSize  Maximum number of offspring for a single individual
#'  (a generous safety margin is advised).
#' @param DummyPrefix character vector of length 2 with prefixes for dummy
#'   dams (mothers) and sires (fathers); maximum 20 characters each.
#' @param Complexity Either "full" (default), "simp" (no explicit consideration
#'   of inbred relationships) or "mono" (monogamous breeding system)
#' @param FindMaybeRel  Identify pairs of non-assigned likely relatives after
#'   pedigree reconstruction. Can be time-consuming in large datasets.
#' @param CalcLLR  Calculate log-likelihood ratios for all assignments. Can be
#'  time-consuming in large datasets.
#' @return A 1-row dataframe with parameter values

SeqPrep <- function(GenoM = NULL,
                     LifeHistData = NULL,
                     nAgeClasses = 1,
                     MaxSibIter = 5,
                     Err = 0.0001,
                     MaxMismatch = 3,
                     Tfilter = -2.0,
                     Tassign = 0.5,
            				 MaxSibshipSize = 100,
            				 DummyPrefix = c("F", "M"),
            				 Complexity = "full",
            				 FindMaybeRel = TRUE,
            				 CalcLLR = TRUE)

{
  if (!is.matrix(GenoM)) stop("'GenoM' should be a numeric matrix")
  if (!all(GenoM %in% c(0,1,2,-9))) stop("'GenoM' not in <seq> format, please use GenoConvert")
  nIndG <- nrow(GenoM)
  IDs_geno <- rownames(GenoM)
  nSnp <- ncol(GenoM)

  if (Err<0 || Err>1 || !is.double(Err)) stop("invalid value for 'Err'")
  if (MaxMismatch<0 || !is.wholenumber(MaxMismatch))  stop("invalid value for 'MaxMismatch'")
  if (!is.double(Tfilter))  stop("invalid value for 'Tfilter'")
  if (!is.double(Tassign) || Tassign<0)  stop("invalid value for 'Tassign'")
  if (MaxSibshipSize<0 || !is.wholenumber(MaxSibshipSize)) {
    stop("invalid value for MaxSibshipSize")
  }
  if (!Complexity %in% c("full", "simp", "mono")) stop("invalid value for 'Complexity'")
  if (!FindMaybeRel %in% c(TRUE, FALSE)) stop("invalid value for 'FindMaybeRel'")
  if (!CalcLLR %in% c(TRUE, FALSE)) stop("invalid value for 'CalcLLR'")
  if (MaxSibIter<0 || !is.wholenumber(MaxSibIter))  stop("invalid value for 'MaxSibIter'")

  if (!is.null(LifeHistData)) {
    nIndLH <- nrow(LifeHistData)
    if (nIndLH >1 & any(LifeHistData$BY >= 0)) {
      nAgeClasses <- with(LifeHistData, diff(range(BY[BY >= 0 & ID %in% IDs_geno],
                                                 na.rm = TRUE))) + 1
    }
    if (nAgeClasses > 100) {
      nAgeClasses <- length(table(LifeHistData$BY[LifeHistData$BY >= 0]))
       if (nAgeClasses > 100) stop("Cannot handle >100 cohorts!")
    }
  }

  Specs <- data.frame(NumberIndivGenotyped = nIndG,
             NumberSnps = nSnp,
             GenotypingErrorRate = Err,
             MaxMismatch = MaxMismatch,
             Tfilter = Tfilter,
             Tassign = Tassign,
      			 nAgeClasses = nAgeClasses,
             MaxSibshipSize = MaxSibshipSize,
      			 MaxSibIter = MaxSibIter,
             DummyPrefixFemale = DummyPrefix[1],
             DummyPrefixMale = DummyPrefix[2],
             Complexity = Complexity,
      			 FindMaybeRel = FindMaybeRel,
      			 CalcLLR = CalcLLR,
      			 stringsAsFactors = FALSE)
  Specs
}


#============================================================================
#============================================================================

#' Check data for duplicates.
#'
#' Check the genotype and life history data for duplicate IDs (not permitted)
#' and duplicated genotypes (not advised), and count how many individuals in
#' the genotype data are not included in the life history data (permitted). The
#' order of IDs in the genotype and life history data is not required to be
#' identical.
#'
#' @param Specs The 1-row dataframe with parameter values
#' @param GenoM matrix with genotype data, size nInd x nSnp
#' @param LhIN  life history data
#' @param quiet suppress messages
#'
#' @return A list with one or more of the following elements:
#' \item{DupGenoID}{Dataframe, rownumbers of duplicated IDs in genotype data.
#'   Please do remove or relabel these to avoid downstream confusion.}
#' \item{DupGenotype}{Dataframe, duplicated genotypes (with or without
#'   identical IDs). The specified number of maximum mismatches is allowed,
#'   and this dataframe may include pairs of closely related individuals.}
#' \item{DupLifeHistID}{Dataframe, rownumbers of duplicated IDs in life
#'   history data}
#' \item{NoLH}{Vector, IDs (in genotype data) for which no life history data is
#'  provided}
#'
#' dataframe with the numbers of duplicated genotypes, duplicated IDs in the
#' genotype data, number of duplicated IDs in the life history data, and number
#'  of individuals with unknown sex.
#'
#' @useDynLib sequoia, .registration = TRUE
# @useDynLib sequoia duplicates

SeqDup <- function(Specs=NULL, GenoM = NULL, LhIN=NULL, quiet=FALSE)
{
  Ng <- FacToNum(Specs[,"NumberIndivGenotyped"])
  SpecsInt <- c(nSnp = FacToNum(Specs[,"NumberSnps"]),
                MaxMis = FacToNum(Specs[,"MaxMismatch"]),
                quiet = as.integer(quiet))
  gID <- rownames(GenoM)
  GenoV <- as.integer(GenoM)

  DUP <- .Fortran("duplicates",
                  Ng = as.integer(Ng),
                  SpecsInt = as.integer(SpecsInt),
                  GenoV = as.integer(GenoV),
                  nDupGenos = as.integer(0),
                  DupGenosFR = integer(2*Ng),  # N x 2 matrix
                  CntMism = integer(Ng))
#                  PACKAGE = "sequoia")

  Duplicates <- list()
  if (any(duplicated(gID))) {
    r1 <- which(duplicated(gID))
    r2 <- which(duplicated(gID, fromLast=TRUE))
    Duplicates$DupGenoID <- data.frame(row1 = r1,
                                           row2 = r2,
                                           ID = gID[r1])
  }
  if(DUP$nDupGenos>0) {
    tmp <- VtoM(DUP$DupGenosFR, DUP$nDupGenos)
    Duplicates$DupGenotype <- data.frame(row1 = tmp[, 1],
                                       row2 = tmp[, 2],
                                       ID1 = gID[tmp[, 1]],
                                       ID2 = gID[tmp[, 2]],
                                       Mismatch = DUP$CntMism[1:DUP$nDupGenos])
  }
  if (!is.null(LhIN)) {
    names(LhIN) <- c("ID", "Sex", "BY")
    LhIN$ID <- as.character(LhIN$ID)
    if (any(duplicated(LhIN[,1]))) {
      r1 <- which(duplicated(LhIN[,1]))
      r2 <- which(duplicated(LhIN[,1], fromLast=TRUE))
      Duplicates$DupLifeHistID <- data.frame(row1 = r1,
                                  row2 = r2,
                                  ID = LhIN[r1, "ID"],
                                  Sex1 = LhIN[r1, "Sex"],
                                  Sex2 = LhIN[r2, "Sex"],
                                  BY1 = LhIN[r1, "BY"],
                                  BY2 = LhIN[r2, "BY"])
    }
    NoLH <- setdiff(gID, LhIN$ID)
    if (length(NoLH)>0) Duplicates$NoLH <- NoLH
  }

  # print warnings
  if (!quiet) {
    if (any(duplicated(gID))) message("duplicate IDs found in genotype data, please remove to avoid confusion")
    if (DUP$nDupGenos>0 && DUP$nDupGenos > sum(duplicated(gID))) {
      message("likely duplicate genotypes found, consider removing")
    }
    if (any(duplicated(LhIN[,1]))) message("duplicate IDs found in lifehistory data, first entry will be used")
  }

  Duplicates
}


#=====================================================================
#=====================================================================

#' Fortran wrapper
#'
#' Call Fortran part and convert its output to a list with dataframes.
#'
#' @param ParSib either "par" to call parentage assignment, or "sib" to call
#'   the rest of the algorithm.
#' @param Specs a named vector with parameter values, as generated by
#'   \code{\link{SeqPrep}}.
#' @param GenoM matrix with genotype data, size nInd x nSnp
#' @param LhIN  life history data: ID - sex - birth year
#' @param AgePriors matrix with agepriors, size Specs["nAgeClasses"] by 8.
#' @param Parents  matrix with rownumbers of assigned parents, size nInd by 2
#' @param quiet suppress messages
#'
#' @return A list with
#' \item{PedigreePar or Pedigree}{the pedigree}
#' \item{MaybeParent or MaybeRel}{Non-assigned likely relatives}
#' \item{DummyIDs}{Info on dummies}
#' \item{TotLikParents or TotLikSib}{Total log-likelihood per iteration}
#'
#' For a detailed description of the output see \code{\link{sequoia}}
#'
#' @useDynLib sequoia, .registration = TRUE
# @useDynLib sequoia makeped deallocall

SeqParSib <- function(ParSib = "par",
                      Specs = NULL,
                      GenoM = NULL,
                      LhIN = NULL,
                      AgePriors = NULL,
                      Parents = NULL,
                      quiet = FALSE)
{
  on.exit(.Fortran("deallocall")) #, PACKAGE = "sequoia"))

  Ng <- FacToNum(Specs[,"NumberIndivGenotyped"])
  SMax <- FacToNum(Specs[,"MaxSibshipSize"])
  gID <- rownames(GenoM)
  GenoV <- as.integer(GenoM)
  LHL <- orderLH(LhIN[LhIN$Sex %in% c(1:3), ], gID)
  if (!is.null(Parents)) {
    PedPar <- IDToNum(Parents, gID)[, 2:3]
    PedPar <- c(as.matrix(PedPar))
    if (length(PedPar) != Ng*2) stop("'PedPar' wrong length")
  } else {
    PedPar <- rep(0, Ng*2)
  }
  Complex <- switch(Specs[,"Complexity"], full = 2, simp = 1, mono = 0)
  PrSb <- switch(ParSib, par = 1, sib = 2)
  nAmbMax <- 7*Ng  # max no. of non-assigned relative pairs to return

  SpecsInt <- c(ParSib = as.integer(PrSb),        # 1
                MaxSibIter = FacToNum(Specs[,"MaxSibIter"]), # 2
                nSnp = FacToNum(Specs[,"NumberSnps"]),       # 3
                MaxMis = FacToNum(Specs[,"MaxMismatch"]),    # 4
                SMax = SMax,                      # 5
                nAgeCl = FacToNum(Specs[,"nAgeClasses"]),    # 6
                Complx = as.integer(Complex),               # 7
                FindMaybe = as.integer(as.logical(Specs[,"FindMaybeRel"])),  # 8
                CalcLLR = as.integer(as.logical(Specs[,"CalcLLR"])),  # 9
                quiet = as.integer(quiet),        # 10
                nAmbMax = as.integer(nAmbMax))       # 11
  SpecsDbl <- c(Er = FacToNum(Specs[,"GenotypingErrorRate"]),
                TF = FacToNum(Specs[,"Tfilter"]),
                TA = FacToNum(Specs[,"Tassign"]))
  if (length(as.integer(AgePriors)) != 8*FacToNum(Specs[,"nAgeClasses"])) {
    stop("'AgePriors' matrix should have size nAgeClasses * 8")
  }

  TMP <- .Fortran("makeped",
                  Ng = as.integer(Ng),
                  SpecsInt = as.integer(SpecsInt),
                  SpecsDbl = as.double(SpecsDbl),
                  GenoV = as.integer(GenoV),
                  Sex = as.integer(LHL$Sex),
                  BY = as.integer(LHL$BY),
                  AP = as.double(AgePriors),
                  parentsRF = as.integer(PedPar),
                  LrRF = double(3*Ng),
                  OhRF = integer(2*Ng),
                  nAmb = as.integer(0),
                  AmbigID = integer(2*nAmbMax),
                  AmbigSex = integer(2*nAmbMax),
                  AmbigAgeDif = integer(nAmbMax),
                  AmbigRel = integer(2*nAmbMax),
                  AmbigLR = double(2*nAmbMax),
                  AmbigOH = integer(nAmbMax),
                  Nd = integer(2),
                  DumParRF = integer(2*Ng),
                  DumLrRF = double(3*Ng),
                  DumBYRF = integer(3*Ng),
                  DumNoff = integer(Ng),
                  DumOff = integer(SMax*Ng),
                  TotLik = double(42))
#                  PACKAGE = "sequoia")

  TMP$LrRF[abs(TMP$LrRF - 999) < 0.1] <- NA
  TMP$AmbigLR[abs(TMP$AmbigLR - 999) < 0.1] <- NA
  TMP$DumLrRF[abs(TMP$DumLrRF - 999) < 0.1] <- NA
  TMP$LrRF <- round(TMP$LrRF, 2)
  TMP$AmbigLR <- round(TMP$AmbigLR, 2)
  TMP$DumLrRF <- round(TMP$DumLrRF, 2)
  TMP$OhRF[TMP$OhRF < 0] <- NA
  TMP$AmbigOH[TMP$AmbigOH < 0] <- NA

  DumPfx <- Specs[,c("DummyPrefixFemale", "DummyPrefixMale")]
  dID <- cbind(paste0(DumPfx[1], formatC(1:Ng, width=4, flag=0)),
               paste0(DumPfx[2], formatC(1:Ng, width=4, flag=0)))

  #=========================
  # pedigree
  Pedigree <- data.frame(id = gID,
                         VtoM(TMP$parentsRF),
                         VtoM(TMP$LrRF, nc=3),
                         stringsAsFactors=FALSE)
  names(Pedigree) <- c("id", "dam", "sire", "LLRdam", "LLRsire", "LLRpair")
  Pedigree$LLRdam[is.na(Pedigree$dam)]  <- NA  # not sure why sometimes '0' - TODO CHECK
  Pedigree$LLRsire[is.na(Pedigree$sire)]  <- NA
  for (k in 1:2) Pedigree[, k+1] <- NumToID(Pedigree[, k+1], k, gID, dID)

  if (grepl("par", ParSib)) {
    Pedigree <- cbind(Pedigree,
                      VtoM(TMP$OhRF))
    names(Pedigree)[7:8] <- c("OHdam", "OHsire")
  }

  if (any(LhIN$Sex==4)) {  # hermaphrodites
    Pedigree <- herm_unclone_Ped(Pedigree, LH=LhIN, herm.suf=c("f", "m"))
  }

  if (!quiet) {
    if (grepl("par", ParSib)) {
      message("assigned ", sum(!is.na(Pedigree$dam)), " dams and ",
           sum(!is.na(Pedigree$sire)), " sires to ", nrow(Pedigree), " individuals")
   } else {
     message("assigned ", sum(!is.na(Pedigree$dam)), " dams and ",
           sum(!is.na(Pedigree$sire)), " sires to ", Ng, " + ", sum(TMP$Nd),
           " individuals (real + dummy)")
   }
  }


  #=========================
  # non-assigned probable relatives
  if (TMP$nAmb > 0) {
    RelName <- c("PO", "FS", "HS", "GP", "FA", "HA", "U ", "XX", "2nd")
    Na <- TMP$nAmb
    TMP$AmbigID <- NumToID(TMP$AmbigID, 0, gID, dID)
    AmbigRel <- factor(TMP$AmbigRel, levels=1:9, labels=RelName)

    MaybeRel <- data.frame(VtoM(TMP$AmbigID, Na),
                           VtoM(TMP$AmbigSex, Na),
                           TMP$AmbigAgeDif[1:Na],
                           VtoM(AmbigRel, Na),
                           VtoM(TMP$AmbigLR, Na),
                           stringsAsFactors=FALSE)
    names(MaybeRel) <- c("ID1", "ID2", "Sex1", "Sex2", "AgeDif",
                         "Relx", "TopRel", "LLR_Rx_U", "LLR")
    MaybeRel <- MaybeRel[,-which(names(MaybeRel) %in% c("Relx", "LLR_Rx_U"))]  # drop; confusing.
    for (i in 1:Na) {
      if (MaybeRel$AgeDif[i] < 0) {
        tmp <- MaybeRel[i,]
        tmp$AgeDif <- abs(tmp$AgeDif)
        MaybeRel[i,] <- tmp[,c(2,1,4,3,5:7)]
      }
    }
    if (grepl("par", ParSib)) {
      MaybeRel$OH <-  TMP$AmbigOH[1:Na]
#      MaybeRel <- with(MaybeRel, MaybeRel[which(TopRel=="PO"), ])
    } else {
      MaybeRel <- with(MaybeRel, MaybeRel[TopRel %in% c("PO", "FS", "HS","GG") |
                                          LLR > FacToNum(Specs[,"Tassign"]), ])
    }
    MaybeRel <- MaybeRel[order(ordered(MaybeRel$TopRel, levels=RelName),
                                     -MaybeRel$LLR),]
    MaybeRel$AgeDif[MaybeRel$AgeDif==999] <- NA

    if (any(LhIN$Sex==4)) {  # hermaphrodites
      MaybeRel <- herm_unclone_MaybeRel(MaybeRel, Pedigree, LH=LhIN, herm.suf=c("f", "m"))
    }

   if (!quiet && nrow(MaybeRel)>0) {
     if (grepl("par", ParSib)) {
     message("there are  ", sum(MaybeRel$TopRel=="PO"),
             "  likely parent-offspring pairs and ", nrow(MaybeRel)-sum(MaybeRel$TopRel=="PO"),
             " other pairs of likely relatives  which are not assigned, ",
             "perhaps due to unknown birth year(s), please see 'MaybeParent'")
     } else {
       message("there are  ", nrow(MaybeRel), "  non-assigned pairs of possible relatives, ",
               "including  ", sum(MaybeRel$TopRel=="PO"), "  likely parent-offspring pairs; ",
               "please see 'MaybeRel'")
     }
  }
  } else  MaybeRel <- NULL


  #=========================
  # dummies
  if (grepl("sib", ParSib) && any(TMP$Nd>0)) {
    nd <- TMP$Nd
    NgOdd <- Ng%%2==1
    for (k in 1:2) if (nd[k]==0)  nd[k] <- 1  # easier; remove below
    mso <- max(TMP$DumNoff)
    TMP$DumOff <- NumToID(TMP$DumOff, 0, gID, dID)  # contains no dummies currently

    DummyIDs <- data.frame(id=c(dID[1:nd[1], 1], dID[1:nd[2], 2]),
                           VtoM(TMP$DumParRF, sum(nd), 2, NgOdd),
                           VtoM(TMP$DumLrRF, sum(nd), 3, NgOdd),
                           sex=rep(1:2, nd),
                           VtoM(TMP$DumBYRF, sum(nd),3, NgOdd),
                           TMP$DumNoff[1:sum(nd)],
                           VtoM(TMP$DumOff, sum(nd), SMax, NgOdd)[, 1:mso],
                           stringsAsFactors=FALSE)
    names(DummyIDs) <-  c("id", "dam", "sire", "LLRdam", "LLRsire", "LLRpair",
                          "sex", "BY.est", "BY.min", "BY.max", "NumOff",
                          paste0("O", 1:mso))
    if (TMP$Nd[1]==0) DummyIDs <- DummyIDs[-1, ]
    if (TMP$Nd[2]==0) DummyIDs <- DummyIDs[-nrow(DummyIDs), ]
    for (k in 1:2) DummyIDs[, k+1] <- NumToID(DummyIDs[, k+1], k, gID, dID)

    Pedigree <- rbind(Pedigree, DummyIDs[, 1:6])
  } else  DummyIDs <- NULL


  #=========================
  # output
  if (grepl("par", ParSib)) {
    return(list(PedigreePar = Pedigree,
        MaybeParent = MaybeRel,
        TotLikParents = TMP$TotLik[1:sum(TMP$TotLik!=0)]))
  } else if (grepl("sib", ParSib)) {
    return(list(Pedigree = Pedigree,
                DummyIDs = DummyIDs,
                MaybeRel = MaybeRel,
                TotLikSib = TMP$TotLik[1:sum(TMP$TotLik!=0)]))
  }
}


#=====================================================================
#=====================================================================


#===============================
#' Order lifehistory data
#'
#' Order lifehistory data to match order of IDs in genotype data,
#' filling in gaps with missing values
#'
#' @param LH dataframe with lifehistory information:
#' \itemize{
#'  \item{ID: }{max. 30 characters long,}
#'  \item{Sex: }{1 = females, 2 = males, other numbers = unkown,}
#'  \item{Birth Year: }{(or hatching year) Use negative numbers to denote
#'  missing values.}}
#' @param gID character vector with IDs in genotype data, in order of occurence
#'
#' @return
#' \item{BY}{Numeric vector with birth years, of same length as gID}
#' \item{Sex}{Numeric vector with genders, of same length as gID}

orderLH <- function(LH=NULL, gID=NULL) {
  if (!is.null(LH)) {
    names(LH) <- c("lhID", "Sex", "BY")
    LH$lhID <- as.character(LH$lhID)
    LH <- LH[!duplicated(LH$lhID), ]
    rownames(LH) <- LH$lhID
    LH <- LH[gID, ]
    Sex <- LH$Sex
    Sex[is.na(Sex)] <- 3
    Sex[Sex<1 | Sex>2] <- 3
    BY <- LH$BY
    BY[is.na(BY)] <- -999
    BY[BY<0] <- -999
  } else {
    Sex <- rep(3, length(gID))
    BY <- rep(-999, length(gID))
  }
  list(Sex=Sex, BY=BY)
}

