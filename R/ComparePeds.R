#============================================================================
#============================================================================
#' @title Compare two Pedigrees
#'
#' @description Compare an inferred pedigree (Ped2) to a previous or simulated
#'   pedigree (Ped1), including comparison of sibship clusters and sibship
#'   grandparents.
#'
#' @details The comparison is divided into different classes of `assignable'
#'   parents. This includes cases where the focal individual and parent
#'   according to Ped1 are both Genotyped (G-G), as well as cases where the
#'   non-genotyped parent according to Ped1 can be lined up with a sibship Dummy
#'   parent in Ped2 (G-D), or where the non-genotyped focal individual in Ped1
#'   can be matched to a dummy individual in Ped2 (D-G and D-D). If SNPd is NULL
#'   (the default), and DumPrefix is set to NULL, the intersect between the IDs
#'   in Pedigrees 1 and 2 is taken as the vector of genotyped individuals.
#'
#' @param  Ped1 original pedigree, dataframe with columns id-dam-sire; only the
#'   first 3 columns will be used.
#' @param  Ped2 infered pedigree, e.g. \code{SeqOUT$Pedigree}, with columns
#'   id-dam-sire.
#' @param  na1  the value for missing parents in Ped1 (assumed NA in Ped2).
#' @param  DumPrefix  character vector of length 2 with the dummy prefices in
#'   Pedigree 2; all IDs not starting with the Dummy prefix are taken as
#'   genotyped.
#' @param SNPd character vector with IDs of genotyped individuals.
#'
#' @return A list with \item{Counts}{A 7 x 5 x 2 named numeric array with the
#'   number of matches and mismatches}
#' \item{MergedPed}{A side-by-side comparison of the two pedigrees}
#' \item{ConsensusPed}{A consensus pedigree, with Pedigree 2 taking priority
#'   over Pedigree 1}
#' \item{DummyMatch}{Dataframe with all dummy IDs in Pedigree 2 (id), and the
#'  best-matching individual in Pedigree 1 (id.r)}
#' \item{Mismatch}{A subset of MergedPed with mismatches between Ped1 and Ped2,
#'  as defined below. The two additional columns are Cat (category, 'GG', 'GD',
#'  'DG' or 'DD', as described below) and Parent ('dam' or 'sire' indicating
#'  which is mismatching)}
#' \item{Ped1only}{as Mismatches, with parents in Ped1 that were not assigned
#'   in Ped2}
#' \item{Ped2only}{as Mismatches, with parents in Ped2 that were missing in
#'   Ped1}
#'
#' The first dimension of \code{Counts} denotes the following categories:
#' \item{GG}{Genotyped individual, assigned a genotyped parent in either
#'   pedigree}
#' \item{GD}{Genotyped individual, assigned a dummy parent, or at least 1
#'   genotyped sibling or a genotyped grandparent in Pedigree 1)}
#' \item{GT}{Genotyped individual, total}
#' \item{DG}{Dummy individual, assigned a genotyped parent (i.e., grandparent
#'    of the sibship in Pedigree 2)}
#' \item{DD}{Dummy individual, assigned a dummy parent (i.e., avuncular
#'   relationship between sibships in Pedigree 2)}
#' \item{DT}{Dummy total}
#' \item{TT}{Total total, includes all genotyped individuals, plus
#'   non-genotyped individuals in Pedigree 1, plus non-replaced dummy
#'   individuals (see below) in Pedigree 2}
#'
#' The dummy individual count includes all non-genotyped individuals in
#' Pedigree 1 who have, according to either pedigree, at least 2 genotyped
#' offspring, or at least one genotyped offspring and a genotyped parent.
#'
#' The second dimension of \code{Counts} gives the outcomes:
#' \item{Total}{The total number of individuals with a parent assigned in
#'    either or both pedigrees}
#' \item{Match}{The same parent is assigned in both pedigrees (non-missing).
#'      For dummy parents, it is considered a match if the inferred sibship
#'      which contains the most offspring of a non-genotyped parent, consists
#'      for more than half of this individual's offspring.}
#' \item{Mismatch}{Different parents assigned in the two pedigrees. When
#'    a sibship according to Pedigree 1 is split over two sibships in Pedigree
#'    2, the smaller fraction is included in the count here.}
#' \item{P1only}{Parent in Pedigree 1 but not 2; includes non-assignable
#'     parents (e.g. not genotyped and no genotyped offspring).}
#' \item{P2only}{Parent in Pedigree 2 but not 1.}
#'
#' The third dimension \code{Counts} separates between maternal and paternal
#' assignments, where e.g. paternal 'DR' is the assignment of fathers to both
#' maternal and paternal sibships.
#'
#' 'MergedPed' provides the following columns:
#' \item{id}{All ids in both Pedigree 1 and 2}
#' \item{dam.1, sire.1}{parents in Pedigree 1}
#' \item{dam.2, sire.2}{parents in Pedigree 2}
#' \item{id.r, dam.r, sire.r}{when in Pedigree 2 a dummy parent is assigned,
#'   this column gives the best-matching non-genotyped individual according to
#'   Pedigree 1, or "nomatch". If a sibship in Pedigree 1 is divided over 2
#'   sibships in Pedigree 2, the smaller one will be denoted as "nomatch"}
#'
#' In 'ConsensusPed', the priority used is parent.r (if not "nomatch) >
#'   parent.2 > parent.1. The columns 'dam.cat' and 'sire.cat' give a 2-letter
#'   code denoting whether the focal individual (first letter) and its assigned
#'   parent (2nd letter) are
#'   \item{G}{Genotyped}
#'   \item{D}{Dummy individual (in Pedigree 2)}
#'   \item{R}{Dummy individual in pedigree 2 replaced by best matching
#'    non-genotyped individual in pedigree 1}
#'   \item{U}{Ungenotyped (in Pedigree 1, with no dummy match)}
#'   \item{X}{No parent in either pedigree}
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{DyadCompare}, \link{sequoia}, \link{EstConf}}.
#'
#' @examples
#' \dontrun{
#' data(Ped_HSg5, SimGeno_example, LH_HSg5, package="sequoia")
#' SeqOUT <- sequoia(GenoM = SimGeno_example, LifeHistData = LH_HSg5)
#' compare <- PedCompare(Ped1=Ped_HSg5, Ped2=SeqOUT$Pedigree)
#' compare$Counts   # 2 mismatches, due to simulated genotyping errors
#' head(compare$MergedPed)
#'
#' PedM <- compare$MergedPed
#' # find mismatching mothers:
#' with(PedM, PedM[which(dam.1!=dam.2 & dam.1!=dam.r),])
#'
#' # find mothers in Ped1 which are genotyped but not assigned in Ped2:
#' with(PedM, PedM[which(is.na(dam.2) & !is.na(dam.1) &
#'                        !is.na(id) & dam.1 %in% id),])
#' }
#' @export

PedCompare <- function(Ped1 = NULL,
                       Ped2 = NULL,
                       na1 = c(NA, "0"),
                       DumPrefix = c("F0", "M0"),
                       SNPd = NULL)
{
  if(is.null(Ped1)) stop("No 'Ped1' provided")
  if(is.null(Ped2)) stop("No 'Ped2' provided'")
  names(Ped1)[1:3] <- c("id", "dam.1", "sire.1")
  names(Ped2)[1:3] <- c("id", "dam.2", "sire.2")
  for (i in 1:3) {
    Ped1[, i] <- as.character(Ped1[, i])
    Ped1[Ped1[, i] %in% na1, i] <- NA
  }
  for (i in 1:3) Ped2[, i] <- as.character(Ped2[, i])
  if (!any(Ped2$id %in% Ped1$id))  stop("no common IDs in Ped1 and Ped2")
  Ped1 <- Ped1[!is.na(Ped1$id), 1:3]
  Ped2 <- Ped2[!is.na(Ped2$id), 1:3]
  Ped1 <- AddParPed(Ped1, ZeroToNA=TRUE)
  Ped2 <- AddParPed(Ped2, ZeroToNA=TRUE)
  if (is.null(DumPrefix) & is.null(SNPd)) {
    SNPd <- intersect(Ped2$id, Ped1$id)
  } else if (is.null(SNPd)) {
    DPnc <- nchar(DumPrefix)
    SNPd <- Ped2$id[substr(Ped2$id,1,DPnc[1])!=DumPrefix[1] &
                      substr(Ped2$id,1,DPnc[2])!=DumPrefix[2]]
  }
  if (sum(SNPd %in% Ped2$id)==0)  stop("none of 'SNPd' in Ped2")
  SNPd <- stats::na.exclude(SNPd)
  PedX <- merge(Ped1, Ped2[Ped2$id %in% SNPd, ], all.y=TRUE)
  DumPed <- Ped2[!Ped2$id %in% SNPd, ]
  Dummies <- list(DumPed$id[DumPed$id %in% Ped2$dam],
                  DumPed$id[DumPed$id %in% Ped2$sire])
  Par <- matrix(c("dam.1", "dam.2", "dam.r",
                  "sire.1", "sire.2", "sire.r"), 2,byrow=TRUE)

  #===  match dummies to non-genotyped parents  ===
  Sibs1 <- list()
  SibScore <- list()
  DumReal <- list()
  for (p in 1:2) {
    Sibs1[[p]] <- split(PedX$id, PedX[, Par[p,1]])
    Sibs1[[p]] <- Sibs1[[p]][!names(Sibs1[[p]]) %in% SNPd]
  }
  NoDummies <- with(PedX, all(dam.2 %in% SNPd | is.na(dam.2)) &
                      all(sire.2 %in% SNPd | is.na(sire.2)))
  if (!NoDummies) {
    for (p in 1:2) {
      Sibs2 <- split(PedX$id, PedX[, Par[p,2]])
      Sibs2 <- Sibs2[!names(Sibs2) %in% SNPd]

      if (length(Sibs1[[p]])>0) {
        SibScore[[p]] <- tdf(sapply(Sibs1[[p]], SibMatch, Sibs2, SNPd))  # slow!
        if (length(stats::na.exclude(SibScore[[p]][, "Best"])) >
            length(unique(stats::na.exclude(SibScore[[p]][, "Best"])))) {
          SibScore[[p]] <- SibScore[[p]][order(SibScore[[p]][, "OK"],
                                               decreasing=TRUE), ]
          dups <- duplicated(SibScore[[p]][, "Best"], incomparables=NA)
          BestS.d <- unique(SibScore[[p]]$Best[dups])
          if (sum(dups) > 1) {
            SibScore[[p]][dups, "Er"] <- rowSums(SibScore[[p]][dups, c("OK", "Er")])
          } else {
            SibScore[[p]][dups, "Er"] <- sum(SibScore[[p]][dups, c("OK", "Er")])
          }
          SibScore[[p]][dups, c("Best", "OK")] <- NA
          for (s in 1:length(BestS.d)) {
            tmp <- SibScore[[p]][which(SibScore[[p]]$Best==BestS.d[s]), ]
            if (tmp$OK==1 & tmp$Er>=1) {
              SibScore[[p]][which(SibScore[[p]]$Best==BestS.d[s]),
                            c("Best", "OK")] <- NA
            }
          }
        }
        tmp <- SibScore[[p]][!is.na(SibScore[[p]][, "Best"]), ]
        DumReal[[p]] <- as.data.frame(cbind(real = rownames(tmp),
                              dummy = tmp$Best), stringsAsFactors=FALSE)
        DumReal[[p]] <- merge(DumReal[[p]], DumPed[DumPed$id %in% Dummies[[p]],],
                              by.x="dummy", by.y="id", all=TRUE)[,c("real", "dummy")]
        DumReal[[p]][,"real"][is.na(DumReal[[p]][,"real"])] <- "nomatch"
      } else {
        SibScore[[p]] <- NA
        DumReal[[p]] <- data.frame(real=NA, dummy=NA)
      }
    }
    for (p in 1:2) {
      tmp <- stats::setNames(DumReal[[p]], c(Par[p,3], Par[p,2]))
      PedX <- merge(PedX, tmp, all.x=TRUE)
      DumPed <- merge(DumPed, DumReal[[p]], by.x="id", by.y="dummy", all.x=TRUE,
                      suffixes = c(".x",".y"))
    }
    DumPed$id.r <- apply(DumPed[,c("real.x", "real.y")], 1, function(x)
      ifelse(all(is.na(x)), NA, min(x, na.rm=T)))
    DumPed <- DumPed[, c("id", "dam.2", "sire.2", "id.r")]
    DumPed <- merge(DumPed, Ped1, all.x=TRUE, by.x="id.r", by.y="id",
                    suffixes=c(".2", ".1"))
    for (p in 1:2) {
      DumPed <- merge(DumPed, stats::setNames(DumReal[[p]], c(Par[p,3], Par[p,2])),
                      all.x=TRUE)
    }
  } else {
    PedX$dam.r <- NA
    PedX$sire.r <- NA
  }

  #===  Combined pedigree  ===
  PedX$id.r <- PedX$id
  PedY <- merge(PedX, DumPed, all=TRUE, sort=FALSE)  # NA's for id.r = "nomatch"
  PedY <- merge(PedY, stats::setNames(Ped1, c("id.r", "dam.1", "sire.1" )),
                all=TRUE, sort=FALSE)

  for (p in 1:2) {
    PedTmp <- PedY[which(!PedY[,Par[p,2]] %in% SNPd &
                           !is.na(PedY[,Par[p,2]]) &
                           PedY[,Par[p,3]]!="nomatch"),]
    tbl <- table(PedTmp[,Par[p,2]])
    npar1 <- plyr::daply(PedTmp, Par[p,2],
                   function(x) length(unique(stats::na.exclude(x[,Par[p,1]]))))
    unmatch <- names(which(c(tbl)==npar1 & npar1>1))
    if (length(unmatch)>0) {
      PedY[which(PedY[,Par[p,2]] %in% unmatch), Par[p,3]] <- "nomatch"
    }
  }
  PedY <- PedY[, c("id", "dam.1", "sire.1", "dam.2", "sire.2",
                 "id.r", "dam.r", "sire.r")]

  Founders <- list()
  NGpar <- list()
  for (p in 1:2) {
    Founders[[p]] <- unique(unlist(PedY[is.na(PedY[,Par[p,1]]) & is.na(PedY[,Par[p,2]]),
                                 c("id", "id.r")]))
    Sibs1[[p]] <- Sibs1[[p]][sapply(Sibs1[[p]],length)>1  |
                               !names(Sibs1[[p]]) %in% Founders[[p]]]
    NGpar[[p]] <- names(Sibs1[[p]])  # potential dummies
  }
  GD <- list(G = SNPd, D = c(unlist(NGpar), unlist(Dummies)) )

  Score <- array(0, dim=c(7, 5, 2),
                 dimnames=list(c("GG", "GD", "GT", "DG", "DD", "DT", "TT"),
            c("Total", "Match", "Mismatch", "P1only", "P2only"),
            c("dam", "sire")))
  NoMatch <- list("Mismatch" = cbind(PedY[1:2, ], Cat = NA, Parent = NA)[0,],
                  "P1only" = cbind(PedY[1:2, ], Cat = NA, Parent = NA)[0,],
                  "P2only" = cbind(PedY[1:2, ], Cat = NA, Parent = NA)[0,])
  ID <- c(G="id", D="id.r")

  for (p in 1:2) {
    for (i in c("G", "D")) {  # focal
      for (j in c("G", "D")) {  # parent
        if ((i=="D" | j=="D") & NoDummies) break
        ij <- paste0(i,j)
        k <- ifelse(j=="G", 2, 3)
        PedTmp <- PedY[which(PedY[,ID[i]] %in% GD[[i]] &
                          (PedY[,Par[p,1]] %in% GD[[j]] | PedY[,Par[p,2]] %in% GD[[j]])), ]
        Score[ij, "Total", p] <- nrow(PedTmp)
        Score[ij, "Match", p] <- sum(PedTmp[,Par[p,1]] == PedTmp[,Par[p,k]], na.rm=T)

        NoM <- list()
        NoM[["Mismatch"]] <- PedTmp[which(PedTmp[,Par[p,1]] != PedTmp[,Par[p,k]]), ]
        NoM[["P1only"]] <- PedTmp[which(!is.na(PedTmp[,Par[p,1]]) & PedTmp[,Par[p,1]]
                                         %in% GD[[j]] & is.na(PedTmp[,Par[p,2]])), ]
        NoM[["P2only"]] <- PedTmp[which(!is.na(PedTmp[,Par[p,2]]) & PedTmp[,Par[p,2]]
                                         %in% GD[[j]] & is.na(PedTmp[,Par[p,1]])), ]

        for (x in c("Mismatch", "P1only", "P2only")) {
          Score[ij, x, p] <- nrow(NoM[[x]])
          if (nrow(NoM[[x]])>0) {
            NoMatch[[x]] <- rbind(NoMatch[[x]],
                                cbind(NoM[[x]], Cat = ij, Parent=c("dam", "sire")[p]))
          }
        }

      }
    }
  }

  for (p in 1:2) {  # Totals
    for (i in c("G", "D", "T")) {
      if (i=="D" & NoDummies) next
      ij <- paste0(i,"T")
      if (i != "T") {
        PedTmp <- PedY[which(PedY[,ID[i]] %in% GD[[i]] &
                            (!is.na(PedY[,Par[p,1]]) | !is.na(PedY[,Par[p,2]]))), ]
      } else {
        PedTmp <- PedY[!is.na(PedY[,Par[p,1]]) | !is.na(PedY[,Par[p,2]]), ]
      }
      Score[ij, "Total", p] <- nrow(PedTmp)
      Score[ij, "Match", p] <- sum(PedTmp[,Par[p,1]] == PedTmp[,Par[p,2]] |
                                        PedTmp[,Par[p,1]] == PedTmp[,Par[p,3]], na.rm=T)

      NoM <- list()
      NoM[["Mismatch"]] <- PedTmp[which(PedTmp[,Par[p,1]] != PedTmp[,Par[p,2]] &
                  (PedTmp[,Par[p,1]] != PedTmp[,Par[p,3]] | is.na(PedTmp[,Par[p,3]]))), ]
      NoM[["P1only"]] <- PedTmp[which(!is.na(PedTmp[,Par[p,1]]) & is.na(PedTmp[,Par[p,2]])), ]
      NoM[["P2only"]] <- PedTmp[which(is.na(PedTmp[,Par[p,1]]) & !is.na(PedTmp[,Par[p,2]])), ]

      for (x in c("Mismatch", "P1only", "P2only")) {
        Score[ij, x, p] <- nrow(NoM[[x]])
        if (nrow(NoM[[x]])>0) {
          NoMatch[[x]] <- rbind(NoMatch[[x]],
                              cbind(NoM[[x]], Cat = ij, Parent=c("dam", "sire")[p]))
        }
      }
    }
  }
  for (x in c("Mismatch", "P1only", "P2only")) {
    NoMatch[[x]] <- NoMatch[[x]][!duplicated(NoMatch[[x]][, names(NoMatch[[x]])!="Cat"]), ]
  }

  for (p in 1:2) {
    ParCat <- ifelse(is.na(PedY[,Par[p,2]]) & !is.na(PedY[,Par[p,1]]), "U",
                ifelse(PedY$id %in% SNPd,
                 ifelse(PedY[,Par[p,2]]  %in% SNPd, "GG",
                  ifelse(PedY[,Par[p,2]]  %in% Dummies[[p]],
                    ifelse(PedY[,Par[p,3]]=="nomatch", "GD", "GR"),
                    NA)),
                ifelse(PedY$id %in% unlist(Dummies),
                  ifelse(PedY$id.r == "nomatch",
                    ifelse(PedY[,Par[p,2]]  %in% SNPd, "DG",
                      ifelse(PedY[,Par[p,3]]=="nomatch", "DD", "DR")),
                    ifelse(PedY[,Par[p,2]]  %in% SNPd, "RG",
                      ifelse(PedY[,Par[p,3]]=="nomatch", "RD", "RR"))),
                  NA)))
    if (p==1)  PedY$dam.cat <- ParCat
    if (p==2)  PedY$sire.cat <- ParCat
  }
  PedY$id.r[which(PedY$id==PedY$id.r)] <- NA

  PedC <- with(PedY, data.frame(
    id = ifelse(!is.na(id.r) & id.r!="nomatch", id.r, id),
    dam = ifelse(!is.na(dam.r) & dam.r!="nomatch", dam.r,
                 ifelse(!is.na(dam.2), dam.2, dam.1)),
    sire = ifelse(!is.na(sire.r) & sire.r!="nomatch", sire.r,
                  ifelse(!is.na(sire.2), sire.2, sire.1)),
    dam.cat = as.character(dam.cat),
    sire.cat = as.character(sire.cat),
    stringsAsFactors = FALSE))
  PedC$dam.cat[is.na(PedC$dam.cat)] <- "X"
  PedC$sire.cat[is.na(PedC$sire.cat)] <- "X"


  #==============================
  # out
  OUT <- list(Counts = Score, MergedPed = PedY[, 1:8], ConsensusPed = PedC)
  if (!NoDummies) {
    OUT[["DummyMatch"]] <-  DumPed[order(DumPed$id), c("id", "id.r")]
  }
  c(OUT, NoMatch)
}


#============================================================================
#============================================================================
# Utils functions for comparisons


#' Compare two vectors
#'
#' Compare a vector with inferred sibs to a vector of `true' sibs
#'
#' @param  Infrd  vector of inferred sibs
#' @param  Simld  vector of true sibs
#' @param  SNPd character vector with IDs of genotyped individuals
#'
#' @return a named numeric vector of length 4, with the total length of Simld,
#'   the length of the intersect of the two vectors, the number occurring in
#'   Infrd but not Simld ('err'), and the number occuring in Simld but not
#'   Infrd ('missed').
#'
#' @keywords internal

Vcomp <- function(Infrd, Simld, SNPd)
{
  out <- c("total" = length(Simld),
           "both" = length(intersect(Infrd, Simld)),
           "err" = length(setdiff(Infrd[Infrd %in% SNPd], Simld)),
           "missed" = length(setdiff(Simld, Infrd)))
  out
}


#======================================================================
#' Find the closest matching inferred sibship to a true sibship
#'
#' @param SimX  a vector with the IDs in the true sibship
#' @param Infrd  a list of vectors with the IDs in the inferred sibships
#' @param SNPd character vector with IDs of genotyped individuals
#'
#' @return a named numeric vector with the number of matches ('NumMatch'),
#'   the position of the best match ('Best'), the inferred sibship size of
#'   this best match ('Tot'), the number of matching IDs ('OK'), and the
#'   number of mismatches ('err').
#'
#' @keywords internal

SibMatch <- function(SimX, Infrd, SNPd)
{
  VC <- sapply(Infrd, Vcomp, SimX, SNPd)
  mtch <- which(VC["both",]>0)
  if (length(mtch)>1) {
    if (VC["total",1]==2) {
      mtch <- NULL
    } else {
      mtch <- which.max(VC["both",])[1]  # which inferred sibship has most members of true sibship
    }
  }
  if (length(mtch)==1) {
    if (VC["err", mtch] > VC["both", mtch]) {  # was: >=
      mtch <- NULL
    }
  }
  if (length(mtch)==1) {
    OUT <- data.frame(NumMatch = sum(VC[2,]>0),
             Best = colnames(VC)[mtch],  # as.numeric(substr(colnames(VC)[mtch],3,5)),
             Tot = VC["total", mtch],
             OK = VC["both", mtch],
             Er = VC["err", mtch], stringsAsFactors=FALSE)
  } else {
    OUT <- data.frame(NumMatch=0, Best=NA,  Tot= length(SimX), OK = NA, Er = NA)
  }
  OUT
}


#======================================================================
# transpose matrix to data.frame  (sapply returns wrong format)

# adapted from Hmisc::all.is.numeric
all.is.numeric <- function (x, what = c("test", "vector"))
{
  what <- match.arg(what)
  isnum <- suppressWarnings(!any(is.na(as.numeric(stats::na.exclude(x)))))
  if (what == "test")
    isnum
  else if (isnum)
    as.numeric(x)
  else x
}


tdf <- function(M)
{
  DF <- matrix(NA, ncol(M), nrow(M))
  for (r in 1:nrow(M)) {
    DF[, r] <- unlist(M[r, ])
  }
  DF <- as.data.frame(DF, stringsAsFactors = FALSE, row.names=colnames(M))
  for (r in 1:ncol(DF)) {
    if (all.is.numeric(DF[, r]))  DF[, r] <- as.numeric(DF[, r])
  }
  names(DF) <- rownames(M)
  DF
}


#============================================================================
#============================================================================
#' Find siblings

#' @param x  an ID
#' @param Ped  a pedigree with columns id - dam - sire
#'
#' @return The individuals which are full or half siblings to x, as a
#'   three-column matrix with column names id1 (x), id2 (the siblings), and
#'   RC (the relatedness category, 'FS' or 'HS').
#'
#' @keywords internal

rc <- function(x, Ped) {
  names(Ped) <- c("id", "dam", "sire")
  RelCat <- with(Ped,
                 ifelse(id == id[x], "S",
                        ifelse(eqv(dam[x],dam,FALSE) & eqv(sire[x], sire,FALSE), "FS",
                               ifelse(eqv(dam[x],dam,FALSE) |  eqv(sire[x], sire,FALSE), "HS",
                                      NA))))
  out <- cbind(id1 = Ped$id[x],
               id2 = Ped$id[!is.na(RelCat)],
               RC = stats::na.exclude(RelCat))
  out <- out[out[,"RC"] != "S", ]
  out
}


#============================================================================
#' @title Compare dyads
#'
#' @description Count the number of half and full sibling pairs correctly and
#'   incorrectly assigned
#'
#' @param  Ped1 Original pedigree, dataframe with 3 columns: id-dam-sire
#' @param  Ped2 Second (inferred) pedigree
#' @param  na1  the value for missing parents in Ped1.
#'
#' @return A 3x3 table with the number of pairs assigned as full siblings (FS),
#'   half siblings (HS) or unrelated (U, including otherwise related) in the two
#'   pedigrees, with the classification in Ped1 on rows and that in Ped2 in
#'   columns
#'
#' @seealso \code{\link{PedCompare}}
#'
#' @examples
#' \dontrun{
#' data(Ped_HSg5, SimGeno_example, LH_HSg5, package="sequoia")
#' SeqOUT <- sequoia(GenoM = SimGeno_example, LifeHistData = LH_HSg5,
#'                   MaxSibIter = 0)
#' DyadCompare(Ped1=Ped_HSg5, Ped2=SeqOUT$Pedigree)
#' }
#'
#' @export

DyadCompare <- function(Ped1 = NULL,
                       Ped2 = NULL,
                       na1 = c(NA, "0"))
{
  if(is.null(Ped1)) stop("No 'Ped1' provided")
  if(is.null(Ped2)) stop("No 'Ped2' provided'")
  names(Ped1)[1:3] <- c("id", "dam.1", "sire.1")
  names(Ped2)[1:3] <- c("id", "dam.2", "sire.2")
  for (i in 1:3) {
    Ped1[, i] <- as.character(Ped1[, i])
    Ped1[Ped1[, i] %in% na1, i] <- NA
  }
  for (i in 1:3) Ped2[, i] <- as.character(Ped2[, i])
  if (!any(Ped2$id %in% Ped1$id))  stop("no common IDs in Ped1 and Ped2")
  Ped1 <- Ped1[!is.na(Ped1$id), 1:3]
  Ped2 <- Ped2[!is.na(Ped2$id), 1:3]
  Ped1 <- AddParPed(Ped1)
  Ped2 <- AddParPed(Ped2)

  # note: each pair is counted double
  RCT <- matrix(NA, 0, 3)
  for (x in 1:nrow(Ped1)) {
    RCT <- rbind(RCT, rc(x, Ped1))
  }

  RCI <- matrix(NA, 0, 3)
  for (x in 1:nrow(Ped2)) {
    RCI <- rbind(RCI, rc(x, Ped2))
  }

  RCTI <- merge(as.data.frame(RCT, stringsAsFactors=FALSE),
                as.data.frame(RCI, stringsAsFactors=FALSE),
                by=c("id1", "id2"), all=TRUE, suffixes = c(".1", ".2"))
  RCTI <- RCTI[RCTI$id1 %in% Ped1$id & RCTI$id2 %in% Ped1$id &
                 RCTI$id1 %in% Ped2$id & RCTI$id2 %in% Ped2$id, ]
  RCTI$RC.1[is.na(RCTI$RC.1)] <- "U"
  RCTI$RC.2[is.na(RCTI$RC.2)] <- "U"
  RCTI$RC.1 <- factor(RCTI$RC.1, levels=c("FS", "HS", "U"))
  RCTI$RC.2 <- factor(RCTI$RC.2, levels=c("FS", "HS", "U"))

  tbl <- with(RCTI, Table(RC.1, RC.2))/2  # pairs included double
  tbl["U", "U"] <- nrow(Ped2) * (nrow(Ped2)-1)/2 - sum(tbl)
  tbl
  #  sweep(tbl, 1, rowSums(tbl), "/")
}


#============================================================================
#============================================================================
