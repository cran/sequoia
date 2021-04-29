#' @title Assign Family IDs
#'
#' @description Add a column with family IDs (FIDs) to a pedigree, with each
#'  number denoting a cluster of connected individuals.
#'
#' @details This function repeatedly finds all ancestors and all descendants of
#'  each individual in turn, and ensures they all have the same Family ID. Not
#'  all connected individuals are related, e.g. all grandparents of an
#'  individual will have the same FID, but will typically be unrelated.
#'
#' When UseMaybeRel = TRUE, probable relatives are added to existing family
#'  clusters, or existing family clusters may be linked together. Currently no
#'  additional family clusters are created.
#'
#' @param  Ped dataframe with columns id - parent1 - parent2; only the
#'   first 3 columns will be used.
#' @param  SeqList list as returned by \code{\link{sequoia}}. If 'Ped' is not
#'  provided, the element 'Pedigree' from this list will be used if present,
#'  and element 'Pedigreepar' otherwise.
#' @param  UseMaybeRel use \code{SeqList$MaybeRel}, the dataframe with probable
#' but non-assigned relatives, to assign additional family IDs?
#'
#' @return A dataframe with the provided pedigree, with a column 'FID' added.
#'
#' @export

FindFamilies <- function(Ped=NULL, SeqList=NULL, UseMaybeRel=FALSE) {
  if (is.null(Ped) & is.null(SeqList)) {
    stop("please provide either Ped or SeqList")
  } else if (is.null(Ped)) {
    if (any(names(SeqList)=="Pedigree")) {
      Ped <- SeqList[["Pedigree"]]
    } else if (any(names(SeqList)=="PedigreePar")) {
      Ped <- SeqList[["PedigreePar"]]
    } else {
      stop("please provide either Ped or SeqList with element 'PedigreePar' or 'Pedigree'")
    }
  } else {
    if (!(is.data.frame(Ped) | is.matrix(Ped)) || nrow(Ped)<2 || ncol(Ped)<3) {
      stop("'Ped' must be a dataframe with at least columns id - dam - sire")
    }
  }
  Ped <- Ped[, 1:3]

  Ped$FID <- 0
  for (i in 1:nrow(Ped)) {
    if (Ped$FID[i]!=0 | !is.na(Ped[i,2]) | !is.na(Ped[i,3]))  next
    Ped$FID[i] <- i
    AL <- list()
    DL <- list()
    AL[[1]] <- unique(unlist(GetAncest(i, Ped)))
    DL[[1]] <- unique(unlist(GetDesc(i, Ped)))
    for (x in 1:20) {

      if (length(AL[[x]]) > 0) {
        DXL <- list()
        for (a in seq_along(AL[[x]])) {
          DXL[[a]] <- unique(unlist(GetDesc(which(Ped[,1] == AL[[x]][a]), Ped)))
        }
        DL[[x+1]] <- unique(unlist(DXL))
        if (x>1)  DL[[x+1]] <- setdiff(DL[[x+1]], unlist(c(DL[1:x], AL[1:(x-1)])))
      } else {
        DL[[x+1]] <- numeric()
      }

      if (length(DL[[x]]) > 0) {
        AXL <- list()
        for (d in seq_along(DL[[x]])) {
          AXL[[d]] <- unique(unlist(GetAncest(which(Ped[,1] == DL[[x]][d]), Ped)))
        }
        AL[[x+1]] <- unique(unlist(AXL))
        if (x>1)  AL[[x+1]] <- setdiff(AL[[x+1]], unlist(c(AL[1:x], DL[1:(x-1)])))
      } else {
        AL[[x+1]] <- numeric()
      }
    }
    Ped$FID[Ped[,1] %in% unlist(c(AL, DL))] <- i
  }
  Ped$FID <- as.numeric(factor(Ped$FID))

  if (UseMaybeRel) {
    if (any(names(SeqList)=="MaybeRel")) {
      MR <- SeqList$MaybeRel
    } else {
      MR <- SeqList$MaybePar
    }

    PedX <- Ped[, c("id", "FID")]
    Singles <- as.numeric(names(table(PedX$FID)[table(PedX$FID)==1]))
    PedX$FID[PedX$FID %in% Singles] <- 0
    MR <- merge(MR, stats::setNames(PedX, c("ID1", "FID1")), all.x=TRUE)
    MR <- merge(MR, stats::setNames(PedX, c("ID2", "FID2")), all.x=TRUE)

    NewLinks <- with(MR, MR[which(FID1 != FID2), ])
    if (nrow(NewLinks)>0) {
      for (i in 1:nrow(NewLinks)) {
        if (NewLinks$FID1[i]==0 & NewLinks$FID2[i]!=0) {
          Ped$FID[Ped$id == NewLinks$ID1[i]] <- NewLinks$FID2[i]
        } else if (NewLinks$FID1[i]!=0 & NewLinks$FID2[i]==0) {
          Ped$FID[Ped$id == NewLinks$ID2[i]] <- NewLinks$FID1[i]
        } else {
          newFID <- min(NewLinks[i, c("FID1", "FID2")])
          oldFID <- max(NewLinks[i, c("FID1", "FID2")])
          Ped$FID[which(Ped$FID == oldFID)] <- newFID
        }
      }
    }
  }

  Ped
}


#======================================================================
# get a list with all ancestors of individual on row i of Ped
GetAncest <- function(i, Ped) {
  PL <- list()
  PL[[1]] <- unique(stats::na.exclude(unlist(Ped[i, 2:3])))
  for (g in 1:100) {   # assuming Ped < 100 generations
    if (length(PL[[g]]) == 0)  break
    PL[[g+1]] <- unique(stats::na.exclude(unlist(Ped[Ped[,1] %in% PL[[g]], 2:3])))  # next generation
  }
  PL
}


#======================================================================
# get a list with all descendants of individual on row i of Ped
GetDesc <- function(i, Ped) {
  PL <- list()
  PL[[1]] <- unlist(Ped[which(Ped[,2]==Ped[i,1] | Ped[,3]==Ped[i,1]), 1])
  for (g in 1:100) {   # assuming Ped < 100 generations
    if (length(PL[[g]]) == 0)  break
    PL[[g+1]] <- unlist(Ped[which(Ped[,2] %in% PL[[g]] | Ped[,3] %in% PL[[g]]), 1])
  }
  PL
}


#======================================================================
#' @title Back-transform IDs
#'
#' @description Reverse the joining of FID and IID in
#' \code{\link{GenoConvert}} and \code{\link{LHConvert}}
#'
#' @details Note that the family IDs are the ones provided, and not
#'  automatically updated. New, numeric ones can be obtained with
#'   \code{\link{FindFamilies}}.
#'
#' @param Ped pedigree as returned by sequoia (e.g. \code{SeqOUT$Pedigree}).
#' @param FIDsep characters inbetween FID and IID in composite-ID.
#'
#' @return A pedigree with 6 columns
#' \item{FID}{family ID of focal individual (offspring).}
#' \item{id}{within-family of focal individual}
#' \item{dam.FID}{original family ID of assigned dam}
#' \item{dam}{within-family of dam}
#' \item{sire.FID}{original family ID of assigned sire}
#' \item{sire}{within-family of sire}
#'
#' @export

PedStripFID <- function(Ped, FIDsep="__") {
  PedL <- list()
  for (i in 1:3) {
    PedL[[i]] <- StripFam(Ped[, i], FIDsep=FIDsep)
  }
  stats::setNames(cbind(PedL[[1]], PedL[[2]], PedL[[3]]),
           c("FID", "id", "dam.FID", "dam", "sire.FID", "sire"))
}


#==========================================
# IN: vector with FID_IID
# OUT: data.frame with columns FID, IID
StripFam <- function(V, FIDsep="__") {
  y <- plyr::ldply(strsplit(V, split=FIDsep), function(x) as.data.frame(t(x)))
}
