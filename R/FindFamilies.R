#' @title Assign Family IDs
#'
#' @description Find clusters of connected individuals in a pedigree, and assign
#'   each cluster a unique family ID (FID).
#'
#' @details This function repeatedly finds all ancestors and all descendants of
#'  each individual in turn, and ensures they all have the same Family ID. Not
#'  all connected individuals are related, e.g. all grandparents of an
#'  individual will have the same FID, but will typically be unrelated.
#'
#' When \code{UseMaybeRel = TRUE}, probable relatives are added to existing
#' family clusters, or existing family clusters may be linked together.
#' Currently no additional family clusters are created.
#'
#' @param  Pedigree dataframe with columns id - parent1 - parent2; only the
#'   first 3 columns will be used.
#' @param  SeqList list as returned by \code{\link{sequoia}}. If \code{Pedigree}
#'   is not provided, the element \code{Pedigree} from this list will be used if
#'   present, and element \code{Pedigreepar} otherwise.
#' @param  MaybeRel Output from \code{\link{GetMaybeRel}}, a dataframe with
#'   probable but non-assigned relatives.
#'
#' @return A numeric vector with length equal to the number of unique
#'   individuals in the pedigree (i.e. number of rows in pedigree after running
#'   \code{\link{PedPolish}} on \code{Pedigree}).
#'
#' @seealso \code{\link{GetAncestors}, \link{GetDescendants},
#'   \link{getGenerations}}
#'
#' @export
#'
#' @examples
#'
#' PedG <- SeqOUT_griffin$PedigreePar[,1:3]
#' FID_G <- FindFamilies(PedG)
#' PedG[FID_G==4,]


FindFamilies <- function(Pedigree=NULL, SeqList=NULL, MaybeRel=NULL) {
  if (is.null(Pedigree) & is.null(SeqList)) {
    stop("please provide either Ped or SeqList")
  } else if (is.null(Pedigree)) {
    if (any(names(SeqList)=="Pedigree")) {
      Ped <- SeqList[["Pedigree"]][,1:3]
    } else if (any(names(SeqList)=="PedigreePar")) {
      Ped <- SeqList[["PedigreePar"]][,1:3]
    } else {
      stop("please provide either Ped or SeqList with element 'PedigreePar' or 'Pedigree'")
    }
  } else {
    if (!(is.data.frame(Pedigree) | is.matrix(Pedigree)) || nrow(Pedigree)<2 || ncol(Pedigree)<3) {
      stop("'Ped' must be a dataframe with at least columns id - dam - sire")
    }
    Ped <- PedPolish(Pedigree, KeepAllColumns=FALSE, StopIfInvalid=FALSE)
  }

  FID <- rep(0, nrow(Ped))
  for (i in 1:nrow(Ped)) {
    if (FID[i]!=0 | !is.na(Ped[i,2]) | !is.na(Ped[i,3]))  next
    FID[i] <- i
    AL <- list()
    DL <- list()
    AL[[1]] <- unique(unlist(GetAnc(i, Ped)))
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
          AXL[[d]] <- unique(unlist(GetAnc(which(Ped[,1] == DL[[x]][d]), Ped)))
        }
        AL[[x+1]] <- unique(unlist(AXL))
        if (x>1)  AL[[x+1]] <- setdiff(AL[[x+1]], unlist(c(AL[1:x], DL[1:(x-1)])))
      } else {
        AL[[x+1]] <- numeric()
      }
    }
    FID[Ped[,1] %in% unlist(c(AL, DL))] <- i
  }
  FID <- as.numeric(factor(FID))  # re-number to sequential numbers


  MR <- NULL
  MR_errmsg <- "'MaybeRel' must be the list returned by GetMaybeRel(), or a dataframe with columns 'ID1' and 'ID2'"
  if (is.data.frame(MaybeRel)) {
    if ('ID1' %in% colnames(MaybeRel) & 'ID2' %in% colnames(MaybeRel)) {
      MR <- MaybeRel[, c('ID1', 'ID2')]
    } else {
      stop(MR_errmsg)
    }
  } else if (inherits(MaybeRel, 'list')) {
    if ('MaybeRel' %in% names(MaybeRel)) {
      MR <- MaybeRel$MaybeRel[, c('ID1', 'ID2')]
    } else if ('MaybePar' %in% names(MaybeRel)) {
      MR <- MaybeRel$MaybePar[, c('ID1', 'ID2')]
    } else {
      stop(MR_errmsg)
    }
  } else if (!is.null(MaybeRel)) {
    stop(MR_errmsg)
  }

  if (!is.null(MR)) {
    PedX <- data.frame(id=Ped$id, FID)
    Singles <- as.numeric(names(table(PedX$FID)[table(PedX$FID)==1]))
    PedX$FID[PedX$FID %in% Singles] <- 0
    MR <- merge(MR, stats::setNames(PedX, c("ID1", "FID1")), all.x=TRUE)
    MR <- merge(MR, stats::setNames(PedX, c("ID2", "FID2")), all.x=TRUE)

    NewLinks <- with(MR, MR[which(FID1 != FID2), ])
    if (nrow(NewLinks)>0) {
      for (i in 1:nrow(NewLinks)) {
        if (NewLinks$FID1[i]==0 & NewLinks$FID2[i]!=0) {
          FID[Ped$id == NewLinks$ID1[i]] <- NewLinks$FID2[i]
        } else if (NewLinks$FID1[i]!=0 & NewLinks$FID2[i]==0) {
          FID[Ped$id == NewLinks$ID2[i]] <- NewLinks$FID1[i]
        } else {
          newFID <- min(NewLinks[i, c("FID1", "FID2")])
          oldFID <- max(NewLinks[i, c("FID1", "FID2")])
          FID[which(FID == oldFID)] <- newFID
        }
      }
    }
  }

  return( FID )
}


#======================================================================
# get a list with all ancestors of individual on row i of Ped
GetAnc <- function(i, Ped) {
  PL <- list()
  PL[[1]] <- unique(stats::na.exclude(unlist(Ped[i, 2:3])))
  for (g in 1:100) {   # assuming Ped < 100 generations
    if (length(PL[[g]]) == 0)  break
    PL[[g+1]] <- unique(stats::na.exclude(unlist(Ped[Ped[,1] %in% PL[[g]], 2:3])))  # next generation
  }
  PL
}


#' @title Get ancestors
#'
#' @description get all ancestors of an individual
#'
#' @param id  id of the individual
#' @param Pedigree dataframe with columns id - parent1 - parent2; only the
#'   first 3 columns will be used.
#'
#' @return a list with as first element \code{id}, second parents, third
#'  grandparents, etc.. Each element is a vector with ids, the first three
#'  elements are named, the rest numbered. Ancestors are unsorted within each
#'  list element.
#'
#' @export
#'
#' @examples
#' Anc_i200  <- GetAncestors('i200_2010_F', Ped_griffin)
#'
#'

GetAncestors <- function(id, Pedigree) {
  Ped <- PedPolish(Pedigree, KeepAllColumns=FALSE, StopIfInvalid=FALSE)
  row_i <- which(Ped[,1] == id)
  if (length(row_i)==0)  stop('id not in Pedigree')
  Anc <- GetAnc(row_i, Ped)
  if (any(unlist(Anc) == id)) {
    loopsize <- min(which(sapply(Anc, function(v) id %in% v)))
    cli::cli_alert_warning('Individual {id} is its own ancestor {loopsize} generations back')
  }
  c(list('id' = id,
         'parents' = Anc[[1]],
         'grandparents' = Anc[[2]]),
         Anc[3:(length(Anc)-1)])   # last element is always character(0)
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


#' @title Get descendants
#'
#' @description get all descendants of an individual
#'
#' @param id  id of the individual
#' @param Pedigree dataframe with columns id - parent1 - parent2; only the
#'   first 3 columns will be used.
#'
#' @return a list with as first element \code{id}, second offspring, third
#'  grand-offspring, etc.. Each element is a vector with ids, the first three
#'  elements are named, the rest numbered.
#'
#' @export

GetDescendants <- function(id, Pedigree) {
  Ped <- PedPolish(Pedigree, KeepAllColumns=FALSE, StopIfInvalid=FALSE)
  row_i <- which(Ped[,1] == id)
  if (length(row_i)==0)  stop('id not in Pedigree')
  Desc <- GetDesc(row_i, Ped)
  if (any(unlist(Desc) == id)) {
    loopsize <- min(which(sapply(Desc, function(v) id %in% v)))
    cli::cli_alert_warning('Individual {id} is its own ancestor {loopsize} generations back')
  }
  c(list('id' = id,
         'offspring' = Desc[[1]],
         'grandoffspring' = Desc[[2]]),
    Desc[3:(length(Desc)-1)])
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
