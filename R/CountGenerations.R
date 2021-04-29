#' @title Count Generations
#'
#' @description For each individual in a pedigree, count the number of
#'   generations since its most distant pedigree founder.
#'
#' @param Ped  dataframe, pedigree with the first three columns being id - dam -
#'   sire. Column names are ignored, as are additional columns.
#' @param StopIfInvalid if a pedigree loop is detected, stop with an error
#'   (TRUE, default) or return the Pedigree, to see where the problem(s) occur.
#'
#' @return  A vector with the generation number for each individual, starting at
#'   0 for founders. NA indicates a pedigree loop where an individual is its own
#'   ancestor (or that the pedigree has >1000 generations). Returned invisibly
#'   to be a part of QC.
#'
#' @export

getGenerations <- function(Ped, StopIfInvalid=TRUE) {
  for (p in 2:3) {
    Ped[which(Ped[,p]==0), p] <- NA
  }
  Ped <- as.data.frame(Ped)

  if (!all(na.exclude(unlist(Ped[,2:3])) %in% Ped[,1]))
    stop("Some parents do not occur in first column; please call PedPolish() first")

  Ped$gen <- NA  # individual's generation
  Ped$gen.dam <- NA
  Ped$gen.sire <- NA
  Ped$gen[is.na(Ped[,2]) & is.na(Ped[,3])] <- 0
  for (x in 0:1000) {
    Ped$gen.dam[is.na(Ped$gen.dam) & Ped[,2] %in% Ped[which(Ped$gen<=x), 1]] <- x
    Ped$gen.sire[is.na(Ped$gen.sire) & Ped[,3] %in% Ped[which(Ped$gen<=x), 1]] <- x
    Ped$gen[which(is.na(Ped$gen) &
                   (Ped$gen.dam<=x | is.na(Ped[,2])) &
                   (Ped$gen.sire<=x | is.na(Ped[,3])))] <- x+1
    if (!any(is.na(Ped$gen)))  break
  }
  if (any(is.na(Ped$gen))) {
    msg <- "An individual is its own ancestor, or >1000 generations."
    if (StopIfInvalid) {
      stop(msg,
           "\n Use getGenerations(, StopIfInvalid=FALSE) to find the culprit",
           call.=FALSE)
    } else {
      warning(msg, immediate.=TRUE)
    }
  }
  invisible( setNames(Ped$gen, Ped[,1]) )
}
