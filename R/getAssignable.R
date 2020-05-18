#' @title Assignability of reference pedigree
#'
#' @description Identify which individuals are genotyped, and which can
#'   potentially be substituted by a dummy individual. 'Dummifiable' are those
#'   non-genotyped individuals with at least 2 genotyped offspring, or at least
#'   1 genotyped offspring and 1 genotyped parent.
#'
#' @details It is assumed that all individuals in \code{Genotyped} have been
#'   genotyped for a sufficient number of SNPs. To identify samples with a
#'   too-low call rate, use \code{\link{CheckGeno}}. To calculate the call rate
#'   for all samples, see the examples below.
#'
#'   Some parents indicated here as assignable may never be assigned by sequoia,
#'   for example parent-offspring pairs where it cannot be determined which is
#'   the older of the two, or grandparents that are indistinguishable from full
#'   avuncular (i.e. genetics inconclusive because the candidate has no parent
#'   assigned, and ageprior inconclusive).
#'
#' @param Pedigree  dataframe with columns id-dam-sire. Reference pedigree.
#' @param Genotyped  character vector with ids of genotyped individuals.
#'
#' @return the \code{Pedigree} dataframe with 2 additional columns,
#'   \code{dam.cat} and \code{sire.cat}, with coding similar to
#'   that used by \code{\link{PedCompare}}:
#' \item{GG}{Genotyped individual, genotyped parent}
#' \item{GD}{Genotyped individual, Dummy parent; i.e. 'id' has at least 1
#' genotyped sibling or a genotyped grandparent}
#' \item{DG}{Dummy individual, Genotyped parent; i.e. 'id' has at least 1
#' genotyped offspring, and parent is assignable as grandparent of the
#' dummy-substituted-individual's offspring}
#' \item{DD}{Dummy individual, Dummy parent}
#' \item{X}{Either or both id and parent is/are not genotyped, and has/have no
#' genotyped offspring, and therefore the parent- offspring link cannot be
#' assigned.}
#' \item{NA}{No parent in Pedigree}
#'
#'
#' @examples
#' data(Ped_HSg5, SimGeno_example, package="sequoia")
#' PedA <- getAssignCat(Ped_HSg5, rownames(SimGeno_example))
#' table(PedA$dam.cat, PedA$sire.cat, useNA="ifany")
#'
#' # calculate call rate
#' \dontrun{
#' CallRates <- apply(MyGenotypes, MARGIN=1,
#'                    FUN = function(x) sum(x!=-9)) / ncol(MyGenotypes)
#' hist(CallRates, breaks=50, col="grey")
#' GoodSamples <- rownames(MyGenotypes)[ CallRates > 0.8]
#' threshold depends on total number of SNPs, genotyping errors, proportion of
#' candidate parents that are SNPd (sibship clustering is more prone to false
#' positives).
#' PedA <- getAssignCat(MyOldPedigree, rownames(GoodSamples))
#' }
#
#' @export

getAssignCat <- function(Pedigree, Genotyped) {
  # check input
  if(is.null(Pedigree)) stop("Must provide 'Pedigree'")
  if (is.null(Genotyped))  stop("Must provide 'Genotyped'")
  Pedigree <- PedPolish(Pedigree[, 1:3], ZeroToNA=TRUE)
  if (length(intersect(Pedigree$id, Genotyped)) == 0)  stop("'Pedigree' and 'Genotyped' have no IDs in common")

  #~~~~~~~~~~~~~~
  Dummifiable <- unlist(GetDummifiable(Pedigree, Genotyped))

  for (x in c("id", "dam", "sire")) {
    Pedigree[, paste0(x, ".cat")] <- ifelse(Pedigree[,x] %in% Genotyped, "G",
                                            ifelse(Pedigree[,x] %in% Dummifiable, "D",
                                                   "X"))
  }
  return( Pedigree )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GetDummifiable <- function(Pedigree, Genotyped) {
  names(Pedigree) <- c("id", "dam", "sire")
  for (x in 1:3)  Pedigree[,x] <- as.character(Pedigree[,x])
  PedG <- Pedigree[Pedigree$id %in% Genotyped, ]
  UniqueParents <- with(PedG, list(dam = unique(na.exclude(dam[!dam %in% Genotyped])),
                             sire = unique(na.exclude(sire[!sire %in% Genotyped]))))
  Dummifiable <- list()
  for (p in c("dam", "sire")) {
    NumOff <- table(factor(PedG[,p], levels=UniqueParents[[p]]))
    PedP <- Pedigree[match(UniqueParents[[p]], Pedigree$id), ]
    HasGP <- PedP$dam %in% Genotyped | PedP$sire %in% Genotyped
    Dummifiable[[p]] <- UniqueParents[[p]][NumOff > 1 | HasGP]
  }

  return( Dummifiable )
}
