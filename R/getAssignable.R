#' @title Assignability of Reference Pedigree
#'
#' @description Identify which individuals are SNP genotyped, and which can
#'   potentially be substituted by a dummy individual ('Dummifiable').
#'
#' @details It is assumed that all individuals in \code{SNPd} have been
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
#' @param SNPd  character vector with ids of genotyped individuals.
#' @param minSibSize  minimum requirements to be considered 'dummifiable':
#'   \itemize{
#'      \item '1sib' : sibship of size 1, i.e. the non-genotyped individual has
#'        at least 1 genotyped offspring. If there is no sibship-grandparent
#'        this isn't really a sibship, but can be useful in some situations.
#'        Used by \code{\link{CalcOHLLR}}.
#'      \item '1sib1GP': sibship of size 1 with at least 1 genotyped
#'        grandparent. The minimum to be potentially assignable by
#'        \code{\link{sequoia}}.
#'      \item '2sib': at least 2 siblings, with or without grandparents. Used
#'         by \code{\link{PedCompare}}.
#'  }.
#'
#' @return The \code{Pedigree} dataframe with 3 additional columns,
#'   \code{id.cat}, \code{dam.cat} and \code{sire.cat}, with coding similar to
#'   that used by \code{\link{PedCompare}}:
#' \item{G}{Genotyped}
#' \item{D}{Dummy or 'dummifiable'}
#' \item{X}{Not genotyped and not dummifiable, or no parent in pedigree}
#'
#' @examples
#' data(Ped_HSg5, SimGeno_example, package="sequoia")
#' PedA <- getAssignCat(Ped_HSg5, rownames(SimGeno_example))
#' tail(PedA)
#' table(PedA$dam.cat, PedA$sire.cat, useNA="ifany")
#'
#' # calculate call rate
#' \dontrun{
#' CallRates <- apply(MyGenotypes, MARGIN=1,
#'                    FUN = function(x) sum(x!=-9)) / ncol(MyGenotypes)
#' hist(CallRates, breaks=50, col="grey")
#' GoodSamples <- rownames(MyGenotypes)[ CallRates > 0.8]
#' # threshold depends on total number of SNPs, genotyping errors, proportion
#' # of candidate parents that are SNPd (sibship clustering is more prone to
#' # false positives).
#' PedA <- getAssignCat(MyOldPedigree, rownames(GoodSamples))
#' }
#
#' @export

getAssignCat <- function(Pedigree, SNPd, minSibSize = "1sib1GP") {
  # check input
  if (is.null(SNPd))  stop("Must provide 'SNPd'")
  Pedigree <- PedPolish(Pedigree, ZeroToNA=TRUE, NullOK = FALSE, StopIfInvalid=FALSE)
  if (length(intersect(Pedigree$id, SNPd)) == 0)
    stop("'Pedigree' and 'SNPd' have no IDs in common")

  #~~~~~~~~~~~~~~
  Dummifiable <- unlist(GetDummifiable(Pedigree, SNPd, minSibSize))

  for (x in c("id", "dam", "sire")) {
    Pedigree[, paste0(x, ".cat")] <- ifelse(Pedigree[,x] %in% SNPd, "G",
                                            ifelse(Pedigree[,x] %in% Dummifiable, "D",
                                                   "X"))
  }
  return( Pedigree )
}



#=============================================================================
#=============================================================================
#' @title Dummifiable IDs
#'
#' @description Get the dummifiable individuals, using various possible
#'   criteria
#'
#' @param Pedigree dataframe with id - dam - sire.
#' @param gID  vector with IDs of SNP-genotyped individuals.
#' @param minSibSize  minimum requirements to be considered dummifiable:
#'   \itemize{
#'      \item '1sib' : sibship of size 1, with or without grandparents. The
#'      latter aren't really a sibship, but can be useful in some situations.
#'      \item '1sib1GP': sibship of size 1 with at least 1 grandparent
#'      \item '2sib': at least 2 siblings, with or without grandparents. Used
#'         by \code{\link{PedCompare}}
#'  }.
#'
#' @return A length-2 list (dams, sires) with each element a vector with
#'   dummifiable ids
#'
#' @details  values of minSibSize used by calling functions
#'   \describe{
#'     \item{1sib}{CalcOHLLR, CalcPairLL}
#'     \item{1sib1GP}{getAssignCat (default when user called)}
#'     \item{2sib}{PedCompare}
#'   }
#'
#' @keywords internal

GetDummifiable <- function(Pedigree, gID, minSibSize) {
  names(Pedigree) <- c("id", "dam", "sire")
  if (any(Pedigree == 0, na.rm=TRUE)) {
    for (x in 1:3)  Pedigree[which(Pedigree[,x]==0), x] <- NA
  }
  for (x in 1:3)  Pedigree[,x] <- as.character(Pedigree[,x])
  PedG <- Pedigree[Pedigree$id %in% gID, ]

  UniqueParents <- with(PedG, list(dam = unique(na.exclude(dam[!dam %in% gID])),
                                   sire = unique(na.exclude(sire[!sire %in% gID]))))

  Dummifiable <- list(dam = c(), sire=c())
  for (p in c("dam", "sire")) {
    NumOff <- table(factor(PedG[,p], levels=UniqueParents[[p]]))
    PedP <- Pedigree[match(UniqueParents[[p]], Pedigree$id), ]
    HasGP <- !is.na(PedP$dam) | !is.na(PedP$sire)
    HasGP.G <- PedP$dam %in% gID | PedP$sire %in% gID

    if (minSibSize == "1sib") {
      # to be replaced by fake dummies; CalcOHLLR etc.
      # TODO: count dummifiable offspring too, i.e. make iterative.
      Dummifiable[[p]] <- UniqueParents[[p]]  #[NumOff > 1 | HasGP]
    } else if (minSibSize == "1sib1GP") {
      # potentially assignable by sequoia
      Dummifiable[[p]] <- UniqueParents[[p]][NumOff > 1 | HasGP.G]
    } else if (minSibSize == "2sib") {
      # matcheable by PedCompare
      Dummifiable[[p]] <- UniqueParents[[p]][NumOff > 1]
    } else {
      stop("'minSibSize' must be '1sib', '2sib', or '1sib1GP'")
    }
  }

  return( Dummifiable )
}
