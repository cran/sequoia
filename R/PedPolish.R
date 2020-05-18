#======================================================================
#' @title Pedigree fix
#'
#' @description Ensure all parents & all genotyped individuals are included,
#'   remove duplicates, rename columns, and replace 0 by NA or v.v.
#'
#' @param Ped  dataframe where the first 3 columns are id, dam, sire
#' @param GenoNames  character vector with ids of genotyped individuals
#'   (rownames of genotype matrix)
#' @param ZeroToNA  logical, replace 0's for missing values by NA's (defaults to
#'   \code{TRUE})
#' @param NAToZero  logical, replace NA's for missing values by 0's. If
#'   \code{TRUE}, ZeroToNA is automatically set to \code{FALSE}
#' @param DropNonSNPd  logical, remove any non-genotyped individuals (but keep
#'   non-genotyped parents), & sort pedigree in order of \code{GenoNames}
#' @param FillParents logical, for individuals with only 1 parent assigned, set
#'   the other parent to a dummy (without assigning siblings or grandparents).
#'   Makes the pedigree compatible with R packages and software that requires
#'   individuals to have either 2 or 0 parents, such as
#'   \code{\link[kinship2]{kinship}}.
#'
#' @details recognized column names are any that contain:
#'  \describe{
#'   \item{dam}{"dam", "mother", "mot", "mom", "mum", "mat"}
#'   \item{sire}{"sire", "father", "fat", "dad", "pat"}}
#' \code{sequoia} requires the column order id - dam - sire; columns 2 and 3 are
#' swapped if necessary.
#'
#' @export

PedPolish <- function(Ped,
                      GenoNames=NULL,
                      ZeroToNA=TRUE,
                      NAToZero=FALSE,
                      DropNonSNPd = TRUE,
                      FillParents = FALSE) {
  if (NAToZero)  ZeroToNA <- FALSE
#   if (ZeroToNA & NAToZero)  stop("ZeroToNA and NAToZero can't both be TRUE")
  Ped <- as.data.frame(unique(Ped))
  DamNames <- c("dam", "mother", "mot", "mom", "mum", "mat")
  SireNames <-  c("sire", "father", "fat", "dad", "pat")
  DamCol <- unlist(sapply(DamNames, function(x) grep(x, names(Ped), ignore.case=TRUE)))
  SireCol <- unlist(sapply(SireNames, function(x) grep(x, names(Ped), ignore.case=TRUE)))
  if (length(DamCol)==0 | length(SireCol)==0) {
    stop("Pedigree column names not recognized. Must be id - dam - sire")
  } else if (DamCol[[1]] ==3 & SireCol[[1]] == 2) {
    if (ncol(Ped)==3) {
      Ped <- Ped[, c(1, 3, 2)]
    } else {
      Ped <- Ped[, c(1, 3, 2, 4:ncol(Ped))]
    }
    message("Changed pedigree column order")
  }

  names(Ped)[1:3] <- c("id", "dam", "sire")
  for (x in c("id", "dam", "sire")) {
    if (ZeroToNA) Ped[which(Ped[,x]==0), x] <- NA
    Ped[,x] <- as.character(Ped[,x])
    if (NAToZero) Ped[is.na(Ped[,x]), x] <- 0
  }
  Ped <- unique(Ped[!is.na(Ped[,1]), ])
  UID <- stats::na.exclude(unique(c(unlist(Ped[,1:3]),
                                    GenoNames)))
  if (NAToZero) UID <- UID[UID != 0]

  if (length(UID) > nrow(Ped)) {
    Ped <- merge(data.frame(id = setdiff(UID, Ped$id),
                            dam = ifelse(NAToZero, 0, NA),
                            sire = ifelse(NAToZero, 0, NA),
                            stringsAsFactors=FALSE),
                 Ped,
                 all = TRUE)
  }

  if (FillParents) {
    PP <- c("dam", "sire")
    for (p in 1:2) {
      NeedsPar <- is.na(Ped[, PP[p]]) & !is.na(Ped[, PP[3-p]])
      if (any(NeedsPar)) {
        NewPar <- paste0("X", c("F","M")[p], formatC(1:sum(NeedsPar), width=4, flag=0))
        Ped[NeedsPar, PP[p]] <- NewPar
        Ped <- merge(data.frame(id = NewPar,
                                dam = ifelse(NAToZero, 0, NA),
                                sire = ifelse(NAToZero, 0, NA),
                                stringsAsFactors=FALSE),
                     Ped,
                     all = TRUE)
      }
    }
  }

  if (DropNonSNPd && !is.null(GenoNames))  Ped <- Ped[match(GenoNames, Ped$id), ]

  return( Ped )
}

#============================================================================
#============================================================================
