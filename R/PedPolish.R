#======================================================================
#' @title Fix Pedigree
#'
#' @description Ensure all parents & all genotyped individuals are included,
#'   remove duplicates, rename columns, and replace 0 by NA or v.v..
#'
#' @param Pedigree  dataframe where the first 3 columns are id, dam, sire.
#' @param gID  character vector with ids of genotyped individuals
#'   (rownames of genotype matrix).
#' @param ZeroToNA  logical, replace 0's for missing values by NA's (defaults to
#'   \code{TRUE}).
#' @param NAToZero  logical, replace NA's for missing values by 0's. If
#'   \code{TRUE}, ZeroToNA is automatically set to \code{FALSE}.
#' @param DropNonSNPd  logical, remove any non-genotyped individuals (but keep
#'   non-genotyped parents), & sort pedigree in order of \code{gID}.
#' @param FillParents logical, for individuals with only 1 parent assigned, set
#'   the other parent to a dummy (without assigning siblings or grandparents).
#'   Makes the pedigree compatible with R packages and software that requires
#'   individuals to have either 2 or 0 parents, such as
#'   \code{\link[kinship2]{kinship}}.
#' @param NullOK  logical, is it OK for Ped to be NULL? Then NULL will be
#'   returned.
#' @param LoopCheck  logical, check for invalid pedigree loops by calling
#'   \code{\link{getGenerations}}.
#' @param StopIfInvalid if a pedigree loop is detected, stop with an error
#'   (TRUE, default).
#'
#' @details Recognized column names are any that contain:
#'  \describe{
#'   \item{dam}{"dam", "mother", "mot", "mom", "mum", "mat"}
#'   \item{sire}{"sire", "father", "fat", "dad", "pat"}}
#' \code{sequoia} requires the column order id - dam - sire; columns 2 and 3 are
#' swapped if necessary.
#'
#' @examples
#' \dontrun{
#' # To get the output pedigree into kinship2 compatible format:
#' PedP <- sequoia::PedPolish(SeqOUT$Pedigree, DropNonSNPd=FALSE,
#'                            FillParents = TRUE)
#' PedP$Sex <- with(PedP, ifelse(id %in% dam, "female",  "male"))
#' # default to 'male' to avoid warning: "More than 25% of the gender values are
#' #  'unknown'"
#'
#' Ped.fix <- with(PedP, kinship2::fixParents(id=id, dadid=sire, momid=dam,
#'                                            sex=Sex))
#' Ped.k <- with(Ped.fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))
#' }
#'
#' @export

PedPolish <- function(Pedigree,
                      gID=NULL,
                      ZeroToNA=TRUE,
                      NAToZero=FALSE,
                      DropNonSNPd = TRUE,
                      FillParents = FALSE,
                      NullOK = FALSE,
                      LoopCheck = TRUE,
                      StopIfInvalid = TRUE)
{
  PedName <- as.list(match.call())[["Pedigree"]]
  if (is.null(Pedigree)) {
    if (NullOK) {
      return()
    } else {
      stop("Must provide ", PedName, call. = FALSE)
    }
  }
  if (NAToZero)  ZeroToNA <- FALSE
  #   if (ZeroToNA & NAToZero)  stop("ZeroToNA and NAToZero can't both be TRUE")
  if (!class(Pedigree) %in% c("data.frame", "matrix"))
    stop(PedName, " should be a dataframe or matrix", call. = FALSE)
  if (ncol(Pedigree) < 3)
    stop(PedName, " should be a dataframe with at least 3 columns (id - dam - sire)",
         call. = FALSE)
  Ped <- as.data.frame(unique(Pedigree))

  DamNames <- c("dam", "mother", "mot", "mom", "mum", "mat")
  SireNames <-  c("sire", "father", "fat", "dad", "pat")
  IdNames <- c("id", "iid", "off")
  IdCol <- unlist(sapply(IdNames, function(x) grep(x, names(Ped), ignore.case=TRUE)))
  DamCol <- unlist(sapply(DamNames, function(x) grep(x, names(Ped), ignore.case=TRUE)))
  SireCol <- unlist(sapply(SireNames, function(x) grep(x, names(Ped), ignore.case=TRUE)))
  if ((length(DamCol)==0 | length(SireCol)==0) ||
      (length(IdCol)==0 & (DamCol[[1]]==1 | SireCol[[1]]==1))) {
    stop("Pedigree column names not recognized. Must be id - dam - sire",
         call. = FALSE)
  } else {
    Ped <- Ped[, c(IdCol[[1]], DamCol[[1]], SireCol[[1]],
                   setdiff(seq_len(ncol(Ped)), c(IdCol, DamCol, SireCol)))]
  }

  if (!is.null(gID)) {
    n.shared.ids <- length(intersect(Ped[,1], as.character(gID)))
    if (n.shared.ids==0) {
      stop("GenoM and ", PedName, " do not share any common individuals", call. = FALSE)
    } else if (n.shared.ids < length(gID)/10 && n.shared.ids < nrow(Ped)) {
      warning("GenoM and ", PedName, " share few common individuals", call. = FALSE)
    }
  }

  names(Ped)[1:3] <- c("id", "dam", "sire")
  for (x in c("id", "dam", "sire")) {
    if (any(Ped[,x]=="", na.rm=TRUE)) {  # common issue with read.table
      Ped[which(Ped[,x]==""), x] <- NA
    }
    if (ZeroToNA) Ped[which(Ped[,x]==0), x] <- NA
    Ped[,x] <- as.character(Ped[,x])
    if (NAToZero) Ped[is.na(Ped[,x]), x] <- 0
  }
  Ped <- unique(Ped[!is.na(Ped[,1]), ])
  UID <- stats::na.exclude(unique(c(unlist(Ped[,1:3]),
                                    gID)))
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

  if (LoopCheck)  getGenerations(Ped, StopIfInvalid)   # checks that no individual is its own ancestor

  if (DropNonSNPd && !is.null(gID))  Ped <- Ped[match(gID, Ped$id), ]

  return( Ped )
}

#============================================================================
#============================================================================
