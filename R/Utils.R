# Various miscelaneous functions

#======================================================================
# Convert factor to numeric
FacToNum <- function(x) as.numeric(as.character(x))


#======================================================================
# transform vector to matrix
VtoM <- function(V, nr=NULL, nc=2, Ng_odd=FALSE) {
  if(Ng_odd) {
    V <- V[1:((length(V)/nc-1)*nc)]
  }
  M <- matrix(V, length(V)/nc, nc)
  if(!is.null(nr)) M <- M[1:nr, , drop=FALSE]
  M
}


#======================================================================
# Data input
ReadTable <- function(..., sep="\t") utils::read.table(...,
                                                header=TRUE,
                                                stringsAsFactors=FALSE,
                                                sep=sep,
                                                na.strings=c("", NA))


#======================================================================
#' table
#'
#' sets UseNA to 'ifany'.
#'
#' @param ... one or more objects which can be interpreted as factors
#'  (including character strings), or a list (or data frame) whose components
#'   can be so interpreted.
Table <- function(...) table(..., useNA="ifany")


#======================================================================
# Value Matching
"%ina%" <- function(x, y) ifelse(!is.na(x), match(x, y, nomatch = 0) > 0, NA)


#======================================================================
# paste directory and file name
pasteD <- function(folder, fileName) paste(folder, fileName, sep="/")


#======================================================================
# create a table, and ensure that the levels TRUE, FALSE and NA are always all
tbl.logic <- function(x) table(factor(x, levels=c(TRUE, FALSE, NA)),
                               useNA="always")


#======================================================================
# Comparison
eqv <- function(x, V, xNA=FALSE) {
  if (length(x)==1) {
    if (!is.na(x)) {
      y <- ifelse(!is.na(V), x==V, FALSE)
    } else {
      y <- rep(xNA, length(V))
    }
  } else if (length(x)==length(V)) {
    y <- ifelse(!is.na(x) & !is.na(V), x==V,
              ifelse(is.na(x) & is.na(V), xNA, FALSE))
  } else {
    stop("unequal lengths")
  }
  y
}


#======================================================================
# function in Examples from integer {base}
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#======================================================================
# Partial pedigree fix
AddParPed <- function(PedIN) {
  Ped <- unique(PedIN)
  UID <- unique(c(as.character(Ped[,1]),
                  as.character(Ped[,2]),
                  as.character(Ped[,3])))
  UID <- stats::na.exclude(UID)
  if (length(UID) > nrow(Ped)) {
    AddPed <- data.frame(id=setdiff(UID, Ped[,1]),
                         dam=NA,
                         sire=NA,
                         stringsAsFactors=FALSE)
    Ped <- merge(AddPed, PedIN, all=TRUE)  # presume ancestors
  }
  Ped
}


#===============================
# replace numbers by names
NumToID <- function(x, k=NULL, GenoNames=NULL, dumID=NULL) {
  type <- ifelse(x>0, "r", ifelse(x<0, "d", "n"))
  ID <- rep(NA, length(x))
  ID[type=="r"] <- GenoNames[x[type=="r"]]
  ID[type=="d"] <- dumID[-x[type=="d"], k]
  ID
}


IDToNum <- function(NamePed, GenoNames) {  # GenoNames = rownames(GenoM)
  names(NamePed) <- c("id", "dam", "sire")
  if (!all(GenoNames %in% NamePed$id)) {
    NamePed <- rbind(NamePed, data.frame(id = GenoNames[!GenoNames %in% NamePed$id],
                                         dam = NA, sire = NA))
  }
  rownames(NamePed) <- NamePed$id
  NamePed <- NamePed[GenoNames, ]
  GenoNums <- 1:length(GenoNames)
  names(GenoNums) <- GenoNames
  NumPed <- data.frame(id = GenoNums,
                       dam = GenoNums[NamePed$dam],
                       sire = GenoNums[NamePed$sire])
  for(i in 2:3) NumPed[is.na(NumPed[,i]), i] <- 0
  NumPed
}




