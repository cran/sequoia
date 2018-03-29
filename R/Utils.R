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
#' @title write data to a file column-wise
#'
#' @description write data.frame or matrix to a text file, using white
#' space padding to keep columns aligned as in \code{print}
#'
#' @param x the object to be written, preferably a matrix or data frame.
#'  If not, it is attempted to coerce x to a matrix.
#' @param file a character string naming a file.
#' @param row.names a logical value indicating whether the row names of x are
#'  to be written along with x.
#' @param col.names a logical value indicating whether the column names of x
#'  are to be written along with x
#'
#' @export
#'
writeColumns <- function(x, file="", row.names=TRUE,
                         col.names=TRUE) {
  M <- as.matrix(x)
  if(col.names)  M <- rbind(colnames(x), M)
  if (row.names) {
    if (col.names) M <- cbind(c("", rownames(x)), M)
    if (!col.names) M <- cbind(rownames(x), M)
  }
  write(t(format(M)), ncolumns = ncol(M), file = file)
}


#======================================================================
#' @title special merge
#'
#' @description As regular merge, but combine data from columns with the same
#'  name
#'
#' @param df1  first dataframe (lowest priority if overwrite=TRUE)
#' @param df2  second dataframe (highest priority if overwrite=TRUE)
#' @param by  columns used for merging, required.
#' @param overwrite  If FALSE (the default), NA's in df1 are replaced by
#'   values from df2. If TRUE, all values in df1 are overwritten by values from
#'   df2, except where df2 has NA.
#' @param ...  additional arguments to merge, such as \code{all}.
#'
#' @export

Merge <- function(df1, df2, by, overwrite=FALSE, ...) {
  commonNames <- names(df1)[which(colnames(df1) %in% colnames(df2))]
  commonNames <- commonNames[!commonNames %in% by]
  dfmerged <- merge(df1,df2,by=by,...)
  for(i in commonNames){
    left <- paste0(i, ".x")
    right <- paste0(i, ".y")
    if (!overwrite) {
      dfmerged[is.na(dfmerged[left]),left] <- dfmerged[is.na(dfmerged[left]),right]
    } else {
      dfmerged[!is.na(dfmerged[right]),left] <- dfmerged[!is.na(dfmerged[right]),right]
    }
    dfmerged[right]<- NULL
    colnames(dfmerged)[colnames(dfmerged) == left] <- i
  }
  dfmerged
}


#======================================================================
# table, sets UseNA to 'ifany'
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
eqv <- function(x, V, xNA=NA) {
  if (length(x)==1) {
    if (!is.na(x)) {
      y <- ifelse(!is.na(V), x==V, FALSE)
    } else if (is.na(xNA)) {
      y <- is.na(V)
    } else {
      y <- rep(xNA, length(V))
    }
  } else if (length(x)==length(V)) {
    y <- ifelse(!is.na(x) & !is.na(V),
              x==V,
              ifelse(is.na(x) & is.na(V),
                     ifelse(is.na(xNA),
                          TRUE,
                          xNA),
                     FALSE))
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
    names(AddPed)[1:3] <- names(PedIN)[1:3]
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




