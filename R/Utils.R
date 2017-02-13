# Various miscelaneous functions

#======================================================================
#' Convert factor to numeric
#'
#' Converts factors to their numeric values, via \code{as.character}.
#'
#' @param x a factor
#' @return A numeric vector of the same length as x.
FacToNum <- function(x) as.numeric(as.character(x))


#======================================================================
#' transform vector to matrix
#'
#'
#' @param V a vector
#' @param nr number of non-empty rows in matrix
#' @param nc number of columns in matrix
#' @param Ng_odd boolean, apply correction when number of genotyped individuals
#'   is odd
#'
#' @return A n x ncol matrix
#'
VtoM <- function(V, nr=NULL, nc=2, Ng_odd=FALSE) {
  if(Ng_odd) {
    V <- V[1:((length(V)/nc-1)*nc)]
  }
  M <- matrix(V, length(V)/nc, nc)
  if(!is.null(nr)) M <- M[1:nr, , drop=FALSE]
  M
}


#======================================================================
#' Data input
#'
#' Read data with header row.
#'
#' @param ... parameters for read.table
#' @param sep column seperator
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
#' Value Matching
#'
#' Like \code{\%in\%}, returns a logical vector indicating if there is a match
#'  or not for its left operand, but returns NA for NA's in the left operand.
#'
#' @param x vector: the values to be matched
#' @param y vector: the values to be matched against
#'
#' @return A logical vector of the same length as x.
#'
# #' @examples
# #' X <- c(1:5, NA, NA)
# #' Y <- c(3:10)
# #' X %in% Y
# #' table(X %in% Y, useNA="ifany")
# #' X %ina% Y
# #' table(X %ina% Y, useNA="ifany")
"%ina%" <- function(x, y) ifelse(!is.na(x), match(x, y, nomatch = 0) > 0, NA)


#======================================================================
#' paste directory and file name
#'
#' @param folder foldername
#' @param fileName filename

pasteD <- function(folder, fileName) paste(folder, fileName, sep="/")


#======================================================================
#' create a table, and ensure that the levels TRUE, FALSE and NA are always all
#'  represented
#'
#' @param x  a logical vector

tbl.logic <- function(x) table(factor(x, levels=c(TRUE, FALSE, NA)),
                               useNA="always")


#======================================================================
#' Comparison
#'
#' Test which elements in a vector are equal to x, returning FALSE at missing
#' values in V
#'
#' @param x  a value
#' @param V  a vector of the same type as x
#'
#' @return a logical vector of the same length as v
#'
eqv <- function(x, V) {
  if (!is.na(x)) {
    y <- ifelse(!is.na(V), x==V, FALSE)
  } else {
    y <- logical(length=length(V))
  }
  y
}


#======================================================================
#' function in Examples from integer {base}
#'
#' @param x a number
#' @param tol tolerance
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol



#======================================================================
#' Partial pedigree fix
#'
#' Add rows for all parents in pedigree who are not yet represented in the
#' first column
#'
#' @param PedIN  A dataframe with IDs in the first column and parents in the
#'   second and third columns
#'
#' @return A dataframe with one row per individual occuring in the pedigree


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
    names(AddPed) <- names(PedIN)[1:3]
    Ped <- rbind(AddPed, PedIN[,1:3])  # presume ancestors
  }
  Ped
}
