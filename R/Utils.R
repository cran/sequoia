# Various miscelaneous functions

#======================================================================
# Convert factor to numeric
FacToNum <- function(x) as.numeric(as.character(x))


#======================================================================
fc <- function(x, w=2)  formatC(x, width=w, flag="0")


#=====================================================================
# for subsetting vector V[1:n] : unexpected behaviour when n=0

s <- function(n)  if(n>0) 1:n else 0

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
# test if can be converted to integers/numbers

check.integer <- function(xx) ifelse(is.na(xx), NA,
                                     grepl("^[-]{0,1}[0-9]{1,}$", xx))
check.numeric <- function(xx) ifelse(is.na(xx), NA,
                                     grepl("^[-]{0,1}[0-9]{0,}.{0,1}[0-9]{1,}$", xx))


#======================================================================
#' @title special merge
#'
#' @description As regular merge, but combine data from columns with the same
#'  name
#'
#' @param df1  first dataframe (lowest priority if overwrite=TRUE)
#' @param df2  second dataframe (highest priority if overwrite=TRUE)
#' @param by  columns used for merging, required.
#' @param overwrite  If FALSE (the default), NA's in df1 are replaced by values
#'   from df2. If TRUE, all values in df1 are overwritten by values from df2,
#'   except where df2 has NA.
#' @param ...  additional arguments to merge, such as \code{all}.
#'
#' @keywords internal

MergeFill <- function(df1, df2, by, overwrite=FALSE, ...) {
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


#===============================
# replace numbers by names
NumToID <- function(x, k=NULL, GenoNames=NULL, dumID=NULL) {
  type <- ifelse(is.na(x), "n",
               ifelse(x>0, "r",
                  ifelse(x<0, "d", "n")))
  ID <- rep(NA, length(x))
  if (any(type=="r"))  ID[type=="r"] <- GenoNames[x[type=="r"]]
  if (any(type=="d"))  ID[type=="d"] <- dumID[-x[type=="d"], k]
  return( ID )
}


#===============================
IDToNum <- function(NamePed, GenoNames, DumPrefix=c("F0", "M0")) {  # GenoNames = rownames(GenoM)
  GenoNums <- seq_along(GenoNames)
  names(GenoNums) <- GenoNames
  DumPrefixNchar <- nchar(DumPrefix[1])

  names(NamePed) <- c("id", "dam", "sire")
  for (x in 1:3)  NamePed[, x] <- as.character(NamePed[, x])
  if (!all(GenoNames %in% NamePed$id)) {
    NamePed <- rbind(NamePed, data.frame(id = GenoNames[!GenoNames %in% NamePed$id],
                                         dam = NA,
                                         sire = NA))
  }
  DumNameToID <- function(x, ncp=1) -as.numeric(substr(x,ncp+1,nchar(x)))

  NumPed <- matrix(0, nrow(NamePed), 3,
                   dimnames=list(NamePed$id, c("id", "dam", "sire")))
  for (x in 1:3) {
    type <- getGDO(NamePed[,x], GenoNames, DumPrefix)
    NumPed[type=="Genotyped", x] <- GenoNums[NamePed[type=="Genotyped", x]]
    if (any(type=="Dummy")) {
      NumPed[type=="Dummy", x] <- DumNameToID(NamePed[type=="Dummy", x], DumPrefixNchar)
    }
  }
#  NumPed <- NumPed[order(NumPed[,"id"] < 0, abs(NumPed[,"id"])), ]  # real than dummy
  return( NumPed )
}


#===============================
getGDO <- function(id, SNPd = NULL, DumPrefix = c("F0", "M0")) {
  GDO <- rep("Genotyped", length(id))
  GDO[is.na(id)] <- "None"
  if (!is.null(DumPrefix)) {
    for (p in 1:2) {
      GDO[substr(id,1,nchar(DumPrefix[p])) == DumPrefix[p]] <- "Dummy"
    }
  }
  if (!is.null(SNPd)) {
    GDO[(!id %in% SNPd) & GDO=="Genotyped"] <- "Observed"
  }
  GDO <- factor(GDO, levels=c("Genotyped", "Dummy", "Observed", "None"), ordered=TRUE)
}


#===============================
Replace <- function(V, old, new) {
  # base function 'replace' with match only replaces first match.
  if (length(old) != length(new))  stop("'old' and 'new' must have same length")
  if (!all(old %in% V))  stop("all 'old' must be in V")
  these <- lapply(seq_along(old), function(x, y=V) which(y == old[x]))
  newr <- rep(new, sapply(these, length))
  replace(V, unlist(these), newr)
}



#===============================

ChkOwnAnc <- function(Ped) {
  Ped <- as.matrix(Ped)
  for (i in 1:nrow(Ped)) {
    focal <- Ped[i, 1]
    Par <- unlist(unique(na.exclude(Ped[i, 2:3])))
    for (g in 1:100) {   # assuming Ped < 100 generations
      if (length(Par) == 0)  break
#      cat(g, Par, "\n")
      if (focal %in% Par) {
        stop("Own Ancestor! ", g, " ", focal, " ", Par)
      }
      Par <- unlist(unique(na.exclude(c(Ped[Ped[,1] %in% Par, 2:3]))))  # next generation
    }
  }
  cat("ped OK", "\n")
}


#===============================
.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
          sep = "", collapse = " ")
}
