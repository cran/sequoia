#=======================================================================
#' @title Check LifeHistData
#'
#' @description Check that the provided LifeHistData is in the correct format.
#'
#' @param LifeHistData the dataframe with ID - Sex - Birth year, and optionally
#'   BY.min - BY.max.
#' @param gID  character vector with names of genotyped individuals, i.e.
#'   rownames(GenoM).
#' @param sorted  logical, return lifehistdata for genotyped individuals only,
#'   in strictly the same order. Will including padding with 'empty' rows if an
#'   individual in gID was not in the input-LH.
#' @param returnDups  logical, instead of just the (sorted) LifeHistData, return
#'   a list that also includes a dataframe with duplicate entries and/or a
#'   character vector with genotyped IDs not occuring in LifeHistData (as
#'   formerly returned by \code{DuplicateCheck}).
#'
#' @return A dataframe with LifeHistData formatted for use by the Fortran
#'    part of the program, or a list with duplicate and missing entries.
#'
#' @keywords internal

CheckLH <- function(LifeHistData, gID = NA, sorted=TRUE, returnDups = FALSE)
{
  if (is.null(gID))  gID <- NA
  if (is.null(LifeHistData)) {
    LH.OUT <- data.frame(ID = gID,
                       Sex = 3,
                       BirthYear = -999,
                       BY.min = -999,
                       BY.max = -999,
                       stringsAsFactors = FALSE)
    if (returnDups) {
      return( list(DupLifeHistID = NULL, NoLH = NULL, LifeHistData = LH.OUT) )
    } else {
      return( LH.OUT )
    }
  }
  if (all(is.na(LifeHistData)))
    stop("invalid value for LifeHistData, provide NULL or dataframe",
         call.=FALSE)

  # check if LH ID's match genotype data
  names(LifeHistData)[1:3] <- c("ID", "Sex", "BirthYear")
  LifeHistData$ID <- as.character(LifeHistData$ID)
  if (!all(is.na(gID))) {
    gID <- as.character(gID)
    if (length(intersect(LifeHistData$ID, gID))==0)
      stop("None of the genotyped individuals included in lifehistory data",
           call.=FALSE)
  } else {
    if (sorted)  stop("if sorted=TRUE, gID cannot be NA", call.=FALSE)
  }


  # check/add BY.min & BY.max columns ----
  if (ncol(LifeHistData)==3) {
    LifeHistData$BY.min <- NA
    LifeHistData$BY.max <- NA
  } else if (any(colnames(LifeHistData) %in% c("BY.min", "BY.max"))) {
    if (colnames(LifeHistData)[4] == "BY.min" & colnames(LifeHistData)[5] == "BY.max") {
      LifeHistData <- LifeHistData[, 1:5]
    } else {
      if (any(colnames(LifeHistData) == "BY.min")) {
        BYmin <- LifeHistData$BY.min
      } else {
        BYmin <- NA
      }
      if (any(colnames(LifeHistData) == "BY.min")) {
        BYmax <- LifeHistData$BY.max
      } else {
        BYmax <- NA
      }
      LifeHistData <- LifeHistData[, 1:3]
      LifeHistData$BY.min <- BYmin
      LifeHistData$BY.max <- BYmax
    }
  } else if (ncol(LifeHistData)==5) {
    LifeHistData <- LifeHistData[, 1:5]
    names(LifeHistData)[4:5] <- c("BY.min", "BY.max")
  } else {
    stop("Confused by columns in LifeHistData", call.=FALSE)
  }
  if (any(LifeHistData$BY.max < LifeHistData$BY.min, na.rm=TRUE)) {
    stop("'BY.max' must be greater than or equal to 'BY.min', or NA", call.=FALSE)
  }


  # check if column entries are valid ----
  if (any(grepl(" ", LifeHistData$ID))) {
    stop("LifeHistData IDs must not include spaces", call.=FALSE)
  }
  for (x in c("Sex", "BirthYear", "BY.min", "BY.max")) {
    IsInt <- check.integer(LifeHistData[,x])
    if (any(!IsInt, na.rm=TRUE)) {
#      if (sum(!IsInt, na.rm=TRUE) > sum(!is.na(LifeHistData[,x]))/2) {
        stop("LifeHistData column ", x, " should be integers (whole numbers)",
             call.=FALSE)
#      } else {
#        warning("Converting all values in LifeHistData column ", x, " to integers",
#                immediate. = TRUE)
#      }
    }
    LifeHistData[, x] <- ifelse(IsInt, suppressWarnings(as.integer(as.character(LifeHistData[, x]))), NA)
    if (x=="Sex") {
      if ((sum(LifeHistData$Sex %in% c(1,2,4)) < nrow(LifeHistData)/10) &
          (sum(is.na(LifeHistData$Sex)) + sum(LifeHistData$Sex==3, na.rm=TRUE) < nrow(LifeHistData)/2)) {
        stop("LifeHistData column 2 should contain Sex coded as 1=female, 2=male, 3/NA=unknown, 4=hermaphrodite",
             call.=FALSE)
      }
      LifeHistData$Sex[is.na(LifeHistData$Sex)] <- 3
      LifeHistData$Sex[!LifeHistData$Sex %in% 1:4] <- 3
    } else {
      LifeHistData[is.na(LifeHistData[,x]), x] <- -999
      LifeHistData[LifeHistData[,x] < 0, x] <- -999
    }
  }

  # duplicate check ----
  LifeHistData <- unique(LifeHistData)    # doesn't matter if both entries are identical
  OUT <- list(DupLifeHistID = NULL, NoLH = NULL)
  if (any(duplicated(LifeHistData$ID))) {
    if (returnDups) {
      r1 <- which(duplicated(LifeHistData$ID))
      r2 <- which(duplicated(LifeHistData$ID, fromLast=TRUE))
      OUT$DupLifeHistID <- with(LifeHistData,
                                data.frame(row1 = r1,
                                           row2 = r2,
                                           ID = ID[r1],
                                           Sex1 = ID[r1],
                                           Sex2 = ID[r2],
                                           BirthYear1 = BirthYear[r1],
                                           BirthYear2 = BirthYear[r2],
                                           stringsAsFactors = FALSE))
    }
    message("duplicate IDs found in lifehistory data, first entry will be used")
    LifeHistData <- LifeHistData[!duplicated(LifeHistData$ID), ]
  }
  if (returnDups & !all(is.na(gID))) {
    NoLH <- setdiff(gID, LifeHistData$ID)
    if (length(NoLH)>0)  OUT$NoLH <- NoLH
  }


  # check: compliance with hard coded max. age difference of 100 ----
  # (could be adjusted)
  MaxAgeDif <- with(LifeHistData, suppressWarnings(
    diff(range(BirthYear[BirthYear >= 0 & (ID %in% gID | is.na(gID))], na.rm = TRUE))))
  if (MaxAgeDif > 100) stop("Cannot handle >100 cohorts!", call.=FALSE)


  # return results ----
  if (sorted)
    LifeHistData <- orderLH(LifeHistData, gID)

  if (returnDups) {
    return( c(OUT, list(LifeHistData=LifeHistData)) )
  } else {
    return( LifeHistData )
  }
}


#=======================================================================
#=======================================================================
#' @title Order Lifehistory Data
#'
#' @description Order lifehistory data to match order of IDs in genotype data,
#'   filling in gaps with missing values.
#'
#' @param LH dataframe with lifehistory information:
#' \describe{
#'  \item{ID}{max. 30 characters long,}
#'  \item{Sex}{1 = females, 2 = males, other numbers = unknown,}
#'  \item{Birth Year}{(or hatching year) Use negative numbers to denote
#'  missing values.}
#'  \item{BY.min}{minimum birth year (optional)}
#'  \item{BY.max}{maximum birth year (optional)}}
#' @param gID character vector with IDs in genotype data, in order of
#'   occurrence.
#'
#' @return A dataframe with the same 5 columns, but with individuals in exactly
#'   the same order as gID, including padding with 'empty' rows if an individual
#'   in gID was not in the input-LH. Missing values are recoded to 3 for the
#'   'Sex' column, and -999 for the birth year columns.
#'
#' @keywords internal

orderLH <- function(LH=NULL, gID=NULL) {
  if (!all(gID %in% LH$ID)) {
    LHX <- data.frame(ID = gID[!gID %in% LH$ID],
                      Sex = 3,
                      BirthYear = -999,
                      BY.min = -999,
                      BY.max = -999)
  }
  if (is.null(LH)) {
    LHF <- LHX
  } else {
    if (!all(gID %in% LH$ID)) {
      LH <- merge(LHX, LH, all = TRUE)
    }
    LHF <- LH[match(gID, LH$ID), ]
    # fill in gaps: some gID may not be in LH
    LHF$Sex[is.na(LHF$Sex)] <- 3
    for (y in c("BirthYear", "BY.min", "BY.max")) {
      LHF[is.na(LHF[,y]), y] <- -999
    }
  }
  return( LHF )
}
