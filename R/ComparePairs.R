#============================================================================
#' @title Compare Pairwise Relationships
#'
#' @description Compare, count and identify different types of relative pairs
#'   between two pedigrees, or within one pedigree.
#'
#' @details If \code{Pairs2} is as returned by \code{\link{GetMaybeRel}}
#'   (identified by the additional column names 'LLR' and 'OH'), these
#'   relationship categories are appended with an '?' in the output, to
#'   distinguish them from those derived from \code{Ped2}.
#'
#'   When \code{Pairs2$TopRel} contains values other than the ones listed among
#'   the return values for the combination of \code{patmat} and \code{GenBack},
#'   they are prioritised in decreasing order of factor levels, or in decreasing
#'   alphabetical order, and before the default (\code{ped2} derived) levels.
#'
#'   The matrix returned by \code{\link{DyadCompare}} [Deprecated] is a subset
#'   of the matrix returned here using default settings.
#'
#' @param Ped1 first (e.g. original/reference) pedigree, dataframe with 3
#'   columns: id-dam-sire.
#' @param Ped2 optional second (e.g. inferred) pedigree.
#' @param Pairs2 optional dataframe with relationships categories between pairs
#'   of individuals, instead of or in addition to \code{Ped2}, e.g. as returned
#'   by \code{\link{GetMaybeRel}}. First three columns: ID1-ID2-relationship,
#'   column names and any additional columns are ignored.
#' @param GenBack  number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grandparental,
#'   avuncular and first cousins. GenBack >2 is not implemented.
#' @param patmat  logical, distinguish between paternal versus maternal relative
#'   pairs?
#' @param ExcludeDummies  logical, exclude dummy IDs from output? Individuals
#'   with e.g. the same dummy father will still be counted as paternal halfsibs.
#'   No attempt is made to match dummies in one pedigree to individuals in the
#'   other pedigree; for that use \code{\link{PedCompare}}.
#' @param  DumPrefix  character vector with the prefixes identifying dummy
#'   individuals. Use 'F0' ('M0') to avoid matching to regular individuals with
#'   IDs starting with 'F' ('M'), provided \code{Ped2} has fewer than 999 dummy
#'   females (males).
#' @param Return  return a matrix with \code{Counts} or a \code{Summary} of the
#'   number of identical relationships and mismatches per relationship, or
#'   detailed results as a 2xNxN \code{Array} or as a \code{Dataframe}.
#'   \code{All} returns a list with all four.
#'
#' @return Depending on \code{Return}, one of the following, or a list with all:
#' \item{Counts}{(the default), a matrix with counts, with the classification in
#'   \code{Ped1} on rows and that in \code{Ped2} in columns. Counts for
#'   'symmetrical' pairs ("FS", "HS", "MHS", "PHS", "FC1", "DFC1", "U","X") are
#'   divided by two.}
#'   \item{Summary}{a matrix with one row per relationship type and four columns
#'     , named as if \code{Ped1} is the true pedigree:
#'     \describe{
#'       \item{n}{total number of pairs with that relationship in \code{Ped1},
#'         and occurring in \code{Ped2}}
#'      \item{OK}{Number of pairs with same relationship in \code{Ped2} as in
#'        \code{Ped1}}
#'      \item{hi}{Number of pairs with 'higher' relationship in \code{Ped2} as
#'        in \code{Ped1} (e.g. FS instead of HS; ranking is the order given
#'        below)}
#'      \item{lo}{Number of pairs with 'lower' relationship in \code{Ped2} as in
#'       \code{Ped1}, but not unrelated in \code{Ped2}}
#'   }}
#'   \item{Array}{a 2xNxN array (if \code{Ped2} or \code{Pairs2} is specified)
#'     or a NxN matrix , where N is the total number of individuals occurring in
#'     \code{Ped1} and/or \code{Ped2}.}
#'   \item{Dataframe}{a dataframe with \eqn{N^2} rows and four columns:
#'     \describe{
#'       \item{id.A}{First individual of the pair}
#'       \item{id.B}{Second individual of the pair}
#'       \item{RC1}{the relationship category in \code{Ped1}, as a factor with all
#'         considered categories as levels, including those with 0 count}
#'       \item{RC2}{the relationship category in \code{Ped2}}
#'     }
#'     Each pair is listed twice, e.g. once as P and once as O, or twice as FS.}
#'
#' @section Relationship abbreviations and ranking:
#' By default (\code{GenBack=1, patmat=FALSE}) the following 7 relationships are
#' distinguished:
#' \itemize{
#'    \item \strong{S}: Self (not included in \code{Counts})
#'    \item \strong{MP}: Parent
#'    \item \strong{O}: Offspring (not included in \code{Counts})
#'    \item \strong{FS}: Full sibling
#'    \item \strong{HS}: Half sibling
#'    \item \strong{U}: Unrelated, or otherwise related
#'    \item \strong{X}: Either or both individuals not occurring in both
#'      pedigrees
#' }
#' In the array and dataframe, 'MP' indicates that the second (column)
#' individual is the parent of the first (row) individual, and 'O' indicates the
#' reverse.
#'
#' When \code{GenBack=1, patmat=TRUE} the categories are (S)-M-P-(O)-FS-MHS-PHS-
#' U-X.
#'
#' When \code{GenBack=2, patmat=TRUE}, the following relationships are
#' distinguished:
#' \itemize{
#'    \item \strong{S}: Self (not included in \code{Counts})
#'    \item \strong{M}: Mother
#'    \item \strong{P}: Father
#'    \item \strong{O}: Offspring (not included in \code{Counts})
#'    \item \strong{FS}: Full sibling
#'    \item \strong{MHS}: Maternal half-sibling
#'    \item \strong{PHS}: Paternal half-sibling
#'    \item \strong{MGM}: Maternal grandmother
#'    \item \strong{MGF}: Maternal grandfather
#'    \item \strong{PGM}: Paternal grandmother
#'    \item \strong{PGF}: Paternal grandfather
#'    \item \strong{GO}: Grand-offspring (not included in \code{Counts})
#'    \item \strong{FA}: Full avuncular; maternal or paternal aunt or uncle
#'    \item \strong{HA}: Half avuncular
#'    \item \strong{FN}: Full nephew/niece (not included in \code{Counts})
#'    \item \strong{HN}: Half nephew/niece (not included in \code{Counts})
#'    \item \strong{FC1}: Full first cousin
#'    \item \strong{DFC1}: Double full first cousin
#'    \item \strong{U}: Unrelated, or otherwise related
#'    \item \strong{X}: Either or both individuals not occurring in both pedigrees
#' }
#' Note that for avuncular and cousin relationships no distinction is made
#' between paternal versus maternal, as this may differ between the two
#' individuals and would generate a large number of sub-classes. When a pair is
#' related via multiple paths, the first-listed relationship is returned.
#'
#' When \code{GenBack=2, patmat=FALSE}, MGM, MGF, PGM and PGF are combined
#' into GP, with the rest of the categories analogous to the above.
#'
#' @seealso \code{\link{PedCompare}} for individual-based comparison;
#'   \code{\link{GetRelM}} for a pairwise relationships matrix of a single
#'   pedigree; \code{\link{PlotRelPairs}} for visualisation of relationships
#'   within each pedigree.
#'
#'   To estimate P(actual relationship (Ped1) | inferred relationship (Ped2)),
#'   see examples at \code{\link{EstConf}}.
#'
#' @examples
#' \donttest{
#' data(Ped_griffin, SeqOUT_griffin, package="sequoia")
#' PairsG <- ComparePairs(Ped_griffin, SeqOUT_griffin[["Pedigree"]],
#'                        patmat = TRUE, ExcludeDummies = TRUE, Return = "All")
#' PairsG$Counts
#'
#' # pairwise correct assignment rate:
#' PairsG$Summary[,"OK"] / PairsG$Summary[,"n"]
#'
#' # check specific pair:
#' PairsG$Array[, "i190_2010_M", "i168_2009_F"]
#' # or
#' RelDF <- PairsG$Dataframe   # for brevity
#' RelDF[RelDF$id.A=="i190_2010_M" & RelDF$id.B=="i168_2009_F", ]
#'
#' # Colony-style lists of full sib dyads & half sib dyads:
#' FullSibDyads <- with(RelDF, RelDF[Ped1 == "FS" & id.A < id.B, ])
#' HalfSibDyads <- with(RelDF, RelDF[Ped1 == "HS" & id.A < id.B, ])
#' # Use 'id.A < id.B' because each pair is listed 2x
#' }
#'
#' @export

ComparePairs <- function(Ped1 = NULL,
                         Ped2 = NULL,
                         Pairs2 = NULL,
                         GenBack = 1,
                         patmat = FALSE,
                         ExcludeDummies = TRUE,
                         DumPrefix = c("F0", "M0"),
                         Return = "Counts")
{
  if(is.null(Ped1)) stop("No 'Ped1' provided")
  Return <- .simpleCap(Return)  # capitalise 1st letter
  if (!Return %in% c("Array", "Counts", "Dataframe", "Summary", "All"))
    stop("'Return' must be 'Counts', 'Array', 'Dataframe', 'Summary' or 'All'")
  if (Return == "Summary" & is.null(Ped2) & is.null(Pairs2))
    stop("Return='Summary' not availble for a single pedigree")
  # if no Ped2 & no Pairs2, 'Summary' is identical to 'Counts'.

  if (!GenBack %in% 1:2)  stop("'GenBack' must be 1 or 2")
  if (!patmat %in% c(TRUE,FALSE))  stop("'patmat' must be TRUE or FALSE")

  Ped1 <- PedPolish(Ped1, ZeroToNA=TRUE, NullOK=FALSE, StopIfInvalid=FALSE)
  if (GenBack==2)  Ped1 <- GPcols(Ped1)
  RCM.1 <- GetRelM(Ped1, GenBack=GenBack, patmat=patmat, Return="Matrix")

  lvls <- list(GB1 = list(no = c("S", "MP", "O", "FS", "HS", "U", "X"),
                          yes = c("S","M", "P", "O", "FS", "MHS", "PHS", "U", "X")),
               GB2 = list(no = c("S","MP", "O", "FS", "HS", "GP", "GO", "FA",
                                 "HA", "FN", "HN", "FC1", "DFC1", "U", "X"),
                          yes = c("S","M", "P", "O", "FS", "MHS", "PHS","MGM",
                                  "MGF", "PGM", "PGF", "GO",
                                  "FA", "FN", "HA", "HN", "DFC1", "FC1","U", "X")))

  if (!is.null(Ped2)) {
    Ped2 <- PedPolish(Ped2, ZeroToNA=TRUE, StopIfInvalid=FALSE)[, 1:3]
    if (!any(Ped2$id %in% Ped1$id))  stop("no common IDs in Ped1 and Ped2")
    if (GenBack==2)  Ped2 <- GPcols(Ped2)
    RCM.2 <- GetRelM(Ped2, GenBack=GenBack, patmat=patmat, Return="Matrix")
  }

  RelRank <- c("S", "M", "P", "MP", "O", "PO?",
               "FS","FS?", "MHS", "PHS", "HS", "HS?",
               "MGM", "MGF", "PGM", "PGF", "GP", "GO","GP?",
               "FA", "FN", "FA?", "2nd?", "HA", "HN","HA?",
               "DFC1", "FC1", "XX?", "Q?", "U", "X")
	#@@@@@@@@@@@@@@@@@@@

	if (!is.null(Pairs2)) {
	  if (!class(Pairs2) %in% c("data.frame", "matrix"))  stop("Pairs2 should be a dataframe or matrix")
	  MaybeRelNames <- c("PO", "FS", "HS", "GP", "FA", "2nd", "HA", "Q")
	  if (all(c("ID1", "ID2", "TopRel", "LLR", "OH") %in% names(Pairs2)) &&
	      all(Pairs2$TopRel %in% MaybeRelNames)) {
	    Pairs2$TopRel <- paste0(Pairs2$TopRel, "?")  # output from GetMaybeRel()
	  }
		if (!"TopRel" %in% names(Pairs2))  names(Pairs2)[1:3] <- c("ID1", "ID2", "TopRel")

		MaybeRelNames <- c(paste0(MaybeRelNames, "?"), "U")

		if (is.factor(Pairs2$TopRel)) {
		  lvls.X2 <- levels(Pairs2$TopRel)
		} else {
  		lvls.X2 <- sort(unique(na.exclude(Pairs2$TopRel)))
		}
	  if (all(lvls.X2 %in% MaybeRelNames)) {  # !is.null(Ped2) &
		  lvls.X2 <- MaybeRelNames
	  }
		lvls.Y2 <- lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]

		if (!is.null(Ped2) | all(lvls.X2 %in% lvls.Y2)) {
		  lvls2 <- unique(c(lvls.X2, lvls.Y2))
		} else {
		  lvls2 <- lvls.X2
		}
		lvls2 <- c(RelRank[RelRank %in% lvls2], lvls2[!lvls2 %in% RelRank])


    RCM.X2 <- GetRelM(Pairs = Pairs2[, c("ID1", "ID2", "TopRel")],
                      GenBack=GenBack, patmat=patmat)

		if (is.null(Ped2)) {
			RCM.2 <- RCM.X2
		} else {
#      RCM.X2 <- matrix(paste0(RCM.X2, "?"), dim(RCM.X2))  # enable distinction between Ped2 & Pairs2 origin
			IDs2 <- unique(c(colnames(RCM.2), colnames(RCM.X2)))
			RCM.2i <- inflate(RCM.2, IDs2)
			RCM.X2i <- inflate(RCM.X2, IDs2)
			RCM.2i[is.na(RCM.2i)] <- "X"
			RCM.X2i[is.na(RCM.X2i)] <- "X"
			RCM.P2.f <- factor(RCM.2i, levels=lvls2)
			RCM.X2.f <- factor(RCM.X2i, levels=lvls2)
			RCM.2 <- matrix(factor(pmin(as.numeric(RCM.P2.f), as.numeric(RCM.X2.f), na.rm=TRUE),
			                          levels = seq_along(lvls2), labels=lvls2),
			                   length(IDs2), length(IDs2), dimnames=list(IDs2, IDs2))
		}
  } else {
    if (!is.null(Ped2)) {
		  lvls2 <- lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]
    } else {
      lvls2 <- "X"
      RCM.2 <- NULL
    }
		lvls.X2 <- NULL
	}

	# align the two matrices into an array
  IDs <- unique(c(colnames(RCM.1),  colnames(RCM.2)))
  RCA <- array(dim=c(2, length(IDs), length(IDs)),
               dimnames=list(c("Ped1", "Ped2"), IDs, IDs))
  RCA["Ped1", colnames(RCM.1), colnames(RCM.1)] <- RCM.1
  RCA["Ped2", colnames(RCM.2), colnames(RCM.2)] <- RCM.2
  RCA[is.na(RCA)] <- "X"

	# delete dummies (doing this earlier may cause trouble with 2-generations-back GetRelM)
	if (ExcludeDummies) {
	  DPnc <- nchar(DumPrefix)
    Dummies <- rep(FALSE, length(IDs))
    for (x in seq_along(DumPrefix)) {
      Dummies <- ifelse(Dummies, Dummies, substr(IDs,1,DPnc[x]) == DumPrefix[x])
    }
	  if (sum(Dummies)>0) {
	    RCA <- RCA[, !Dummies, !Dummies]
	    IDs <- IDs[!Dummies]
	  }
	}


  # Plot ----
  # if (Plot) {
    # PlotRelPairs(RCA["Ped1", , ])
    # if (!is.null(Ped2) | !is.null(Pairs2)) {
      # if (interactive()) {
        # inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
      # }
      # PlotRelPairs(RCA["Ped2", , ])
    # }
  # }


  # Array ----
	if (Return == "Array") {
	  if (!is.null(Ped2) | !is.null(Pairs2)) {
		  return( RCA )
	  } else {
	    return( RCA["Ped1", , , drop=TRUE])
	  }
	}

  # Counts ----
	if (Return %in% c("Counts", "Summary", "All")) {
	  lvls.tbl <- list(GB1 = list(no = c("MP", "FS", "HS", "U", "X"),
	                              yes = c("M", "P", "FS", "MHS", "PHS", "U", "X")),
	                   GB2 = list(no = c("MP", "FS", "HS", "GP", "FA", "HA", "FC1", "DFC1", "U", "X"),
	                              yes = c("M", "P", "FS", "MHS", "PHS","MGM", "MGF", "PGM", "PGF",
	                                      "FA", "HA","FC1", "DFC1", "U", "X")))
	  tbl <- table(Ped1 = factor(RCA["Ped1",,], levels=lvls.tbl[[GenBack]][[ifelse(patmat, "yes", "no")]]),
	               Ped2 = factor(RCA["Ped2",,], levels=lvls2[lvls2 != "S"]))
	  dup <- c("FS", "HS", "MHS", "PHS", "FC1", "DFC1", "U","X",   # pairs included double
	           lvls.X2)
	  these.dup <- rownames(tbl) %in% dup
		those.dup <- colnames(tbl) %in% dup
	  tbl[these.dup, those.dup] <- tbl[these.dup, those.dup]/2
	  if (Return == "Counts") return( tbl )

	}

  # Dataframe ----
	if (Return %in% c("Dataframe", "All")) {
	  DF <- data.frame(id.A = rep(IDs, times=length(IDs)),
	                   id.B = rep(IDs, each=length(IDs)),
	                   Ped1 = factor(RCA["Ped1",,],
	                                 levels=lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]),
	                   Ped2 = factor(RCA["Ped2",,], levels=lvls2))
	  if (Return == "Dataframe") return( DF )
	}

  # Summary ----
	if (Return %in% c("Summary", "All") & !is.null(RCM.2)) {
	  # if no Ped2 & no Pairs2, 'Summary' is identical to 'Counts'.
	  tblz <- tbl[rownames(tbl) != "X", colnames(tbl) != "X", drop=FALSE]
    rowR <- as.numeric(factor(rownames(tblz), levels=RelRank))
    colR <- as.numeric(factor(colnames(tblz), levels=RelRank))
    hi <- tblz * outer(rowR, colR, function(x,y) x > y)
    eq <- tblz * outer(rowR, colR, function(x,y) x == y)
    tblzz <- tblz
    if ("U" %in% colnames(tblz))  tblzz[,"U"] <- 0   # not if Pairs2
    lo <- tblzz * outer(rowR, colR, function(x,y) x < y)
    ARER <- cbind(n = rowSums(tblz),
                  OK = rowSums(eq),
                  lo = rowSums(lo),
                  hi = rowSums(hi))
    if (Return == "Summary") return( ARER )
	} else if (Return %in% c("Summary", "All") & is.null(RCM.2)) {
	  ARER <- NULL
	}

	if (Return == "All") {
	  return(list(Array = RCA,
	              Counts = tbl,
	              Dataframe = DF,
	              Summary = ARER))
	}
}


#============================================================================
#============================================================================
#' @title Compare Dyads
#'
#' @description Count the number of half and full sibling pairs correctly and
#'   incorrectly assigned. DEPRECATED - PLEASE USE \code{\link{ComparePairs}}
#'
#' @param  Ped1 original pedigree, dataframe with 3 columns: id-dam-sire.
#' @param  Ped2 second (inferred) pedigree.
#' @param  na1  the value for missing parents in Ped1.
#'
#' @return A 3x3 table with the number of pairs assigned as full siblings (FS),
#'   half siblings (HS) or unrelated (U, including otherwise related) in the two
#'   pedigrees, with the classification in Ped1 on rows and that in Ped2 in
#'   columns.
#'
#' @seealso \code{\link{ComparePairs}} which supersedes this function;
#'   \code{\link{PedCompare}}
#'
#' @examples
#' \donttest{
#' data(Ped_HSg5, SimGeno_example, LH_HSg5, package="sequoia")
#' SeqOUT <- sequoia(GenoM = SimGeno_example, LifeHistData = LH_HSg5,
#'                   Module="par", quiet=TRUE, Plot=FALSE)
#' DyadCompare(Ped1=Ped_HSg5, Ped2=SeqOUT$Pedigree)
#' }
#'
#' @export

DyadCompare <- function(Ped1 = NULL,
                        Ped2 = NULL,
                        na1 = c(NA, "0"))
{
  warning("This function is deprecated, please use ComparePairs()",
          immediate.=TRUE)
  if(is.null(Ped1) || nrow(Ped1)<2) stop("No 'Ped1' provided")
  if(is.null(Ped2) || nrow(Ped2)<2) stop("No 'Ped2' provided'")
  names(Ped1)[1:3] <- c("id", "dam.1", "sire.1")
  names(Ped2)[1:3] <- c("id", "dam.2", "sire.2")
  for (i in 1:3) {
    Ped1[, i] <- as.character(Ped1[, i])
    Ped1[Ped1[, i] %in% na1, i] <- NA
  }
  for (i in 1:3) Ped2[, i] <- as.character(Ped2[, i])
  if (!any(Ped2$id %in% Ped1$id))  stop("no common IDs in Ped1 and Ped2")
  Ped1 <- Ped1[!is.na(Ped1$id), 1:3]
  Ped2 <- Ped2[!is.na(Ped2$id), 1:3]
  Ped1 <- PedPolish(Ped1)
  Ped2 <- PedPolish(Ped2)

  # note: each pair is counted double
  RCT <- matrix(NA, 0, 3)
  for (x in 1:nrow(Ped1)) {
    RCT <- rbind(RCT, rc(x, Ped1))
  }

  RCI <- matrix(NA, 0, 3)
  for (x in 1:nrow(Ped2)) {
    RCI <- rbind(RCI, rc(x, Ped2))
  }

  RCTI <- merge(as.data.frame(RCT, stringsAsFactors=FALSE),
                as.data.frame(RCI, stringsAsFactors=FALSE),
                by=c("id1", "id2"), all=TRUE, suffixes = c(".1", ".2"))
  RCTI <- RCTI[RCTI$id1 %in% Ped1$id & RCTI$id2 %in% Ped1$id &
                 RCTI$id1 %in% Ped2$id & RCTI$id2 %in% Ped2$id, ]
  RCTI$RC.1[is.na(RCTI$RC.1)] <- "U"
  RCTI$RC.2[is.na(RCTI$RC.2)] <- "U"
  RCTI$RC.1 <- factor(RCTI$RC.1, levels=c("FS", "HS", "U"))
  RCTI$RC.2 <- factor(RCTI$RC.2, levels=c("FS", "HS", "U"))

  tbl <- with(RCTI, table(RC.1, RC.2, useNA="ifany"))/2  # pairs included double
  tbl["U", "U"] <- nrow(Ped2) * (nrow(Ped2)-1)/2 - sum(tbl)
  tbl
  #  sweep(tbl, 1, rowSums(tbl), "/")
}

#============================================================================
#============================================================================

#' @title Find siblings

#' @param x  an ID
#' @param Ped  a pedigree with columns id - dam - sire
#'
#' @return The individuals which are full or half siblings to x, as a
#'   three-column matrix with column names id1 (x), id2 (the siblings), and
#'   RC (the relatedness category, 'FS' or 'HS').
#'
#' @keywords internal

rc <- function(x, Ped) {
  names(Ped) <- c("id", "dam", "sire")
  RelCat <- with(Ped,
                 ifelse(id == id[x], "S",
                        ifelse(eqv(dam[x],dam,FALSE) & eqv(sire[x], sire,FALSE), "FS",
                               ifelse(eqv(dam[x],dam,FALSE) |  eqv(sire[x], sire,FALSE), "HS",
                                      NA))))
  out <- cbind(id1 = Ped$id[x],
               id2 = Ped$id[!is.na(RelCat)],
               RC = stats::na.exclude(RelCat))
  out <- out[out[,"RC"] != "S", ]
  out
}

#============================================================================
#============================================================================
