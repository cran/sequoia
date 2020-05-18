#============================================================================
#' @title Comparison of all pairwise relationships in 2 pedigrees
#'
#' @description Compare, count and identify different types of relative pairs
#'   between two pedigrees. The matrix returned by \code{\link{DyadCompare}}
#'   [Deprecated] is a subset of the matrix returned here using default
#'   settings.
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
#' @param Ped1 (Original/reference) pedigree, dataframe with 3 columns:
#'   id-dam-sire
#' @param Ped2 Second (inferred) pedigree
#' @param Pairs2 dataframe with relationships categories between pairs of
#'   individuals, instead of or in addition to Ped2, e.g. as returned by
#'   \code{\link{GetMaybeRel}}. First three columns: ID1-ID2-relationship,
#'   column names and any additional columns are ignored.
#' @param GenBack  Number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grandparental,
#'   avuncular and first cousins. GenBack >2 is not implemented.
#' @param patmat  logical, distinguish between paternal versus maternal relative
#'   pairs?
#' @param DumPrefix character vector of length 2 with the dummy prefixes in Ped1
#'   and/or Ped2. IDs starting with these prefixes will not be excluded, but
#'   individuals with dummy parents are compared. Use \code{\link{GetRelCat}} on
#'   a single pedigree to find relationships with dummies.
#' @param Return  Return a matrix with \code{Counts} or a \code{Summary} of the
#'   number of identical relationships and mismatches per relationship, or
#'   detailed results as a 2xNxN \code{Array} or as a \code{Dataframe}.
#'   \code{All} returns a list with all four.
#'
#' @return a matrix with counts, a 3D array or a 4-column dataframe, depending
#'   on \code{Return}, with by default (\code{GenBack=1, patmat=FALSE}) the
#'   following 7 relationships:
#'    \item{S}{Self (not in counts)}
#'    \item{MP}{Parent}
#'    \item{O}{Offspring (not in counts)}
#'    \item{FS}{Full sibling}
#'    \item{HS}{Half sibling}
#'    \item{U}{Unrelated, or otherwise related}
#'    \item{X}{Either or both individuals not occurring in both pedigrees}
#' Where in the array and dataframe, 'MP' indicates that the second (column)
#' individual is the parent of the first (row) individual, and 'O' indicates the
#' reverse.
#'
#' When \code{GenBack=2, patmat=TRUE}, the following relationships are
#' distinguished:
#'    \item{S}{Self (not in counts)}
#'    \item{M}{Mother}
#'    \item{P}{Father}
#'    \item{O}{Offspring (not in counts)}
#'    \item{FS}{Full sibling}
#'    \item{MHS}{Maternal half-sibling}
#'    \item{PHS}{Paternal half-sibling}
#'    \item{MGM}{Maternal grandmother}
#'    \item{MGF}{Maternal grandfather}
#'    \item{PGM}{Paternal grandmother}
#'    \item{PGF}{Paternal grandfather}
#'    \item{GO}{Grand-offspring (not in counts}
#'    \item{FA}{Full avuncular; maternal or paternal aunt or uncle}
#'    \item{HA}{Half avuncular}
#'    \item{FN}{Full nephew/niece (not in counts}
#'    \item{HN}{Half nephew/niece (not in counts}
#'    \item{FC1}{Full first cousin}
#'    \item{DFC1}{Double full first cousin}
#'    \item{U}{Unrelated, or otherwise related}
#'    \item{X}{Either or both individuals not occurring in both pedigrees}
#' Note that for avuncular and cousin relationships no distinction is made
#' between paternal versus maternal, as this may differ between the two
#' individuals and would generate a large number of subclasses. When a pair is
#' related via multiple paths, the first-listed relationship is returned.
#'
#' When \code{GenBack=1, patmat=TRUE} the categories are (S)-M-P-(O)-FS-MHS-PHS-
#' U-X. When \code{GenBack=2, patmat=FALSE}, MGM, MGF, PGM and PGF are combined
#' into GP, with the rest of the categories analogous to the above.
#'
#' Note that in the dataframe each pair is listed twice, e.g. once as P and once
#' as O, or twice as FS.
#'
#' When \code{Return = "Counts"} (the default), a matrix with counts is
#' returned, with the classification in Ped1 on rows and that in Ped2 in
#' columns. Counts for 'symmetrical' pairs ("FS", "HS", "MHS", "PHS", "FC1",
#' "DFC1", "U","X") are divided by two.
#'
#' When \code{Return = 'Summary'}, the counts table is distilled down into a matrix
#' with four columns, which names assuming \code{Ped1} is the true pedigree:
#'  \item{n}{total number of pairs with that relationship in Ped1}
#'  \item{OK}{Number of pairs with same relationship in Ped2 as in Ped1}
#'  \item{lo}{Number of pairs with 'lower' relationship in Ped2 as in Ped1 (see
#'  ranking above), but not unrelated in Ped2}
#'  \item{hi}{Number of pairs with 'higher' relationship in Ped2 as in Ped1}
#'
#' When \code{Return = "Array"}, the first dimension is 1=Ped1, 2=Ped2,
#' the 2nd and 3rd dimension are the two individuals of the pair.
#'
#' When \code{Return = "Dataframe"}, the columns are
#'   \item{id.A}{First individual of the pair}
#'   \item{id.B}{Second individual of the pair}
#'   \item{RC1}{the relationship category in Ped1, as a factor with all
#'     considered categories as levels, including those with 0 count}
#'   \item{RC2}{the relationship category in Ped2}
#'
#'
#' @seealso \code{\link{PedCompare}} for individual-based comparison;
#'   \code{\link{GetRelCat}} for pairs of relatives within a single pedigree.
#'
#'
#' @examples
#' \dontrun{
#' data(Ped_HSg5, SimGeno_example, LH_HSg5, package="sequoia")
#' SeqOUT <- sequoia(GenoM = SimGeno_example, LifeHistData = LH_HSg5,
#'                   MaxSibIter = 0)
#' ComparePairs(Ped1=Ped_HSg5, Ped2=SeqOUT$Pedigree, Return="Counts")
#' # matrix with counts of pairs
#' RC.A <- ComparePairs(Ped1=Ped_HSg5, Ped2=SeqOUT$Pedigree, Return="Array")
#' RC.A[, "a05017", "b05018"] # check specific pairs
#'
#' RC.DF <- ComparePairs(Ped1=Ped_HSg5, Ped2=SeqOUT$Pedigree,
#'   Return="Dataframe")
#' RC.DF[RC.DF$id.A=="a05017" & RC.DF$id.B=="b05018", ] # check specific pairs
#' table(RC.DF$Ped1, RC.DF$Ped2)
#' # incl. S,O,GO,FN,HN; duplicated counts for FS,HS,FC1,DFC1,U,X
#' Mismatches <- RC.DF[RC.DF$Ped1 != RC.DF$Ped2, ]
#'
#' Maybe <- GetMaybeRel(SimGeno_example, SeqList=SeqOUT, ParSib="sib")
#' cp <- ComparePairs(Ped1=Ped_HSg5, Ped2=SeqOUT$Pedigree,
#'                    Pairs2=Maybe$MaybeRel, Return="All")
#' cp$Counts[, colSums(cp$Counts)>0]
#' cp$Summary[,"OK"] / cp$Summary[,"n"]  # pairwise assignment rate
#'
#' }
#'
#' @export

ComparePairs <- function(Ped1 = NULL,
                         Ped2 = NULL,
                         Pairs2 = NULL,
                         GenBack = 1,
                         patmat = FALSE,
                         DumPrefix = c("F0", "M0"),
                         Return = "Counts")
{
  if(is.null(Ped1)) stop("No 'Ped1' provided")
  if(is.null(Ped2) & is.null(Pairs2)) stop("Must provide 'Ped2' or 'Pairs2'")
  Return <- .simpleCap(Return)  # capitalise 1st letter
	if (!Return %in% c("Array", "Counts", "Dataframe", "Summary", "All")) {
	    stop("'Return' must be 'Counts', 'Array', 'Dataframe', 'Summary' or 'All'")
	}
  if (!GenBack %in% 1:2)  stop("'GenBack' must be 1 or 2")
  if (!patmat %in% c(TRUE,FALSE))  stop("'patmat' must be TRUE or FALSE")

  Ped1 <- PedPolish(Ped1[,1:3], ZeroToNA=TRUE)
  RCM.1 <- t(sapply(1:nrow(Ped1), GetRelCat, Ped1, GenBack=GenBack, patmat=patmat))

	lvls <- list(GB1 = list(no = c("S", "MP", "O", "FS", "HS", "U", "X"),
	                              yes = c("S","M", "P", "O", "FS", "MHS", "PHS", "U", "X")),
	                   GB2 = list(no = c("S","MP", "O", "FS", "HS", "GP", "GO", "FA", "HA", "FC1", "DFC1", "U", "X"),
	                              yes = c("S","M", "P", "O", "FS", "MHS", "PHS","MGM", "MGF", "PGM", "PGF", "GO",
	                                      "FA", "FN", "HA", "HN", "DFC1", "FC1","U", "X")))

  if (!is.null(Ped2)) {
    Ped2 <- PedPolish(Ped2[,1:3], ZeroToNA=TRUE)
    if (!any(Ped2$id %in% Ped1$id))  stop("no common IDs in Ped1 and Ped2")
	  RCM.2 <- t(sapply(1:nrow(Ped2), GetRelCat, Ped2, GenBack=GenBack, patmat=patmat))
#	  rownames(RCM.2) <- colnames(RCM.2)
  }

	#@@@@@@@@@@@@@@@@@@@
	inflate <- function(M, IDnew) {
		Mnew <- matrix(NA, length(IDnew), length(IDnew), dimnames=list(IDnew, IDnew))
		if (is.null(rownames(M)) & nrow(M)==ncol(M))  rownames(M) <- colnames(M)
		Mnew[rownames(M), colnames(M)] <- M
		return( Mnew )
	}
	#@@@@@@@@@@@@@@@@@@@

	if (!is.null(Pairs2)) {
	  MaybeRelNames <- c("PO", "FS", "HS", "GP", "FA", "2nd", "HA", "Q")
	  if (all(c("ID1", "ID2", "TopRel", "LLR", "OH") %in% names(Pairs2)) &&
	      all(Pairs2$TopRel %in% MaybeRelNames)) {
	    Pairs2$TopRel <- paste0(Pairs2$TopRel, "?")  # output from GetMaybeRel()
	  }
		if (!"TopRel" %in% names(Pairs2))  names(Pairs2)[1:3] <- c("ID1", "ID2", "TopRel")

		MaybeRelNames <- c(paste0(MaybeRelNames, "?"), "U")
		RelRank <- c("S", "M", "P", "MP", "PO?", "O",
		             "FS","FS?", "MHS", "PHS", "HS", "HS?",
		             "MGM", "MGF", "PGM", "PGF", "GP", "GO","GP?",
		             "FA", "FN", "FA?", "2nd?", "HA", "HN","HA?",
		             "DFC1", "FC1", "XX?", "U", "X")

		if (is.factor(Pairs2$TopRel)) {
		  lvls.X2 <- levels(Pairs2$TopRel)
		} else {
  		lvls.X2 <- sort(unique(na.exclude(Pairs2$TopRel)))
		}
		if (all(lvls.X2 %in% lvls[[GenBack]][[ifelse(patmat, "yes", "no")]])) {
			lvls2 <- lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]
		} else if (!is.null(Ped2) & all(lvls.X2 %in% MaybeRelNames)) {
		  lvls.X2 <- MaybeRelNames
	    lvls2 <- RelRank
	  } else if (!is.null(Ped2)) {
	    lvls2 <- unique(c(lvls.X2, lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]))
	  } else {
      lvls2 <- lvls.X2
    }

		Pairs2$TopRel <- as.character(Pairs2$TopRel)
		RCM.X2.tmp <- plyr::daply(Pairs2,
             .variables=c("ID1", "ID2"), .fun=function(df) df$TopRel)  # assuming Pairs2 = Maybe$MaybePar
		# make into symmetrical matrix to get consistency in above/below diagonal:
		IDsX2 <- unique(c(as.character(Pairs2$ID1), as.character(Pairs2$ID2)))
		RCM.X2a <- inflate(RCM.X2.tmp, IDsX2)
		RCM.X2b <- inflate(t(RCM.X2.tmp), IDsX2)
		if(any(RCM.X2a != RCM.X2b, na.rm=TRUE)) {
			stop("One or more pairs occur 2x in 'Pairs2', with different relationship")
		}
		RCM.X2 <- RCM.X2a
		RCM.X2[,] <- "U"
		RCM.X2[!is.na(RCM.X2a)] <- RCM.X2a[!is.na(RCM.X2a)]
		RCM.X2[!is.na(RCM.X2b)] <- RCM.X2b[!is.na(RCM.X2b)]

		if (is.null(Ped2)) {
			RCM.2 <- RCM.X2
		} else {
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
		lvls2 <- lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]
		lvls.X2 <- NULL
	}

	# align the two matrices into an array
	IDs <- unique(c(colnames(RCM.1),  colnames(RCM.2)))
	RCA <- array(dim=c(2, length(IDs), length(IDs)),
							 dimnames=list(c("Ped1", "Ped2"), IDs, IDs))
	RCA["Ped1", colnames(RCM.1), colnames(RCM.1)] <- RCM.1
	RCA["Ped2", colnames(RCM.2), colnames(RCM.2)] <- RCM.2
	RCA[is.na(RCA)] <- "X"

	# delete dummies (doing this earlier may cause trouble with 2-generations-back GetRelCat)
	DPnc <- nchar(DumPrefix)
	Dummies <- substr(IDs,1,DPnc[1])==DumPrefix[1] | substr(IDs,1,DPnc[2])==DumPrefix[2]
	if (sum(Dummies)>0) {
	  RCA <- RCA[, !Dummies, !Dummies]
	  IDs <- IDs[!Dummies]
	}

	#~~~~~~~~~~~
	if (Return == "Array") {
		return( RCA )

	}
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
	if (Return %in% c("dataframe", "All")) {
	  DF <- data.frame(id.A = rep(IDs, times=length(IDs)),
	                   id.B = rep(IDs, each=length(IDs)),
	                   Ped1 = factor(RCA["Ped1",,],
	                                 levels=lvls[[GenBack]][[ifelse(patmat, "yes", "no")]]),
	                   Ped2 = factor(RCA["Ped2",,], levels=lvls2))
	  if (Return == "Dataframe") return( DF )
	}
	if (Return %in% c("Summary", "All")) {
	  tblz <- tbl[rownames(tbl) != "X", colnames(tbl) != "X"]
    rowR <- as.numeric(factor(rownames(tblz), levels=RelRank))
    colR <- as.numeric(factor(colnames(tblz), levels=RelRank))
    hi <- tblz * outer(rowR, colR, function(x,y) x > y)
    eq <- tblz * outer(rowR, colR, function(x,y) x == y)
    tblzz <- tblz
    tblzz[,"U"] <- 0
    lo <- tblzz * outer(rowR, colR, function(x,y) x < y)
    ARER <- cbind(n = rowSums(tblz),
                  OK = rowSums(eq),
                  lo = rowSums(lo),
                  hi = rowSums(hi))
    if (Return == "Summary") return( ARER )
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
#' @title Compare dyads
#'
#' @description Count the number of half and full sibling pairs correctly and
#'   incorrectly assigned. DEPRECATED - SEE \code{\link{ComparePairs}}
#'
#' @param  Ped1 Original pedigree, dataframe with 3 columns: id-dam-sire
#' @param  Ped2 Second (inferred) pedigree
#' @param  na1  the value for missing parents in Ped1.
#'
#' @return A 3x3 table with the number of pairs assigned as full siblings (FS),
#'   half siblings (HS) or unrelated (U, including otherwise related) in the two
#'   pedigrees, with the classification in Ped1 on rows and that in Ped2 in
#'   columns
#'
#' @seealso \code{\link{PedCompare}}
#'
#' @examples
#' \dontrun{
#' data(Ped_HSg5, SimGeno_example, LH_HSg5, package="sequoia")
#' SeqOUT <- sequoia(GenoM = SimGeno_example, LifeHistData = LH_HSg5,
#'                   MaxSibIter = 0)
#' DyadCompare(Ped1=Ped_HSg5, Ped2=SeqOUT$Pedigree)
#' }
#'
#' @export

DyadCompare <- function(Ped1 = NULL,
                        Ped2 = NULL,
                        na1 = c(NA, "0"))
{
  message("This function is deprecated, please use ComparePairs()")
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

  tbl <- with(RCTI, Table(RC.1, RC.2))/2  # pairs included double
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
