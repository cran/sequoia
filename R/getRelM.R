#' @title Matrix with Pairwise Relationships
#'
#' @description Generate a matrix or 3D array with all pairwise relationships
#'   from a pedigree or dataframe with pairs.
#'
#' @param Pedigree  dataframe with columns id - dam - sire.
#' @param Pairs  dataframe with columns ID1 - ID2 - Rel.
#' @param GenBack  number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grand-parental,
#'   avuncular and first cousins.
#' @param patmat  logical, distinguish between paternal versus maternal relative
#'   pairs? For avuncular pairs, the distinction is never made.
#' @param Return  'Matrix' or 'Array'. The former returns an N x N matrix with
#'   the closest relationship between each pair, the latter an N x N x R array
#'   with for each of the R considered relationships whether it exists between
#'   the pair (1) or not (0). See Details below.
#'
#' @return  If \code{Return='Matrix'}, an N x N square matrix, with N equal to
#'   the number of rows in \code{Pedigree} (after running
#'   \code{\link{PedPolish}}) or the number of unique individuals in
#'   \code{Pairs}. If \code{Return='Array'}, an N x N x R array is returned,
#'   with R, the number of different relationships, determined by \code{GenBack}
#'   and \code{patmat}.
#'
#'  The following abbreviations are used within the returned \code{Matrix}, or as
#'  names of the 3rd dimension in the \code{Array}:
#'    \item{S}{Self}
#'    \item{M}{Mother}
#'    \item{P}{Father}
#'    \item{MP}{Mother or Father (\code{patmat=FALSE})}
#'    \item{O}{Offspring}
#'    \item{FS}{Full sibling}
#'    \item{MHS}{Maternal half-sibling}
#'    \item{PHS}{Paternal half-sibling}
#'    \item{XHS}{other half-sibling (hermaphrodites)}
#'    \item{HS}{half-sibling (\code{patmat=FALSE})}
#'    \item{MGM}{Maternal grandmother}
#'    \item{MGF}{Maternal grandfather}
#'    \item{PGM}{Paternal grandmother}
#'    \item{PGF}{Paternal grandfather}
#'    \item{GP}{Grandparent (\code{patmat=FALSE})}
#'    \item{GO}{Grand-offspring}
#'    \item{FA}{Full avuncular; maternal or paternal aunt or uncle}
#'    \item{HA}{Half avuncular}
#'    \item{FN}{Full nephew/niece}
#'    \item{HN}{Half nephew/niece}
#'    \item{FC1}{Full first cousin}
#'    \item{DFC1}{Double full first cousin}
#'    \item{U}{Unrelated (or otherwise related)}
#'
#' @details  Double relationships are ignored when \code{Return='Matrix'}, but
#'   not when \code{Return='Array'}. For example, when A and B are both
#'   mother-offspring and paternal siblings (A mated with her father to produce
#'   B), only the mother-offspring relationship will be indicated when
#'   \code{Return='Matrix'}.
#'
#'   Note that full siblings are the exception to this rule: in the \code{Array}
#'   they will be indicated as 'FS' only, and not as 'MHS' or 'PHS'. Similarly,
#'   full avuncular pairs are not indicated as 'HA'.
#'
#'   When \code{Pairs} is provided, \code{GenBack} and \code{patmat} are
#'   ignored, and no check is performed if the abbreviations are compatible with
#'   other functions.
#'
#' @seealso \code{\link{ComparePairs}} for comparing pairwise relationships
#'   between two pedigrees; \code{\link{PlotRelPairs}}.
#'
#' @examples
#' data(Ped_griffin)
#' Rel.griffin <- GetRelM(Ped_griffin, patmat=TRUE, GenBack=2)
#' table(c(Rel.griffin))
#' # turning matrix into vector first makes table() much faster
#' PlotRelPairs(Rel.griffin)
#'
#' @export

GetRelM <- function(Pedigree = NULL,
                    Pairs = NULL,
                    GenBack = 1,
                    patmat = FALSE,
                    Return = "Matrix")
{
  if (!is.null(Pedigree) & !is.null(Pairs))
    stop("Please provide Pedigree or Pairs, not both")
  if (is.null(Pedigree) & is.null(Pairs))
    stop("Please provide Pedigree or Pairs")

  if (!is.null(Pedigree)) {
    if (!(Return %in% c("Matrix", "Array")))   stop("Return must be 'Matrix' or 'Array'")
    Pedigree <- PedPolish(Pedigree, ZeroToNA=TRUE, NullOK=FALSE)
    if (GenBack==2) {
    nGPcols <- length(intersect(c("MGM","MGF","PGM","PGF"), names(Pedigree)))
    if (nGPcols > 0 & nGPcols < 4) {
      stop("Pedigree must either have none of the columns 'MGM', 'MGF', 'PGM', 'PGF', or all 4")
    } else if (nGPcols == 0) {
      Pedigree <- GPcols(Pedigree)
    }
  }
    RelA <- GetRelA(Pedigree, GenBack = GenBack, patmat = patmat)

    if (Return == "Matrix" ) {
      rel.i <- which(RelA == 1, arr.ind=TRUE)
      rel.i <- rel.i[!duplicated(rel.i[,1:2]), ]  # 'highest' relationship only per pair
      RelNames <- dimnames(RelA)[[3]]
      RelM <- matrix("U", dim(RelA)[1], dim(RelA)[2], dimnames = list(Pedigree$id, Pedigree$id))
      RelM[rel.i[,1:2]] <- RelNames[rel.i[,3]]
      return( RelM )
    } else {
      return( RelA )
    }

  } else if (!is.null(Pairs)) {
    if (Return != "Matrix")   stop("When providing Pairs, Return must be 'Matrix'")
    if (!class(Pairs) %in% c("data.frame", "matrix"))
      stop("Pairs should be a dataframe or matrix")
    Pairs <- as.data.frame(Pairs)
    names(Pairs)[1:3] <- c("ID1", "ID2", "Rel")
    for (x in 1:3)  Pairs[,x] <- as.character(Pairs[,x])

    RelM.tmp <- plyr::daply(Pairs, .variables=c("ID1", "ID2"), .fun=function(df) df$Rel)
    # make into symmetrical matrix to get consistency in above/below diagonal:
    IDs <- unique(c(as.character(Pairs$ID1), as.character(Pairs$ID2)))
    RelM.a <- inflate(RelM.tmp, IDs)
    RelM.b <- inflate(t(RelM.tmp), IDs)
    if(any(RelM.a != RelM.b, na.rm=TRUE)) {
      stop("One or more pairs occur 2x in 'Pairs', with different relationship")
    }
    RelM <- RelM.a
    RelM[,] <- "U"
    diag(RelM) <- "S"
    RelM[!is.na(RelM.a)] <- RelM.a[!is.na(RelM.a)]
    RelM[!is.na(RelM.b)] <- RelM.b[!is.na(RelM.b)]
    return( RelM )
  }

}





################################################################################
################################################################################

#' @title Array with Pairwise Relationships
#'
#' @description Generate an array indicating the relationship(s) between all
#'  pairs of individuals according to the pedigree.
#'
#' @param Ped  dataframe with columns id - dam - sire.
#' @param GenBack  number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grand-parental,
#'   avuncular and first cousins.
#' @param patmat  logical, distinguish between paternal versus maternal relative
#'   pairs? For avuncular pairs, the distinction is never made.
#'
#' @return a 3D array indicating if the pair has the specified relationship (1)
#' or not (0). The various relationship considered are in the 3rd dimension:
#'  \item{M}{}
#'  \item{P}{}
#'  \item{FS}{full siblings, including double 'other half sibs'}
#'  \item{MS}{}
#'  \item{PS}{}
#'  \item{XS}{other sibs: mother of A is father of B, or vv}
#' etc.
#'
#' @keywords internal

GetRelA <- function(Ped = NULL, GenBack = 1, patmat = TRUE)
{
  if (is.null(Ped))  stop("Please provide Pedigree")
  nInd <- nrow(Ped)

  if (GenBack == 1) {
    Anc <- c("dam", "sire")
  } else if (GenBack == 2) {
    Anc <- c("dam", "sire", "MGM","MGF","PGM","PGF")
  }
  PedX <- matrix(NA, nInd, length(Anc), dimnames=list(Ped$id, Anc))
  for (p in Anc) {
    PedX[,p] <- sapply(Ped[,p], function(x) ifelse(is.na(x), NA,
                                                        which(Ped[,"id"]==x)))
  }

  Rels <- c("S", "M", "P", "O", "FS", "MHS", "PHS", "XHS")  # XHS: 'cross'-half-sibs
  if (GenBack == 2) {
    Rels <- c(Rels, "MGM", "MGF", "PGM", "PGF", "GO",
    "FA", "FN", "HA", "HN", "DFC1", "FC1")
  }

  RelA <- array(0, dim = c(nInd, nInd, length(Rels)),
                dimnames=list(Ped$id, Ped$id, Rels))

  # self
  RelA[,,"S"][cbind(seq_len(nInd), seq_len(nInd))] <- 1

  # parent-offspring
  RelA[,,"M"][cbind(seq_len(nInd), PedX[,"dam"])] <- 1
  RelA[,,"P"][cbind(seq_len(nInd), PedX[,"sire"])] <- 1
  RelA[,,"O"][cbind(c(PedX[,"dam"], PedX[,"sire"]), rep(seq_len(nInd),2))] <- 1

  # grandparent-grandoffspring
  if (GenBack == 2) {
    for (a in c("MGM","MGF","PGM","PGF")) {
      RelA[,,a][cbind(seq_len(nInd), PedX[,a])] <- 1
      RelA[,,"GO"][cbind(PedX[,a], seq_len(nInd))] <- 1
    }
  }

  # siblings
  RelA[,,"FS"][(outer(PedX[,"dam"], PedX[,"dam"], "==") &
                  outer(PedX[,"sire"], PedX[,"sire"], "==")) |
                 (outer(PedX[,"dam"], PedX[,"sire"], "==") &
                    outer(PedX[,"sire"], PedX[,"dam"], "=="))] <- 1
  RelA[,,"MHS"][outer(PedX[,"dam"], PedX[,"dam"], "==")] <- 1
  RelA[,,"PHS"][outer(PedX[,"sire"], PedX[,"sire"], "==")] <- 1
  RelA[,,"XHS"][outer(PedX[,"dam"], PedX[,"sire"], "==") |
                  outer(PedX[,"sire"], PedX[,"dam"], "==")] <- 1   # hermaphrodites

  if (GenBack == 2) {
    # avuncular
    RelA[,,"FA"][(outer(PedX[,"MGM"], PedX[,"dam"], "==") &
                  outer(PedX[,"MGF"], PedX[,"sire"], "==")) |
                 (outer(PedX[,"PGM"], PedX[,"dam"], "==") &
                    outer(PedX[,"PGF"], PedX[,"sire"], "==")) |
                 (outer(PedX[,"MGM"], PedX[,"sire"], "==") &     # hermaphrodites
                    outer(PedX[,"MGF"], PedX[,"dam"], "==")) |
                 (outer(PedX[,"PGM"], PedX[,"sire"], "==") &    # hermaphrodites
                    outer(PedX[,"PGF"], PedX[,"dam"], "=="))] <- 1
    RelA[,,"FN"] <- t(RelA[,,"FA"])
    RelA[,,"HA"][outer(PedX[,"MGM"], PedX[,"dam"], "==") |
                   outer(PedX[,"MGF"], PedX[,"sire"], "==") |
                   outer(PedX[,"PGM"], PedX[,"dam"], "==") |
                   outer(PedX[,"PGF"], PedX[,"sire"], "==") |
                   outer(PedX[,"MGM"], PedX[,"sire"], "==") |     # hermaphrodites
                   outer(PedX[,"MGF"], PedX[,"dam"], "==") |
                   outer(PedX[,"PGM"], PedX[,"sire"], "==") |
                   outer(PedX[,"PGF"], PedX[,"dam"], "==")] <- 1
    RelA[,,"HN"] <- t(RelA[,,"HA"])

    # (double) full 1st cousins
    RelA[,,"DFC1"][(outer(PedX[,"MGM"], PedX[,"MGM"], "==") &
                      outer(PedX[,"MGF"], PedX[,"MGF"], "==") &
                      outer(PedX[,"PGM"], PedX[,"PGM"], "==") &
                      outer(PedX[,"PGF"], PedX[,"PGF"], "==")) |
                     (outer(PedX[,"MGM"], PedX[,"PGM"], "==") &
                        outer(PedX[,"MGF"], PedX[,"PGF"], "==") &
                        outer(PedX[,"PGM"], PedX[,"MGM"], "==") &
                        outer(PedX[,"PGF"], PedX[,"MGF"], "=="))] <- 1
    RelA[,,"FC1"][(outer(PedX[,"MGM"], PedX[,"MGM"], "==") &
                     outer(PedX[,"MGF"], PedX[,"MGF"], "==")) |
                    (outer(PedX[,"PGM"], PedX[,"PGM"], "==") &
                       outer(PedX[,"PGF"], PedX[,"PGF"], "==")) |
                    (outer(PedX[,"MGM"], PedX[,"PGM"], "==") &
                       outer(PedX[,"MGF"], PedX[,"PGF"], "==")) |
                    (outer(PedX[,"PGM"], PedX[,"MGM"], "==") &
                       outer(PedX[,"PGF"], PedX[,"MGF"], "=="))] <- 1
    # cousins not implemented for hermaphrodites yet.
  }

  # some exceptions
  if (GenBack == 2) {
    RelA[,,"FC1"][RelA[,,"MHS"]==1 | RelA[,,"PHS"]==1 | RelA[,,"XHS"]==1] <- 0
    RelA[,,"HA"][RelA[,,"FA"]==1] <- 0
    RelA[,,"FA"][RelA[,,"M"]==1 | RelA[,,"P"]==1] <- 0
  }
  for (r in c("FS", "MHS", "PHS", "XHS")) {
    RelA[,,r][RelA[,,"S"]==1] <- 0
    if (r=="FS")  next
    RelA[,,r][RelA[,,"FS"]==1] <- 0
  }

  if (!patmat) {
    RelsB <- c("S", "MP", "O", "FS", "HS")
    if (GenBack == 2) {
      RelsB <- c(RelsB, "GP", "GO",
                 "FA", "FN", "HA", "HN", "DFC1", "FC1")
    }
    RelA.B <- array(0, dim = c(nInd, nInd, length(RelsB)),
                    dimnames=list(Ped$id, Ped$id, RelsB))
    RelA.B[,,"MP"] <- RelA[,,"M"] | RelA[,,"P"]
    RelA.B[,,"HS"] <- RelA[,,"MHS"] | RelA[,,"PHS"] | RelA[,,"XHS"]
    if (GenBack == 2) {
      RelA.B[,,"GP"] <- RelA[,,"MGM"] | RelA[,,"MGF"] | RelA[,,"PGM"] | RelA[,,"PGF"]
    }
    for (r in RelsB) {
      if (r %in% c("MP", "HS", "GP"))  next
      RelA.B[,,r] <- RelA[,,r]
    }
    RelA <- RelA.B
  }

  return( RelA )
}


