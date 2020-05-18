#' @title Pairwise relationship
#'
#' @description Determine the relationship between individual X and all other
#'   individuals in the pedigree, going up to 1 or 2 generations back.
#'
#' @param x The focal individual, either its rownumber in the pedigree or ID.
#' @param Pedigree  dataframe columns id - dam - sire.
#' @param GenBack  Number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grandparental,
#'   avuncular and first cousins.
#' @param patmat  logical, distinguish between paternal versus maternal relative
#'   pairs?
#'
#' @return A named vector of length equal to the number of rows in Ped, with for
#'   each ID its relationship to the focal individual:
#'    \item{S}{Self}
#'    \item{M}{Mother}
#'    \item{P}{Father}
#'    \item{O}{Offspring}
#'    \item{FS}{Full sibling}
#'    \item{MHS}{Maternal half-sibling}
#'    \item{PHS}{Paternal half-sibling}
#'    \item{MGM}{Maternal grandmother}
#'    \item{MGF}{Maternal grandfather}
#'    \item{PGM}{Paternal grandmother}
#'    \item{PGF}{Paternal grandfather}
#'    \item{GO}{Grand-offspring}
#'    \item{FA}{Full avuncular; maternal or paternal aunt or uncle}
#'    \item{HA}{Half avuncular}
#'    \item{FN}{Full nephew/niece}
#'    \item{HN}{Half nephew/niece}
#'    \item{FC1}{Full first cousin}
#'    \item{DFC1}{Double full first cousin}
#'    \item{U}{Unrelated (or otherwise related)}
#'
#' @seealso \code{\link{ComparePairs}} to compare pairwise relationships
#'   between 2 pedigrees.
#'
#' @examples
#' data(Ped_griffin)
#' # find all relatives of a specific individual
#' Rel42 <- GetRelCat("i042_2003_F", Ped_griffin)
#' Rel42[Rel42 != "U"]
#'
#' # make NxN matrix with relationship categories:
#' Ped_griffin_sub <- Ped_griffin[Ped_griffin$birthyear<2003,]  # quicker
#' RCM <- sapply(seq_along(Ped_griffin_sub$id), GetRelCat, Ped_griffin_sub)
#' table(RCM)
#' #   M  MHS    O    P    S    U
#' #  10    6   16    6   40 1522
#' # note that sibling & cousin pairs are counted twice!
#' # Parent-offspring pairs are counted directionally:
#' # once as offspring (O), once as mother (M) or father (P)
#'
#' # for large pedigrees, table(factor()) is much faster:
#' table(factor(RCM, levels=c("M","P","FS","MHS","PHS","U")))
#'
#' # list the maternal half-siblings:
#' these <- which(RCM=="MHS", arr.ind=TRUE)
#' data.frame(id1 = Ped_griffin_sub$id[these[,1]],
#'            id2 = Ped_griffin_sub$id[these[,2]])
#'
#'
#' # Get Colony-style lists of full sibs & half sibs dyads:
#' \dontrun{
#' RCM <- sapply(seq_along(MyPedigree$id), GetRelCat, Pedigree = MyPedigree,
#'               GenBack = 1, patmat = FALSE)
#' # rownumbers of pairs of FS & HS
#' FullSibDyads <- which(RCM == "FS", arr.ind=TRUE)
#' HalfSibDyads <- which(RCM == "HS", arr.ind=TRUE)
#'
#' # each pair is listed 2x - fix:
#' FullSibDyads <- FullSibDyads[FullSibDyads[,1] < FullSibDyads[,2], ]
#' HalfSibDyads <- HalfSibDyads[HalfSibDyads[,1] < HalfSibDyads[,2], ]
#'
#' # translate rownumbers into IDs
#' MyPedigree$id <- as.character(MyPedigree$id)
#' FullSibDyads <- cbind(MyPedigree$id[FullSibDyads[,1]],
#'                       MyPedigree$id[FullSibDyads[,2]])
#' HalfSibDyads <- cbind(MyPedigree$id[HalfSibDyads[,1]],
#'                       MyPedigree$id[HalfSibDyads[,2]])
#' }
#'
#' @importFrom stats setNames
#'
#' @export

GetRelCat <- function(x, Pedigree, GenBack=2, patmat=TRUE) {
  Ped <- Pedigree
  if (is.character(x))   x <- which(Ped$id == x)
  if (length(intersect(c("id","dam","sire"), tolower(names(Ped)))) < 3) {
    stop("Expect Pedigree column names to contain 'id', 'dam' and 'sire'")
  } else {
    for (i in 1:ncol(Ped)) {
      if (tolower(names(Ped)[i]) %in% c("id","dam","sire")) {
        names(Ped)[i] <- tolower(names(Ped)[i])
      }
    }
  }
  if (!GenBack %in% c(1,2)) {
    stop(paste("Expecting value '1' or '2' for GenBack (no. generations back), got ", GenBack))
  }
  if (!patmat %in% c(TRUE,FALSE))  stop("'patmat' must be TRUE or FALSE")

  for (i in 1:ncol(Ped))  Ped[,i] <- as.character(Ped[,i])
  if (GenBack==2) {
    nGPcols <- length(intersect(c("MGM","MGF","PGM","PGF"), names(Ped)))
    if (nGPcols > 0 & nGPcols < 4) {
      stop("Pedigree must either have none of the columns 'MGM', 'MGF', 'PGM', 'PGF', or all 4")
    } else if (nGPcols == 0) {
      IDorder <- Ped$id   # merge() ignores sort=FALSE
      Ped <- merge(Ped, setNames(Ped[Ped$id %in% Ped$dam,], c("dam", "MGM", "MGF")), all.x=TRUE)
      Ped <- merge(Ped, setNames(Ped[,c("id","dam","sire")], c("sire", "PGM", "PGF")), all.x=TRUE)
      rownames(Ped) <- Ped$id
      Ped <- Ped[IDorder, ]
    }
  }

  if (patmat) {
    RCV <- with(Ped,
            ifelse(id == id[x], "S",
             ifelse(eqv(dam[x],id,FALSE), "M",
              ifelse(eqv(sire[x], id,FALSE), "P",
               ifelse(eqv(id[x], dam, FALSE) | eqv(id[x], sire, FALSE), "O",
                ifelse(eqv(dam[x],dam,FALSE) & eqv(sire[x], sire,FALSE), "FS",
                 ifelse(eqv(dam[x],dam,FALSE), "MHS",
                  ifelse(eqv(sire[x], sire,FALSE), "PHS",
                    "U"))))))))
  } else {
    RCV <- with(Ped,
            ifelse(id == id[x], "S",
             ifelse(eqv(dam[x],id,FALSE) | eqv(sire[x], id,FALSE), "MP",
               ifelse(eqv(id[x], dam, FALSE) | eqv(id[x], sire, FALSE), "O",
                ifelse(eqv(dam[x],dam,FALSE) & eqv(sire[x], sire,FALSE), "FS",
                 ifelse(eqv(dam[x],dam,FALSE) | eqv(sire[x], sire,FALSE), "HS",
                        "U"))))))
  }

  if (GenBack==2) {
    RCV <- with(Ped,
         ifelse(RCV != "U", RCV,
            ifelse(eqv(MGM[x],id,FALSE), "MGM",
             ifelse(eqv(MGF[x],id,FALSE), "MGF",
              ifelse(eqv(PGM[x],id,FALSE), "PGM",
               ifelse(eqv(PGF[x],id,FALSE), "PGF",
                ifelse(eqv(id[x], MGM, FALSE) | eqv(id[x], MGF, FALSE) |
                         eqv(id[x], PGM, FALSE) | eqv(id[x], PGF, FALSE), "GO",
                ifelse((eqv(MGM[x],dam,FALSE) & eqv(MGF[x], sire, FALSE)) |
                          (eqv(PGM[x],dam,FALSE) & eqv(PGF[x], sire, FALSE)), "FA", # full avuncular
                ifelse((eqv(dam[x], MGM,FALSE) & eqv(sire[x], MGF, FALSE)) |
                          (eqv(dam[x], PGM,FALSE) & eqv(sire[x], PGF, FALSE)), "FN", # full niece/nephew
                 ifelse(eqv(MGM[x],dam,FALSE) | eqv(MGF[x], sire, FALSE) |
                           eqv(PGM[x],dam,FALSE) | eqv(PGF[x], sire, FALSE), "HA",  # half avuncular
                 ifelse((eqv(dam[x], MGM,FALSE) | eqv(sire[x], MGF, FALSE)) |
                          (eqv(dam[x], PGM,FALSE) | eqv(sire[x], PGF, FALSE)), "HN",
                  ifelse((eqv(MGM[x],MGM,FALSE) & eqv(MGF[x], MGF, FALSE) &
                        eqv(PGM[x],PGM,FALSE) & eqv(PGF[x], PGF, FALSE)) |
                         (eqv(MGM[x],PGM,FALSE) & eqv(MGF[x], PGF, FALSE) &
                          eqv(PGM[x],MGM,FALSE) & eqv(PGF[x], MGF, FALSE)), "DFC1",  # double full 1st cousins
                   ifelse((eqv(MGM[x],MGM,FALSE) & eqv(MGF[x], MGF, FALSE)) |
                        (eqv(PGM[x],PGM,FALSE) & eqv(PGF[x], PGF, FALSE)) |
                         (eqv(MGM[x],PGM,FALSE) & eqv(MGF[x], PGF, FALSE)) |
                          (eqv(PGM[x],MGM,FALSE) & eqv(PGF[x], MGF, FALSE)), "FC1",  # full 1st cousins
                          "U")))))))))))))
    if (!patmat) {
      RCV[RCV %in% c("MGM", "MGF", "PGM", "PGF")] <- "GP"
    }
  }

  names(RCV) <- Ped$id
	return( RCV )
}
