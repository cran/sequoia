#' @title Pairwise relationship
#'
#' @description Determine the relationship between individual X and all other
#'   individuals in the pedigree, going up to 1 or 2 generations back.
#'
#' @details if
#'
#' @param x The focal individual, either its rownumber in the pedigreee or ID.
#' @param Ped  dataframe with pedigree; columns id - dam - sire.
#' @param GenBack  Number of generations back to consider; 1 returns
#'   parent-offspring and sibling relationships, 2 also returns grandparental,
#'   avuncular and first cousins.
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
#' @examples
#' \dontrun{
#' # make NxN matrix with relationship categories:
#' RCM <- sapply(1:nrow(Ped), GetRelCat, Ped)
#' # table is considerably faster on a factor than a character vector:
#' table(factor(RCM, levels=c("M","P","FS","MHS","PHS","MGM")))
#' # note that sibling & cousin pairs are counted twice!
#'}
#'
#' @importFrom stats setNames
#'
#' @export

GetRelCat <- function(x, Ped, GenBack=2) {
  for (i in 1:ncol(Ped)) Ped[,i] <- as.character(Ped[,i])
  if (length(intersect(c("id","dam","sire"), names(Ped))) < 3) {
    stop("Expect column names to contain 'id', 'dam' and 'sire'")
  }
  if (!GenBack %in% c(1,2)) {
    stop(paste("Expecting value '1' or '2' for GenBack (no. generations back), got ", GenBack))
  }

  if (GenBack==2) {
    if (ncol(Ped)==3) {
      IDorder <- Ped$id   # merge() ignores sort=FALSE
      Ped <- merge(Ped, setNames(Ped[Ped$id %in% Ped$dam,], c("dam", "MGM", "MGF")), all.x=TRUE)
      Ped <- merge(Ped, setNames(Ped[,c("id","dam","sire")], c("sire", "PGM", "PGF")), all.x=TRUE)
      rownames(Ped) <- Ped$id
      Ped <- Ped[IDorder, ]
    }
    if (length(intersect(c("MGM","MGF","PGM","PGF"), names(Ped))) < 4) {
      stop('Expect column names to contain "MGM","MGF","PGM", and "PGF" when GenBack=2 \n and Ped has more than 3 columns')
    }
  }
  if (is.character(x))   x <- which(Ped$id == x)

  RCV <- with(Ped,
          ifelse(id == id[x], "S",
           ifelse(eqv(dam[x],id,FALSE), "M",
            ifelse(eqv(sire[x], id,FALSE), "P",
             ifelse(eqv(id[x], dam, FALSE) | eqv(id[x], sire, FALSE), "O",
              ifelse(eqv(dam[x],dam,FALSE) & eqv(sire[x], sire,FALSE), "FS",
               ifelse(eqv(dam[x],dam,FALSE), "MHS",
                ifelse(eqv(sire[x], sire,FALSE), "PHS",
                  "U"))))))))

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
  }

  names(RCV) <- Ped$id
	return( RCV )
}
