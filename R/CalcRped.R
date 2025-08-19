#' @title Calculate Pedigree Relatedness
#'
#' @description Morph pedigree into a \pkg{kinship2} compatible format and use
#'   \code{\link[kinship2]{kinship}} to calculate kinship coefficients;
#'   relatedness = 2*kinship.
#'
#' @param Pedigree dataframe with columns id-dam-sire.
#' @param OUT  desired output format, 'M' for matrix or 'DF' for dataframe with
#'   columns IID1 - IID2 - R.ped.
#'
#' @return A matrix or dataframe.
#'
#' @export

CalcRped <- function(Pedigree, OUT="DF")
{
  if (!OUT %in% c("M", "DF"))
    stop("'OUT' must be 'M' (matrix) or 'DF' (data.frame)")

  if (!requireNamespace("kinship2", quietly = TRUE)) {
    if (interactive())  cli::cli_alert_info("Installing package {.pkg kinship2}... ")
    utils::install.packages("kinship2")
  }

  PedP <- sequoia::PedPolish(Pedigree, DropNonSNPd=FALSE, FillParents = TRUE)
  PedP$Sex <- with(PedP, ifelse(id %in% dam, "female",  "male"))
  # default to 'male' to avoid warning: "More than 25% of the gender values are 'unknown'"

  Ped.fix <- with(PedP, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                             sex=Sex))
  Ped.k <- with(Ped.fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))
  kin.M <- kinship2::kinship(Ped.k)
  kin.M <- kin.M[Pedigree$id, Pedigree$id]   # drop dummies created by FillParents=TRUE
  if (OUT == "M") {
    return( 2*kin.M )

  } else {
    kin.df <- data.frame("IID1" = rep(rownames(kin.M), each=ncol(kin.M)),
                         "IID2" = rep(colnames(kin.M), times=nrow(kin.M)),
                         "R.ped" = 2*c(t(kin.M)),
                         stringsAsFactors=FALSE)
    return( kin.df )
  }
}
