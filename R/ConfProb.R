#============================================================================
#============================================================================
#' @title Estimate confidence probability
#'
#' @description Estimate the assignment error rate by repeatedly simulating
#' data from
#' a reference pedigree using \code{\link{SimGeno}}, reconstruction a pedigree
#' from this using \code{\link{sequoia}}, and counting the number of mismatches
#' using \code{\link{PedCompare}}.
#'
#' @details The confidence probability is taken as the number of correct
#' (matching)
#' assignments, divided by all assignments made. A confidence of '1' should be
#' interpreted as '> 1 - 1/(sum(!is.na(Ped$dam)) * nSim)'
#'
#' @param Ped Reference pedigree from which to simulate, dataframe with
#'  columns id-dam-sire. Additional columns are ignored
#' @param LifeHistData Dataframe with id, sex (1=female, 2=male, 3=unknown),
#'  and birth year.
#' @param Specs Parameter values for running sequoia, as named vector.
#' @param Full Full pedigree reconstruction (TRUE) or only parentage assignment
#'   (FALSE)
#' @param nSim number of simulations to perform.
#' @param ParMis proportion of parents assumed to have a fully missing genotype.
#' @param args.sim  list of additional arguments to pass to \code{\link{SimGeno}}
#' @param return.PC  return all \code{\link{PedCompare}} \code{Counts}?
#' @param quiet suppress messages. `very' also suppresses simulation counter
#'
#' @return A 2x2 matrix for parentage assignment, or a 2x7x2 array for
#'  full pedigree reconstruction, with for dams and sires and per category (see
#'  \code{\link{PedCompare}}) the average and minimum number of Match/(Match +
#'  Mismatch + P2only).
#'
#' When return.PC is TRUE, a list is returned with two arrays: ConfProb
#' contains the average confidence probability across simulations, and
#' SimCounts all counts of matches, mismatches, Pedigree1-only and pedigree2-
#' only per simulation.
#'
#' @examples
#' \dontrun{
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#' SeqOUT <- sequoia(GenoM = SimGeno_example,
#'                   LifeHistData = LH_HSg5, MaxSibIter = 0)
#' ConfPr <- EstConf(Ped = SeqOUT$PedigreePar, LifeHistData = LH_HSg5,
#'                  Specs = SeqOUT$Specs, Full = FALSE, nSim = 10)
#' }
#'
#' @export

EstConf <- function(Ped = NULL,
                    LifeHistData = NULL,
                    Specs = NULL,
                    Full = TRUE,
                    # todo: use previous ageprior
                    nSim = 10,
                    ParMis = 0.4,
                    args.sim = NULL,
                    return.PC = FALSE,
                    quiet=TRUE)
{
  PedA <- Ped[,1:3]
  SpecsA <- Specs
  seq.quiet <- ifelse(is.logical(quiet), quiet,
                    ifelse(quiet == "very", TRUE, FALSE))

  if (is.null(Specs))  stop("Please provide Specs")
  if (is.null(Ped))  stop("Please provide Ped")
  if (Full & Specs["MaxSibIter"]==0)  stop("Please specify Specs['MaxSibIter'] >0 or Full=FALSE")
  if (any(substr(PedA$id,1,2) %in% c("F0", "M0"))) {
    PedA$id <- addpref(PedA$id)
    PedA$dam <- addpref(PedA$dam)
    PedA$sire <- addpref(PedA$sire)
  }
  SpecsA["FindMaybeRel"] <- FALSE
  SpecsA["CalcLLR"] <- FALSE
  if (!Full) SpecsA["MaxSibIter"] <- 0

  if (!Full) {
    PC <- array(dim=c(nSim, 5, 2),
                dimnames=list(1:nSim, c("Total", "Match", "Mismatch", "P1only", "P2only"),
                              c("dam", "sire")))
  } else {
    PC <- array(dim=c(nSim, 7, 5, 2),
                  dimnames=list(1:nSim, c("GG", "GD", "GT", "DG", "DD", "DT", "TT"),
                                c("Total", "Match", "Mismatch", "P1only", "P2only"),
                                 c("dam", "sire")))
  }
  for (i in 1:nSim) {
    if (quiet != "very")  cat("\n\n i=", i, "\t", format(Sys.time(), "%H:%M:%S"), "\t\t")
    GM <- SimGeno(Ped = PedA, nSnp = SpecsA[, "NumberSnps"],
                  ParMis=ParMis, quiet=seq.quiet)  # TODO: pass list of arguments
    Seq.i <- sequoia(GenoM = GM,
                     LifeHistData = LifeHistData,
                     SeqList = list(Specs = SpecsA),
                     MaxSibIter = SpecsA[, "MaxSibIter"],
                     quiet = seq.quiet)
    if (!Full) {
      PC[i,,] <- PedCompare(Ped1 = PedA, Ped2 = Seq.i$PedigreePar)$Counts["GG",,]
      cat("AR: ", round(sum(PC[i,"Match",])/sum(PC[i,"Total",]),3), "\t",
          "# ER:", sum(PC[i,c("Mismatch", "P2only"),]), "\n")
    } else {
      PC[i,,,] <- PedCompare(Ped1 = PedA, Ped2 = Seq.i$Pedigree)$Counts
      cat("AR: ", round(sum(PC[i,"TT","Match",])/sum(PC[i,"TT","Total",]),3), "\t",
          "# ER:", sum(PC[i,"TT",c("Mismatch", "P2only"),]), "\n")
    }
  }

  if (!Full) {
    ECP <- matrix(NA,2,2, dimnames=list(c("dam", "sire"), c("mean","min")))
    for (s in 1:2) {
      tmp <- PC[,"Match",s]/apply(PC[,c("Match", "Mismatch", "P2only"),s],1,sum)
      ECP[s, "mean"] <- mean(tmp)
      ECP[s, "min"] <- min(tmp)
    }
    ntot <- sum(PC[,"Total",s])
  } else {
    ECP <- array(dim=c(2,7,2), dimnames=list(c("dam", "sire"),
                                          c("GG","GD","GT","DG","DD","DT","TT"),
                                          c("mean","min")))
    for (s in 1:2) {
      tmp <- PC[,,"Match",s]/apply(PC[,,c("Match", "Mismatch", "P2only"),s],c(1:2),sum)
      ECP[s, , "mean"] <- apply(tmp,2,mean)
      ECP[s, , "min"] <- apply(tmp,2,min)
    }
    ntot <- sum(PC[,"TT","Total",s])
  }
  ECP <- round(ECP, nchar(ntot))
  if (return.PC) {
    list(SimCounts = PC, ConfProb = ECP)
  } else {
    ECP
  }
}


#============================================================================

addpref <- function(x, pf="a") {
  y <- ifelse(is.na(x), x,
          ifelse(substr(x,1,2)=="F0" | substr(x,1,2)=="M0",
                 paste0(pf, x), x))
  y
}
