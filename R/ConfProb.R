#============================================================================
#============================================================================
#' @title Estimate confidence probability
#'
#' @description Estimate the assignment error rate by repeatedly simulating data
#'   from a reference pedigree using \code{\link{SimGeno}}, reconstruction a
#'   pedigree from this using \code{\link{sequoia}}, and counting the number of
#'   mismatches using \code{\link{PedCompare}}.
#'
#' @details The confidence probability is taken as the number of correct
#'   (matching) assignments, divided by all assignments made. A confidence of
#'   '1' should be interpreted as '> 1 - 1/(sum(!is.na(Ped$dam)) * nSim)'
#'
#' @param Ped Reference pedigree from which to simulate, dataframe with columns
#'   id-dam-sire. Additional columns are ignored
#' @param LifeHistData Dataframe with id, sex (1=female, 2=male, 3=unknown), and
#'   birth year.
#' @param args.sim  list of arguments to pass to \code{\link{SimGeno}}, such as
#'   \code{nSnp} (number of SNPs), \code{SnpError} (genotyping error rate) and
#'   \code{ParMis} (proportion of non-genotyped parents). Set to NULL to use all
#'   default values.
#' @param args.seq  list of arguments to pass to \code{\link{sequoia}}, such as
#'   \code{MaxSibIter} (max no. sibship clustering iterations, '0' for parentage
#'   assignment only) and \code{Err} (assumed genotyping error rate). May
#'   include (part of) SeqList, the list of sequoia output (i.e. as a
#'   list-within-a-list). Set to NULL to use all default values.
#' @param nSim number of rounds of simulations to perform.
#' @param return.PC  return all \code{\link{PedCompare}} \code{Counts}?
#' @param quiet suppress messages. `very' also suppresses simulation counter,
#'   TRUE merely runs SimGeno and sequoia quietly.
#'
#' @return When \code{return.PC = FALSE}, a 2x2 matrix for parentage assignment, or a
#'   2x7x2 array for full pedigree reconstruction, with for dams and sires and
#'   per category (see \code{\link{PedCompare}}) the average and minimum number
#'   of Match/(Match + Mismatch + P2only).
#'
#'   When \code{return.PC} is TRUE, a list is returned with:
#'   \item{ConfProb}{Average confidence probability across simulations, as
#'     returned when \code{return.PC = FALSE}.}
#'   \item{SimCounts}{All counts of matches, mismatches, Pedigree1-only and
#'     pedigree2-only, per simulation.}
#'   \item{RunParams}{Current call to EstConf, as well as the default
#'     parameter values for \code{EstConf, SimGeno}, and \code{sequoia}.}
#'   \item{RunTime}{\code{sequoia} runtime per simulation in seconds, as
#'     measured by \code{\link{system.time}()['elapsed']}.}
#'
#' @seealso \code{\link{SimGeno}, \link{sequoia}, \link{PedCompare}}
#'
#' @examples
#' \dontrun{
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#'
#' conf.A <- EstConf(Ped = Ped_HSg5, LifeHistData = LH_HSg5,
#'    args.sim = list(nSnp = 100, SnpError = 5e-3, ParMis=c(0.2, 0.5)),
#'    args.seq = list(MaxSibIter = 0, Err=1e-4, Tassign=0.5),
#'    nSim = 3, return.PC = TRUE)
#'
#' # effect of tweaking AgePriors
#' # (only some effect due to low no. SNPs & high error rate,
#' #  effect of increasing no. SNPs is much larger)
#' AP <- MakeAgePrior(Ped = Ped_HSg5, LifeHistData = LH_HSg5,
#'                    Flatten = FALSE, Smooth = FALSE)
#' conf.B <- EstConf(Ped = Ped_HSg5, LifeHistData = LH_HSg5,
#'    args.sim = list(nSnp = 100, SnpError = 5e-3, ParMis=c(0.2, 0.5)),
#'    args.seq = list(MaxSibIter = 0, Err=1e-4, Tassign=0.5,
#'                    SeqList = list(AgePriors = AP)),
#'    nSim = 3, return.PC = TRUE)
#'
#' # with sibship clustering
#' conf.C <- EstConf(Ped = Ped_HSg5, LifeHistData = LH_HSg5,
#'    args.sim = list(nSnp = 200, SnpError = 5e-3, ParMis=c(0.2, 0.5)),
#'    args.seq = list(MaxSibIter = 10, Err=1e-4, Tassign=0.5),
#'    nSim = 3, return.PC = TRUE)
#' conf.C$ConfProb[,"GG",]  # Genotyped individuals, Genotyped parent
#' conf.C$ConfProb[,"GD",]  # Genotyped individuals, Dummy parent
#' AR <- apply(conf.C$SimCounts, 1, function(M) M["TT","Match", ]/M["TT","Total", ])
#' ER <- apply(conf.C$SimCounts, 1,
#'        function(M) (M["TT","Mismatch", ] + M["TT","P2only", ])/M["TT","Total", ])
#' apply(ER, 1, mean)  # separate error rate dams & sires
#' mean(ER)            # overall error rate
#' }
#'
#' @export

EstConf <- function(Ped = NULL,
                    LifeHistData = NULL,
                    args.sim = list(nSnp = 400, SnpError = 1e-3, ParMis=c(0.4, 0.4)),
                    args.seq = list(MaxSibIter = 10, Err=1e-4, Tassign=0.5),   # todo: pass ageprior
                    nSim = 10,
                    return.PC = FALSE,
                    quiet=TRUE)
{
  if (is.null(Ped))  stop("Please provide Ped")
  if (is.null(LifeHistData))  stop("Please provide LifeHistData")
  if (!is.null(args.sim) & !is.list(args.sim))  stop("args.sim should be a list or NULL")
  if (!is.null(args.seq) & !is.list(args.seq))  stop("args.seq should be a list or NULL")

  PedA <- Ped[,1:3]
  if (any(substr(PedA$id,1,2) %in% c("F0", "M0"))) {
    PedA$id <- addpref(PedA$id)
    PedA$dam <- addpref(PedA$dam)
    PedA$sire <- addpref(PedA$sire)
  }

  args.seq["FindMaybeRel"] <- FALSE
  args.seq["CalcLLR"] <- FALSE
  if ("MaxSibIter" %in% names(args.seq)) {
    ParSib <- ifelse(args.seq$MaxSibIter > 0, "sib", "par")
  } else {
    ParSib <- "sib"   # default MaxSibIter=10
  }
  if (quiet != "very") {
    if (ParSib == "par") {
      message("MaxSibIter=0: Simulating parentage assignment only ...")
      PC <- array(dim=c(nSim, 5, 2),
                dimnames=list(1:nSim, c("Total", "Match", "Mismatch", "P1only", "P2only"),
                              c("dam", "sire")))
    } else {
      message("MaxSibIter>0: Simulating full pedigree reconstruction ...")
      PC <- array(dim=c(nSim, 7, 5, 2),
                  dimnames=list(1:nSim, c("GG", "GD", "GT", "DG", "DD", "DT", "TT"),
                                c("Total", "Match", "Mismatch", "P1only", "P2only"),
                                 c("dam", "sire")))
    }
  }
  utils::flush.console()
  seq.quiet <- ifelse(is.logical(quiet), quiet,
                    ifelse(quiet == "very", TRUE, FALSE))

  #~~~~~~~~~~~
  RunTime <- rep(NA, nSim)
  for (i in 1:nSim) {
    if (quiet != "very")  cat("\n i=", i, "\t", format(Sys.time(), "%H:%M:%S"), "\t\t")

    GM <- do.call(SimGeno, c(list(Ped=Ped), args.sim))

    RunTime[i] <- system.time(Seq.i <- do.call(sequoia, c(list(GenoM = GM,
                                     LifeHistData = LifeHistData,
                                     quiet = seq.quiet),
                                     args.seq) ))["elapsed"]


    if (ParSib=="par") {
      PC[i,,] <- PedCompare(Ped1 = PedA, Ped2 = Seq.i$PedigreePar)$Counts["GG",,]
      if (quiet != "very")  cat("AR: ", round(sum(PC[i,"Match",])/sum(PC[i,"Total",]),3), "\t",
          "# ER:", sum(PC[i,c("Mismatch", "P2only"),]), "\n")
    } else {
      PC[i,,,] <- PedCompare(Ped1 = PedA, Ped2 = Seq.i$Pedigree)$Counts
      if (quiet != "very")  cat("AR: ", round(sum(PC[i,"TT","Match",])/sum(PC[i,"TT","Total",]),3), "\t",
          "# ER:", sum(PC[i,"TT",c("Mismatch", "P2only"),]), "\n")
    }
  }
  #~~~~~~~~~~~

  if (ParSib=="par") {
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

  RunParams <- list(SimGeno_default = formals(SimGeno),
                    sequoia_default = formals(sequoia),
                    EstConf_default = formals(EstConf),
                    EstConf_specified = match.call())

  if (return.PC) {
    list(SimCounts = PC, ConfProb = ECP,
         RunParams = RunParams, RunTime = RunTime)
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
