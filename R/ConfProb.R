#============================================================================
#============================================================================
#' @title Estimate confidence probability
#'
#' @description Estimate confidence and assignment error rate by repeatedly
#'   simulating genotype data from a reference pedigree using
#'   \code{\link{SimGeno}}, reconstruction a pedigree from this using
#'   \code{\link{sequoia}}, and counting the number of mismatches using
#'   \code{\link{PedCompare}}.
#'
#' @details The confidence probability is taken as the number of correct
#'   (matching) assignments, divided by all assignments made in the
#'   \emph{observed} (inferred-from-simulated) pedigree. In contrast, the false
#'   negative & false positive assignment rates are proportions of the number of
#'   parents in the \emph{true} (reference) pedigree. Each rate is calculated
#'   separatedly for dams & sires, and separately for each category
#'   (\strong{G}enotyped/\strong{D}ummy(fiable)/\strong{X} (none)) of
#'   individual, parent and co-parent.
#'
#'  This function does not know which individuals in \code{Pedigree} are
#'  genotyped, so the confidence probabilities need to be added to the Pedigree
#'  by the user as shown in the example at the bottom.
#'
#'  A confidence of `1' assignments on simulated data were correct for that
#'  category-combination. It should be interpreted as (and perhaps modified to)
#'  \eqn{> 1 - 1/N}, where sample size \code{N} is given in the last column of
#'  the \code{ConfProb} and \code{PedErrors} dataframes in the output. The same
#'  applies for a false negative/positive rate of `0'.
#'
#' @section Assumptions: Because the actual true pedigree is (typically)
#'   unknown, the provided reference pedigree is used as a stand-in and assumed
#'   to be the true pedigree, with unrelated founders. It is also assumed that
#'   the probability to be genotyped is equal for all parents; in each
#'   iteration, a new random set of parents (proportion set by \code{ParMis}) is
#'   mimicked to be non-genotyped. In addition, SNPs are assumed to segregate
#'   independently.
#'
#' @param Pedigree Reference pedigree from which to simulate, dataframe with
#'   columns id-dam-sire. Additional columns are ignored
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
#' @param quiet suppress messages. `very' also suppresses simulation counter,
#'   TRUE just runs SimGeno and sequoia quietly.
#'
#' @return a list, with the main results in dataframe \code{ConfProb} and array
#'   \code{PedErrors}. \code{ConfProb} has 7 columns:
#' \item{id.cat, dam.cat, sire.cat}{Category of the focal individual, dam, and
#' sire, in the pedigree inferred based on the simulated data. Coded as
#' G=genotyped, D=dummy, X=none}
#' \item{dam.conf}{Probability that the dam is correct, given the categories of
#' the assigned dam and sire (ignoring whether or not the sire is correct).
#' Rounded to \code{nchar(N)} significant digits}
#' \item{sire.conf}{as dam.conf, for the sire}
#' \item{pair.conf}{Probability that both dam and sire are correct, given their
#' categories}
#' \item{N}{Number of individuals per category-combination, across all
#' \code{nSim} simulations}
#'
#' array \code{PedErrors} has three dimensions:
#' \item{class}{\itemize{
#'   \item FalseNeg(atives): could have been assigned but was not
#' (individual + parent both genotyped or dummyfiable; P1only in
#' \code{PedCompare}).
#'   \item FalsePos(itives): no parent in reference pedigree, but
#' one was assigned based on the simulated data (P2only)
#'   \item Mismatch: different parents between the pedigrees}}
#' \item{cat}{Category of individual + parent, as a two-letter code where the first letter
#' indicates the focal individual and the second the parent; G=Genotyped, D=Dummy, T=Total}
#' \item{parent}{dam or sire}
#'
#' The other list elements are:
#'   \item{Pedigree.reference}{the pedigree from which data was simulated}
#'   \item{Pedigree.inferred}{a list with for each iteration the inferred
#'     pedigree based on the simulated data}
#'   \item{SimSNPd}{a list with for each iteration the IDs of the individuals
#'     simulated to have been genotyped}
#'   \item{RunParams}{a list with the current call to EstConf, as well as the
#'   default parameter values for \code{EstConf, SimGeno}, and \code{sequoia}.}
#'   \item{RunTime}{\code{sequoia} runtime per simulation in seconds, as
#'     measured by \code{\link{system.time}()['elapsed']}.}
#'
#' @seealso \code{\link{SimGeno}, \link{sequoia}, \link{PedCompare}}
#'
#' @examples
#' \dontrun{
#' data(Ped_HSg5, LH_HSg5, package="sequoia")
#'
#' ## Example A: parentage assignment only
#' conf.A <- EstConf(Pedigree = Ped_HSg5, LifeHistData = LH_HSg5,
#'    args.sim = list(nSnp = 100, SnpError = 5e-3, ParMis=c(0.2, 0.5)),
#'    args.seq = list(MaxSibIter = 0, Err=1e-3, Tassign=0.5),
#'    nSim = 2)
#'
#' # parent-pair confidence, per category:
#' conf.A$ConfProb
#'
#' # calculate (correct) assignment rates (ignores co-parent)
#' 1 - apply(conf.A$PedErrors, c(1,3), sum, na.rm=TRUE)
#'
#' ## Example B: with sibship clustering, based on sequoia inferred pedigree
#' RealGenotypes <- SimGeno(Ped = Ped_HSg5, nSnp = 100,
#'                          ParMis=c(0.19,0.53), SnpError = 6e-3)
#' SeqOUT <- sequoia(GenoM = RealGenotypes,
#'                   LifeHistData = LH_HSg5,
#'                   Err=5e-3, MaxSibIter=10)
#'
#' conf.B <- EstConf(Pedigree = SeqOUT$Pedigree,
#'               LifeHistData = LH_HSg5,
#'                args.sim = list(nSnp = 100, SnpError = 5e-3,
#'                                ParMis=c(0.2, 0.5)),
#'               args.seq = list(Err=5e-3, MaxSibIter = 10),
#'               nSim = 3)
#' Ped.withConf <- getAssignCat(Pedigree = SeqOUT$Pedigree,
#'                              Genotyped = rownames(RealGenotypes))
#' Ped.withConf <- merge(Ped.withConf, conf.B$ConfProb, all.x=TRUE)
#' Ped.withConf <- Ped.withConf[, c("id","dam","sire", "dam.conf", "sire.conf",
#'                                  "id.cat", "dam.cat", "sire.cat")]
#' }
#'
#' @export
#'

EstConf <- function(Pedigree = NULL,
                    LifeHistData = NULL,
                    args.sim = list(nSnp = 400, SnpError = 1e-3, ParMis=c(0.4, 0.4)),
                    args.seq = list(MaxSibIter = 10, Err=1e-3, Tassign=0.5),
                    nSim = 10,
                    quiet=TRUE)
{
  if (is.null(Pedigree))  stop("Please provide Pedigree")
  if (is.null(LifeHistData))  stop("Please provide LifeHistData")
  if (!is.null(args.sim) & !is.list(args.sim))  stop("args.sim should be a list or NULL")
  if (!is.null(args.seq) & !is.list(args.seq))  stop("args.seq should be a list or NULL")
  if (!is.wholenumber(nSim) || nSim<1)  stop("nSim must be a positive number")

  if ("Err" %in% names(args.sim)) {
    args.sim[["SnpError"]] <- args.sim[["Err"]]
    args.sim[["Err"]] <- NULL    # common confusion, otherwise fuzy matching with 'ErrorFM'.
  }

  Ped.ref <- Pedigree[,1:3]
  if (any(substr(unlist(Ped.ref),1,6) %in% c("sim_F0", "sim_M0"))) {
    stop("Please don't use 'sim_F' or 'sim_M' in reference pedigree")
  }

  if ("MaxSibIter" %in% names(args.seq)) {
    ParSib <- ifelse(args.seq$MaxSibIter > 0, "sib", "par")
  } else {
    ParSib <- "sib"   # default MaxSibIter=10
  }

  if (quiet != "very") {
    if (ParSib == "par") {
      message("MaxSibIter=0: Simulating parentage assignment only ...")
    } else {
      message("MaxSibIter>0: Simulating full pedigree reconstruction ...")
    }
  }
  utils::flush.console()
  seq.quiet <- ifelse(is.logical(quiet), quiet,
                    ifelse(quiet == "very", TRUE, FALSE))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate genotypes & infer pedigree
  RunTime <- rep(NA, nSim)
  Pedigree.inferred <- list()
  SimSNPd <- list()
  for (i in 1:nSim) {
    if (quiet != "very")  cat("i=", i, "\t", format(Sys.time(), "%H:%M:%S"), "\n")

    GM <- do.call(SimGeno, c(list(Pedigree=Pedigree), args.sim))
    SimSNPd[[i]] <- rownames(GM)
    RunTime[i] <- system.time(Seq.i <- do.call(sequoia, c(list(GenoM = GM,
                                                               LifeHistData = LifeHistData,
                                                               DummyPrefix = c("sim_F", "sim_M"),
                                                               quiet = seq.quiet,
                                                               CalcLLR = FALSE,
                                                               FindMaybeRel = FALSE,
                                                               Plot = FALSE),
                                                          args.seq) ))["elapsed"]
    if (ParSib == "par") {
      Pedigree.inferred[[i]] <- Seq.i$PedigreePar
    } else {
      Pedigree.inferred[[i]] <- Seq.i$Pedigree
    }
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # confidence probabilities
  nSimz <- ifelse(nSim>1, nSim,2)  # else problems w R auto-dropping dimension
  CatNames <- c("G", "D", "X")
  CPP <- array(dim=c(nSimz, 3,3,3,3),
               dimnames = list(1:nSimz, id.cat=CatNames,
                               dam.cat=CatNames, sire.cat=CatNames,
                               Conf=c("dam.conf", "sire.conf", "pair.conf")))
  Ni <- array(0, dim=c(nSimz, 3,3,3),
              dimnames = dimnames(CPP)[1:4])

  for (i in 1:nSim) {
    PC.rev <- PedCompare(Ped1 = Pedigree.inferred[[i]], Ped2 = Ped.ref,
                         SNPd = SimSNPd[[i]], Symmetrical=FALSE)
    CPP[i,c("G","D"),,,] <- CalcPairConf(PC.rev$Counts.detail, ParSib)
    Ni[i,,,] <- apply(PC.rev$Counts.detail, 1:3, sum)
  }
  if (ParSib == "par") {
    Ni[,,"X",] <- Ni[,,"X",] + Ni[,,"D",]
    Ni[,,"D",] <- 0
    Ni[,,,"X"] <- Ni[,,,"X"] + Ni[,,,"D"]
    Ni[,,,"D"] <- 0
  }

  # weighed mean across iterations (or not?)
  Conf.A <- array(dim=c(3,3,3,3), dimnames = c(dimnames(CPP)[2:5]))
  for (x in 1:3) {
    Conf.A[,,,x] <- apply(CPP[,,,,x] * Ni, 2:4, sum, na.rm=TRUE) / apply(Ni, 2:4, sum, na.rm=TRUE)
  }
  Conf.A[,"X",,c("dam.conf", "pair.conf")] <- NA
  Conf.A[,,"X",c("sire.conf", "pair.conf")] <- NA

  Conf.df <- ArrToDF(Conf.A, Ni, ParSib)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # assignment errors (ignores co-parent)
  PedErrors.r <- array(dim=c(nSim, 7,3,2),
                       dimnames = list(1:nSim,
                                       cat = c("GG", "GD", "GT", "DG", "DD", "DT", "TT"),
                                       class = c("FalseNeg", "FalsePos", "Mismatch"),
                                       parent = c("dam", "sire")))
  Nr <- array(dim=c(nSim, 7,2),
              dimnames = dimnames(PedErrors.r)[c(1,2,4)])
  for (i in 1:nSim) {
    PC.fwd <- PedCompare(Ped1 = Ped.ref, Ped2 = Pedigree.inferred[[i]],
                         SNPd = SimSNPd[[i]], Symmetrical=FALSE)
    PedErrors.r[i,,,] <- sweep(PC.fwd$Counts[,c("P1only", "P2only","Mismatch"),], c(1,3),
                               PC.fwd$Counts[,"Total",], "/")
    Nr[i,,] <- PC.fwd$Counts[,"Total",]
  }

  # average across iterations:
  PedErrors <- apply(PedErrors.r, 2:4, mean, na.rm=TRUE)
  Ntot <- apply(Nr, 2:3, sum)
  for (x in 1:3) {
    PedErrors[,x,] <- signif(PedErrors[,x,], digits=nchar(Ntot))
  }
  PedErrors[c("GG", "GD", "DG", "DD"), "FalsePos", ] <- NA


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # out
  RunParams <- list(SimGeno_default = formals(SimGeno),
                    sequoia_default = formals(sequoia),
                    EstConf_default = formals(EstConf),
                    EstConf_specified = match.call())

  return( list(ConfProb = Conf.df,
               PedErrors = PedErrors,
               Pedigree.reference = Ped.ref,
               Pedigree.inferred = Pedigree.inferred,
               SimSNPd = SimSNPd,
               RunParams = RunParams,
               RunTime = RunTime) )
}



#============================================================================
CalcPairConf <- function(PedCompDetails, ParSib="sib") {   # per-iteration
  CP <- array(NA, dim=c(2,3,3,3),
              dimnames=c(list("id.cat"=c("G","D")), dimnames(PedCompDetails)[2:3],
                               list(Conf=c("dam", "sire", "pair"))))
  tots1 <- c("Match", "Mismatch", "P2only")   # total assigned in pedigree 2
  not1 <- c("P1only", "_")
  m <- ifelse(ParSib=="par", 1, 2)
  CD <- PedCompDetails
  for (i in 1:m) {
    # parent-pairs
    for (j in 1:m) {
      for (k in 1:m) {
        CP[i,j,k, "dam"] <- sum(CD[i,j,k,"Match", tots1]) / sum(CD[i,j,k,tots1, tots1])
        CP[i,j,k, "sire"] <- sum(CD[i,j,k,tots1,"Match"]) / sum(CD[i,j,k,tots1,tots1])
        CP[i,j,k, "pair"] <- CD[i,j,k,"Match", "Match"] / sum(CD[i,j,k,tots1, tots1])
      }
    }
    # single dams
    for (j in 1:m) {
      CP[i,j,"X","dam"] <- sum(CD[i,j,,"Match", not1]) / sum(CD[i,j,, tots1, not1])
    }
    # single sires
    for (k in 1:m) {
      CP[i,"X",k,"sire"] <- sum(CD[i,,k,not1,"Match"]) / sum(CD[i,,k, not1,tots1])
    }
  }
  return( CP )
}



#============================================================================

ArrToDF <- function(A.P, A.N, PS) {
  if (PS == "par") {
    DF <- plyr::adply(A.P["G",c("G","X"),c("G","X"),,drop="FALSE"], 1:3)
  } else {
    DF <- plyr::adply(A.P[c("G","D"),,,], 1:3)
  }
  N.df <- plyr::adply(A.N, 2:4, sum)
  names(N.df)[4] <- "N"
  DF <- merge(DF, N.df, all.x=TRUE)
  for (x in 1:3) {
    DF[,x] <- factor(DF[,x], levels=c("G", "D", "X"))
  }
  for (x in 4:6) {
    DF[,x] <- signif(DF[,x], digits = nchar(DF$N))
  }
  DF <- DF[order(DF[,1], DF[,2], DF[,3]), ]
  rownames(DF) <- 1:nrow(DF)
  return( DF )
}

#============================================================================
