#' @title Estimate genotyping error rate
#'
#' @description Estimate the genotyping error rates in SNP data, based on a
#'   pedigree and/or duplicates. Estimates probabilities (observed given
#'   actual) hom|other hom, het|hom, and hom|het. THESE ARE APPROXIMATE VALUES!
#'
#' @param GenoM  Genotype matrix
#' @param Pedigree  data.frame with columns id - dam - sire
#' @param Duplicates  matrix or data.frame with 2 columns, id1 & id2
#' @param Er_start  vector of length 3 with starting values for \code{optim}.
#' @param perSNP  logical, estimate error rate per SNP. WARNING not very
#'   precise, use only as an approximate indicator! Try on simulated data first,
#'   e.g. with \code{\link{SimGeno}}.
#'
#' @return  vector of length 3 with estimated genotyping error rates: the
#'  probabilities that
#' \itemize{
#'  \item hom|hom: an actual homozygote is observed as the other homozygote
#'  \item het|hom: an actual homozygote is observed as heterozygote
#'  \item hom|het: an actual heterozygote is observed as homozygote
#'  }
#'
#'  These are three independent parameters, that define the genotyping error
#'  matrix (see \code{\link{ErrToM}}) as follows:
#'
#' \tabular{lccc}{
#'             \tab \strong{0}       \tab \strong{1}   \tab \strong{2}  \cr
#'  \strong{0} \tab \eqn{1-E_1-E_2} \tab \eqn{E_2}    \tab \eqn{E_1} \cr
#'  \strong{1} \tab \eqn{E_3}        \tab \eqn{1-2E_3} \tab \eqn{E_3}   \cr
#'  \strong{2} \tab \eqn{E_1}       \tab \eqn{E_2}    \tab \eqn{1-E_1-E_2} \cr
#' }
#'
#'  Note that for \code{optim} a lower bound of 1e-6 and upper bound of 0.499
#'  are used; if these values are returned this should be interpreted as
#'  'inestimably small' and 'inestimably large', respectively. PLEASE DO NOT USE
#'  THESE VALUES AS INPUT IN SUBSEQUENT ANALYSIS BUT SUBSITUTE BY A SENSIBLE
#'  VALUE!!
#'
#' @details The result should be interpreted as approximate, ballpark estimates!
#' The estimated error rates from a pedigree will not be as accurate as from
#' duplicate samples. Errors in individuals without parents or offspring will
#' not be counted, and errors in individuals with only few offspring may not be
#' noted either. Deviation of genotype frequencies among founders from
#' Hardy-Weinberg equilibrium may wrongly be attributed to genotyping errors.
#' Last but not least, any pedigree errors will result in higher estimated
#' genotyping errors.
#'
#'
#' @importFrom stats optim setNames na.exclude
#'
#' @examples
#' GenoX <- SimGeno(Ped_griffin, nSnp=400, SnpError=c(0.01,0.07, 0.1),
#'                 ParMis=0.1, CallRate=0.9)
#' EstEr(GenoM=GenoX, Pedigree=Ped_griffin)
#'
#' @export

EstEr <- function(GenoM, Pedigree, Duplicates=NULL, Er_start=c(.05,.05,.05),
                  perSNP=FALSE)
{

  # check genotype data
  GenoM <- sequoia::CheckGeno(GenoM, quiet=TRUE, Plot=FALSE)
  Ng <- nrow(GenoM)
  gID <- rownames(GenoM)

  # check pedigree
  Ped <- sequoia::PedPolish(Pedigree, gID, DropNonSNPd=FALSE, NullOK=FALSE)

  # turn IDs into GenoM rownumbers
  GenoNums <- setNames(seq_along(gID), gID)
  NumPed <- matrix(0, nrow(Ped), 3,
                   dimnames=list(Ped$id, c('id','dam','sire')))
  for (x in 1:3) {
    NumPed[, x] <- GenoNums[Ped[, x]]
  }
  Parents <- NumPed[match(gID, Ped$id), c('dam','sire')]
  Parents[is.na(Parents)] <- 0

  # duplicates to vector with GenoM rownumbers
  DupsV <- rep(0, Ng)
  if (!is.null(Duplicates)) {
    NumDups <- cbind(id1 = GenoNums[Duplicates[,1]],
                     id2 = GenoNums[Duplicates[,2]])
    DupsV[NumDups[,'id1']] <- NumDups[,'id2']
    DupsV[NumDups[,'id2']] <- NumDups[,'id1']
  }

  optim_out <- optim(Er_start, CalcLL,
                     GenoM = GenoM, Parents = Parents, DupsV = DupsV,
                     method = 'L-BFGS-B', lower = rep(1e-6,3), upper = rep(0.499,3))

  Er_hat <- optim_out$par
  names(Er_hat) <- c('hom|hom', 'het|hom', 'hom|het')

  if (!perSNP) {

    return(Er_hat)

  } else {

    Err_per_SNP <- ErrPerSNP(Er_hat, GenoM, Parents, DupsV)
    return( Err_per_SNP )

  }

}



#==============================================================================

#' @title wrapper for Fortran function to calculate total likelihood
#'
#' @description function to be optimised
#'
#' @param Er_hat  length 3 vector with genotyping error rates
#' @param GenoM  Genotype matrix
#' @param Parents  Pedigree, already converted to rownumbers in GenoM
#' @param DupsV  vector with duplicate samples, already converted to rownumbers
#'
#' @return  the negative of the total log10-likelihood
#'
#' @useDynLib sequoia, .registration = TRUE
#'
#' @keywords internal

CalcLL <- function(Er_hat, GenoM, Parents, DupsV)
{
  nSnp <- ncol(GenoM)

  TMP <- .Fortran(ester,
                  ng = as.integer(nrow(GenoM)),
                  nl = as.integer(nSnp),
                  genov = as.integer(GenoM),
                  parentsv = as.integer(Parents),
                  dups = as.integer(DupsV),
                  errin = as.double(Er_hat),  # length 3, IN
                  totll = as.double(0),
                  cntobsact = as.double(rep(0,3*3*nSnp)))

  return( -TMP$totll)
}





#==============================================================================

#' @title wrapper for Fortran function to estimate error rate for each SNP
#'
#' @description WARNING: not very precise, especially not for low error rates.
#'
#' @param Er_hat  length 3 vector with genotyping error rates
#' @param GenoM  Genotype matrix
#' @param Parents  Pedigree, already converted to rownumbers in GenoM
#' @param DupsV  vector with duplicate samples, already converted to rownumbers
#'
#' @return  a matrix with 3 columns: for each SNP the probabilities (observed given
#'   actual) hom|other hom, het|hom, and hom|het.
#'
#' @useDynLib sequoia, .registration = TRUE
#'
#' @keywords internal
#'
#'
ErrPerSNP <- function(Er_hat, GenoM, Parents, DupsV)
{
  nSnp <- ncol(GenoM)

  TMP <- .Fortran(ester,
                  ng = as.integer(nrow(GenoM)),
                  nl = as.integer(nSnp),
                  genov = as.integer(GenoM),
                  parentsv = as.integer(Parents),
                  dups = as.integer(DupsV),
                  errin = as.double(Er_hat),  # length 3, IN
                  totll = as.double(0),
                  cntobsact = as.double(rep(0,3*3*nSnp)))

  # proportion actual (estimated) vs observed, per SNP
  PropsOA <- array(TMP$cntobsact, dim=c(3,3,ncol(GenoM)),
               dimnames = list('actual'=0:2, 'observed'=0:2, 'SNP'=1:nSnp))

  # divide by rowSums (=conditional on actual)
  # TODO: double check; divide by rowSums or not, for best correlation w actual error rate?
  POA <- plyr::aaply(PropsOA, 3, function(M) sweep(M, 1, rowSums(M), '/'))

  Err_per_SNP <- cbind('hom|hom' = sapply(1:nSnp, function(i)  (POA[i,'0','2'] + POA[i,'2','0'])/2),
                       'het|hom' = sapply(1:nSnp, function(i)  (POA[i,'0','1'] + POA[i,'2','1'])/2),
                       'hom|het' = sapply(1:nSnp, function(i)  (POA[i,'1','0'] + POA[i,'1','2'])/2))

  return(Err_per_SNP)
}


