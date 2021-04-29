#' @title SNP Summary Statistics
#'
#' @description Estimate allele frequency (AF), missingness and Mendelian
#' errors per SNP.
#'
#' @details Calculation of these summary statistics can be done in PLINK, and
#'   SNPs with low minor allele frequency or high missingness should be filtered
#'   out prior to pedigree reconstruction. This function is provided as an aid
#'   to inspect the relationship between AF, missingness and genotyping error to
#'   find a suitable combination of SNP filtering thresholds to use.
#'
#'   For pedigree reconstruction, SNPs with zero or one copies of the alternate
#'   allele in the dataset (MAF \eqn{\le 1/2N}) are considered fixed, and
#'   excluded.
#'
#' @section Estimated genotyping error:
#' The error rate is estimated from the number of opposing homozygous cases (OH,
#' parent is AA and offspring is aa) Mendelian errors (ME, e.g. parents AA and
#' aa, but offspring not Aa) in parent-parent-offspring trios, and OH cases for
#' offspring with a single genotyped parent.
#'
#' The estimated error rates will not be as accurate as from duplicate samples.
#' A single error in an individual with many offspring will be counted as many
#' times, potentially resulting in non-sensical values of 'Err.hat' close to 1.
#' On the other hand, errors in individuals without parents or offspring will
#' not be counted at all. Moreover, a high error rate may interfere with
#' pedigree reconstruction, and successful assignment will be biased towards
#' parents with lower error count. Nonetheless, it may provide a ballpark
#' estimate for the average error rate, which can be useful for subsequent
#' (rerun of) pedigree reconstruction.
#'
#' @param GenoM  genotype matrix, in sequoia's format: 1 column per SNP, 1 row
#'   per individual, genotypes coded as 0/1/2/-9, and rownames giving individual
#'   IDs.
#' @param Pedigree  dataframe with 3 columns: ID - parent1 - parent2.
#'   Additional columns and non-genotyped individuals are ignored. Used to
#'   estimate the error rate.
#' @param ErrFlavour function that takes the genotyping error rate \code{Err} as
#'   input, and returns a 3x3 matrix of observed (columns) conditional on actual
#'   (rows) genotypes, or choose from inbuilt ones as used in sequoia
#'   'version2.0', 'version1.3', or 'version1.1'. See \code{\link{ErrToM}}.
#' @param Plot  show histograms of the results?
#'
#' @return A matrix with a number of rows equal to the number of SNPs
#'  (=number of columns of GenoM), and when no Pedigree is provided 2 columns:
#' \item{AF}{Allele frequency of the 'second allele' (the one for which the
#'   homozygote is coded 2)}
#' \item{Mis}{Proportion of missing calls}
#' When a Pedigree is provided, there are 7 additional columns:
#' \item{n.dam, n.sire, n.pair}{Number of dams, sires, parent-pairs succesfully
#'   genotyped for the SNP}
#' \item{OHdam, OHsire}{Count of number of opposing homozygous cases}
#' \item{MEpair}{Count of Mendelian errors, includes opposing homozygous cases}
#' \item{Err.hat}{Error rate, as estimated from the joined offspring-parent
#'   (-parent) genotypes and the presumed error structure (\code{ErrFlavour})}
#'
#' @seealso  \code{\link{GenoConvert}} to convert from various data formats;
#'   \code{\link{CheckGeno}} to check the data is in valid format for sequoia
#'   and exclude monomorphic SNPs etc., \code{\link{CalcOHLLR}} to calculate OH
#'   & ME per individual.
#'
#' @examples
#' \donttest{
#' data(Ped_HSg5)
#' Genotypes <- SimGeno(Ped_HSg5, nSnp=400, CallRate = runif(400, 0.2, 0.8),
#'   SnpError = 0.05)
#' SNPstats <- SnpStats(Genotypes, Pedigree=Ped_HSg5)
#' }
#'
#' @importFrom graphics plot points legend
#'
#' @export

SnpStats <- function(GenoM,
                     Pedigree = NULL,
                     ErrFlavour = "version2.0",
                     Plot = TRUE)
{
  Mis <- apply(GenoM, 2, function(x) sum(x==-9))/nrow(GenoM)
  AF <- apply(GenoM, 2, function(x) sum(x[x!=-9])/(2*sum(x!=-9)))

  if (is.logical(Pedigree)) {
    Plot <- Pedigree
    Pedigree <- NULL
  }
  if (!is.null(Pedigree)) {
    Par <- Pedigree[,1:3]
    for (i in 1:3) {
      Par[,i] <- as.character(Par[,i])
      Par[!Par[,i] %in% rownames(GenoM), i] <- NA
    }
    Par <- Par[!is.na(Par[,1]), ]
    ER <- EstErr(GenoM, Par, ErrFlavour)
    OUT <- cbind(AF, Mis, ER)

  } else {
    OUT <- cbind(AF, Mis)
  }

  if (Plot) {
    oldpar <- par(no.readonly = TRUE)
    oldpar <- oldpar[!names(oldpar) %in% c("pin", "fig")]
    if (is.null(Pedigree)) {
      par(mfrow=c(1,3), mai=c(.9,.8,.2,.1))
    } else {
      par(mfrow=c(2,3), mai=c(.9,.8,.2,.1))
    }
    OK <- tryPlot(hist,
                  OUT[,"AF"], breaks=ncol(GenoM)/5, col="grey", main="",
                  xlab="Frequency '1' allele", cex.lab=1.3, ylab="",
#                  ErrMsg = "SnpStats: Plotting area too small",
                  oldpar = oldpar)
    if (OK) {  # if first fits, rest fits, and v.v.
      hist(OUT[,"Mis"], breaks=ncol(GenoM)/5, col="grey", main="",
           xlab="Missingness", cex.lab=1.3, ylab="")
      MAF <- ifelse(OUT[,"AF"] <= 0.5, OUT[,"AF"], 1-OUT[,"AF"])
      if (!is.null(Pedigree)) {
        hist(OUT[,"Err.hat"], breaks=ncol(GenoM)/5, col="grey",
             main="", xlab="Error rate", cex.lab=1.3, ylab="")
      }
      plot(MAF, OUT[,"Mis"], pch=16, cex=1.2, xlim=c(0,0.5),
           xlab="Minor Allele Frequency", ylab="Missingness", cex.lab=1.3)
      if (!is.null(Pedigree)) {
        q95.e <- OUT[,"Err.hat"] > stats::quantile(OUT[,"Err.hat"], prob=0.95)
        points(MAF[q95.e], OUT[q95.e, "Mis"], pch=16, col="red")
        legend("topleft", "5% highest error", pch=16, col="red", inset=.01)

        plot(MAF, OUT[,"Err.hat"], pch=16, cex=1.2, xlim=c(0,0.5),
           xlab="Minor Allele Frequency", ylab="Error rate", cex.lab=1.3)
        q95.m <- OUT[,"Mis"] > stats::quantile(OUT[,"Mis"], prob=0.95)
        points(MAF[q95.m], OUT[q95.m, "Err.hat"], pch=16, col="red")
        legend("topleft", "5% highest missingness", pch=16, col="red", inset=.01)

        plot(OUT[,"Mis"], OUT[,"Err.hat"], pch=16, cex=1.2,
           xlab="Missingness", ylab="Error rate", cex.lab=1.3)
        q95.maf <- MAF < stats::quantile(MAF, prob=0.05)
        points(OUT[q95.maf,"Mis"], OUT[q95.maf,"Err.hat"], pch=16, col="grey")
        legend("topright", "5% lowest MAF", pch=16, col="grey", inset=.01)
      }
      par(oldpar)  # restore old par settings
    }
  }

  rownames(OUT) <- paste0("SNP", formatC(1:nrow(OUT),
                                         width=ifelse(nrow(OUT)<1000, 3, 4),
                                         flag="0"))
  invisible( OUT )
}



#==============================================================================
#' @title Estimate Genotyping Error Rate
#'
#' @description Estimate genotyping error rate from Mendelian errors per SNP.
#'
#' @param GenoM  genotype matrix, in sequoia's format: 1 column per SNP, 1 row
#'   per individual, genotypes coded as 0/1/2/-9, and rownames giving individual
#'   IDs.
#' @param Par  pedigree dataframe, only genotyped parents are used.
#' @param ErrFlavour function that takes the genotyping error rate \code{Err} as
#'   input, and returns a 3x3 matrix of observed (columns) conditional on actual
#'   (rows) genotypes, or choose from inbuilt ones as used in sequoia
#'   'version2.0', 'version1.3', or 'version1.1'. See \code{\link{ErrToM}}.
#'
#' @return A dataframe with columns:
#' \item{Err.hat}{Error rate, as estimated from the joined offspring-parent
#'   (-parent) genotypes and the presumed error structure (\code{ErrFlavour})}
#' \item{n.dam, n.sire, n.pair}{Number of dams, sires, parent-pairs succesfully
#'   genotyped for the SNP}
#' \item{OHdam, OHsire}{Count of number of opposing homozygous cases}
#' \item{MEpair}{Count of Mendelian errors, includes opposing homozygous cases}
#'
#' @seealso \code{\link{SnpStats}}.
#'
#' @keywords internal

EstErr <- function(GenoM, Par, ErrFlavour = "version2.0") {
	GenoMx <- GenoM
  GenoMx[GenoMx==-9] <- 3
  GenoMx <- rbind(GenoMx, "NA" = 3)
  GenoMx <- GenoMx +1
  Par$RowI <- sapply(Par$id, function(x, y) which(y == x), y = rownames(GenoMx))
  Par$RowD <- with(Par, sapply(dam, function(x, y) ifelse(is.na(x), nrow(GenoMx),
                                               which(y == x)), y = rownames(GenoMx)))
  Par$RowS <- with(Par, sapply(sire, function(x, y) ifelse(is.na(x), nrow(GenoMx),
                                            which(y == x)), y = rownames(GenoMx)))

	Obs.OO.all <- array(dim=c(nrow(Par), 3, ncol(GenoMx)))  # 4 = NA
  for (i in 1:nrow(Par)) {
    Obs.OO.all[i,,] <- rbind(GenoMx[Par$RowI[i], ],
                             GenoMx[Par$RowD[i], ],
                             GenoMx[Par$RowS[i], ])
  }
  OO.trio <- plyr::aaply(Obs.OO.all, 3, function(M) table(factor(M[,1], levels=1:4),
                                                    factor(M[,2], levels=1:4),
                                                    factor(M[,3], levels=1:4)))

	AF <- apply(GenoM, 2, function(x) sum(x[x!=-9])/(2*sum(x!=-9)))

	ErrF <- ErrToM(flavour = ErrFlavour, Return = "function")

	E.hat <- numeric(ncol(GenoM))  # based on trio where possible
	for (l in 1:ncol(GenoM)) {
		E.hat[l] <- stats::optimise(CalcChi2, interval=c(0,1), q=AF[l],
		                     A.obs=OO.trio[l,,,], ErrF=ErrF)$minimum
	}

	# mendelian errors
  MER <- array(0, dim=c(4,4,4))  # offspr - mother - father
  MER[1:3,,1] <- matrix(c(0,1,2, 0,0,1, 1,0,1, 0,0,1), 3,4)  # 0/1/2/NA
  MER[1:3,,2] <- matrix(c(0,0,1, 0,0,0, 1,0,0, 0,0,0), 3,4)
  MER[1:3,,3] <- matrix(c(1,0,1, 1,0,0, 2,1,0, 1,0,0), 3,4)
  MER[1:3,,4] <- matrix(c(0,0,1, 0,0,0, 1,0,0, 0,0,0), 3,4)

  tmpdam <- apply(OO.trio, c(1,2,3), sum)  # sum over sire genotypes
  tmpsire <- apply(OO.trio, c(1,2,4), sum)
  Counts <- cbind(n.dam = apply(tmpdam, 1, function(M) sum(M[1:3, 1:3])),
                  OHdam = apply(tmpdam, 1, function(M) M[1,3] + M[3,1]),
                  n.sire = apply(tmpsire, 1, function(M) sum(M[1:3, 1:3])),
                  OHsire = apply(tmpsire, 1, function(M) M[1,3] + M[3,1]),
                  n.pair = apply(OO.trio, 1, function(A) sum(A[1:3, 1:3, 1:3])),
                  MEpair = sapply(1:ncol(GenoM), function(l) sum(MER * OO.trio[l,,,])))

	return( data.frame(Err.hat = E.hat, Counts) )
}



#==============================================================================
#' @title Chi-square Test on Observed vs Expected Genotypes
#'
#' @description For one SNP and all offspring-parent-parent trios or single
#'   parent-offspring pairs, calculate the expected genotype frequencies given
#'   the allele frequency, genotyping error rate, and error flavour, and perform
#'   a chi-square test.
#'
#' @param E  presumed genotyping error rate.
#' @param q  allele frequency.
#' @param A.obs  array of dim 4x4x4 with counts of joined
#'   offspring-parent-parent at the SNPs
#' @param ErrF  ErrFlavour; function that takes the genotyping error rate
#'   \code{Err} as input, and returns a 3x3 matrix of observed (columns)
#'   conditional on actual (rows) genotypes, or choose from inbuilt ones as used
#'   in sequoia 'version2.0', 'version1.3', or 'version1.1'. See
#'   \code{\link{ErrToM}}.
#'
#' @return The chisquare value of the test.
#'
#' @seealso \code{\link{EstErr}, \link{SnpStats}}.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' E.hat <- numeric(ncol(GenoM))  # based on trio where possible
#' for (l in 1:ncol(GenoM)) {
#' 		E.hat[l] <- stats::optimise(CalcChi2, interval=c(0,1), q=AF[l],
#' 		                     A.obs=OO.trio[l,,,], ErrF=ErrF)$minimum
#' } }

CalcChi2 <- function(E, q, A.obs, ErrF) {
  A.obs[4,,] <- NA   # id not SNPd: no information
  A.obs[,4,4] <- NA  # neither parent SNPd: no information
  mis <- c(dam = sum(A.obs[,4,], na.rm=T),
           sire = sum(A.obs[,,4], na.rm=T)) / sum(A.obs, na.rm=T)

  # TODO?: calc q from M.obs & adjust iteratively for E
  AHWE <- c((1-q)^2, 2*q*(1-q), q^2)

  # offspring - parent pair
  AKAP <- matrix(c(1-q, (1-q)/2, 0,
                       q,   1/2,     1-q,
                       0,   q/2,     q), nrow=3, byrow=TRUE)
  P.AA <- sweep(AKAP, 2, AHWE, "*")  # joined prob actual genotypes

  ErrM <- ErrF(E)
  P.OO <- t(ErrM) %*% P.AA %*% ErrM     # joined prob obs genotypes

  # offspring - dam - sire trio
  AKA2P <- array(dim=c(3,3,3))  # offspr - mother - father
  AKA2P[1,,] <- matrix(c(1,.5,0, .5,.25, 0,  0, 0,0), nrow=3)
  AKA2P[2,,] <- matrix(c(0,.5,1, .5, .5,.5,  1,.5,0), nrow=3)
  AKA2P[3,,] <- matrix(c(0, 0,0,  0,.25,.5,  0,.5,1), nrow=3)

  P.AAA <- sweep( sweep(AKA2P, 2, AHWE, "*"), 3, AHWE, "*")
  P.OOA <- array(dim=c(3,3,3))
  for (i in 1:3) {
    P.OOA[,,i] <- t(ErrM) %*% P.AAA[,,i]  %*% ErrM
  }

  P.OOO <- array(NA, dim=c(4,4,4))
  for (j in 1:3) {
    P.OOO[1:3,j,1:3] <- P.OOA[,j,] %*% ErrM
  }
  P.OOO <- P.OOO * (1 - mis["dam"]) * (1 - mis["sire"])
  P.OOO[1:3, 4, 1:3] <- mis["dam"] * P.OO
  P.OOO[1:3, 1:3, 4] <- mis["sire"] * P.OO
  P.OOO <- P.OOO / sum(P.OOO, na.rm=T)  # rounding errors

  A.exp <- P.OOO * sum(A.obs, na.rm=T)

  Chi2 <- sum((A.obs - A.exp)^2 / A.exp, na.rm=TRUE)
  return( Chi2 )
}
