#' @title SNP summary statistics
#'
#' @description Estimate allele frequency (AF), missingness and Mendelian
#' errors per SNP.
#'
#' @details Calculation of these summary statistics can be done in PLINK, and
#' SNPs with low minor allele freuqency or high missigness should be filtered
#' out using PLINK prior to pedigree reconstruction. This function is merely
#' provided as an aid to inspect the relationship between AF, missingness
#' and error to find a suitable combination of thresholds to use.
#'
#' The error count includes both the number of parent-offspring pairs that are
#' opposing homozygotes (parent is AA and offspring is aa), as Mendelian errors
#' in parent-parent-offspring trios (e.g. parents AA and aa, but offspring not
#' Aa).
#'
#' The underlying genotyping error can not be easily estimated from the number
#' of Mendelian errors, as many errors may go undetected and a single error in
#'  a prolific individual can result in a high number of Mendelian errors.
#'  Moreover, a high error rate may interfere with pedigree reconstruction, and
#'  succesful assignment will be biased towards parents with lower error count.
#'
#' @param GenoM  Genotype matrix, in sequoia's format: 1 column per SNP, 1 row
#' per individual, genotypes coded as 0/1/2/-9, and rownames giving individual
#' IDs.
#' @param Ped  a dataframe with 3 columns: ID - parent1 - parent2. Additional
#'  columns and non-genotyped individuals are ignored. Only used to estimate
#'  the error rate.
#' @param Plot  show histograms of the results?
#'
#' @return a matrix with a number of rows equal to the number of SNPs
#'  (=number of columns of GenoM) and 2 or 3 columns:
#' \item{AF}{Allele frequency of the 'second allele' (the one for which the
#'   homozygote is coded 2)}
#' \item{Mis}{Proportion of missing calls}
#' \item{ER}{(only when Ped provided) number of Mendelian errors in parent-
#'  offspring pairs (i.e. the number of opposing homozygotes, 'OHdam' &
#'   'OHsire' in pedigree) and parent-parent-offspring trios ('MEpairs' in
#'   pedigree).}
#'
#' @seealso  \code{\link{GenoConvert}}
#'
#' @export

SnpStats <- function(GenoM,
                     Ped = NULL,
                     Plot = TRUE)
{
  Mis <- apply(GenoM, 2, function(x) sum(x==-9))/nrow(GenoM)
  AF <- apply(GenoM, 2, function(x) sum(x[x!=-9])/(2*sum(x!=-9)))

  if (!is.null(Ped)) {
    Par <- Ped[,1:3]
    for (i in 1:3) {
      Par[,i] <- as.character(Par[,i])
      Par[!Par[,i] %in% rownames(GenoM), i] <- NA
    }
    Par <- Par[!is.na(Par[,1]), ]
    ER <- OHP(GenoM, Par)
    OUT <- cbind(AF, Mis, ER)

  } else {
    OUT <- cbind(AF, Mis)
  }

  if (Plot) {
    if (is.null(Ped)) {
      op <- par(mfrow=c(1,3), mai=c(.8,.8,.5,.1))
    } else {
      op <- par(mfrow=c(2,2), mai=c(.8,.8,.5,.1))
    }
    hist(OUT[,"AF"], breaks=ncol(GenoM)/5, col="grey", main="Frequency '1' allele",
         xlab="", cex.main=1.2)
    hist(OUT[,"Mis"], breaks=ncol(GenoM)/5, col="grey", main="Missingness",
         xlab="", cex.main=1.2)
    MAF <- ifelse(OUT[,"AF"] <= 0.5, OUT[,"AF"], 1-OUT[,"AF"])
    if (!is.null(Ped)) {
      hist(OUT[,"ER"], breaks=c(0:(max(OUT[,"ER"])+2))-.5, col="grey", main="Mendelian Errors",
         xlab="", cex.main=1.2)
    }
    plot(MAF, OUT[,"Mis"], pch=16, cex=1.2, xlim=c(0,0.5), xlab="Minor Allele Frequency",
         ylab="Missingness", cex.lab=1.3)
    if (!is.null(Ped)) {
      if (any(OUT[,"ER"] > 0)) {
        q95 <- OUT[,"ER"] > stats::quantile(OUT[,"ER"], prob=0.95)
        points(MAF[q95], OUT[q95, "Mis"], pch=16, col="red")
        legend("topleft", "5% highest ER", pch=16, col="red", inset=.01)
      }
    }
    par(op)
  }

  rownames(OUT) <- paste0("SNP", formatC(1:nrow(OUT), width=ifelse(nrow(OUT)<1000, 3, 4), flag="0"))
  return( OUT )
}


#===========

OHP <- function(GenoM, Par) {
  GenoMx <- GenoM
  GenoMx[GenoMx==-9] <- 3
  GenoMx <- rbind(GenoMx, "NA" = 3)
  GenoMx <- GenoMx +1
  Par$RowI <- sapply(Par$id, function(x, y) which(y == x), y = rownames(GenoMx))
  Par$RowD <- with(Par, sapply(dam, function(x, y) ifelse(is.na(x), nrow(GenoMx),
                                               which(y == x)), y = rownames(GenoMx)))
  Par$RowS <- with(Par, sapply(sire, function(x, y) ifelse(is.na(x), nrow(GenoMx),
                                            which(y == x)), y = rownames(GenoMx)))

  # mendelian errors
  MER <- array(0, dim=c(4,4,4))  # offspr - mother - father
  MER[1:3,,1] <- matrix(c(0,1,2, 0,0,1, 1,0,1, 0,0,1), 3,4)  # 0/1/2/NA
  MER[1:3,,2] <- matrix(c(0,0,1, 0,0,0, 1,0,0, 0,0,0), 3,4)
  MER[1:3,,3] <- matrix(c(1,0,1, 1,0,0, 2,1,0, 1,0,0), 3,4)
  MER[1:3,,4] <- matrix(c(0,0,1, 0,0,0, 1,0,0, 0,0,0), 3,4)

  ERcount <- matrix(0, nrow(Par), ncol(GenoMx))
  for (i in 1:nrow(Par)) {
    ERcount[i, ] <- MER[cbind(GenoMx[Par$RowI[i], ],
                          GenoMx[Par$RowD[i], ],
                          GenoMx[Par$RowS[i], ])]
  }
  apply(ERcount, 2, sum, na.rm=T)
}
