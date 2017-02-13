#' Convert genotype file
#'
#'  Convert a genotype file from PLINK's .ped or .RAW, or Colony's 2-column-
#'  per-marker format, to sequoia's 1-column-per-marker format.
#'
#'  The following formats are currently supported, specified by 'InFormat' and
#'   'OutFormat':
#'  \itemize{
#'    \item{ped: }{No header row, 6 descriptive columns, genotypes are coded as
#'    A, C, T, G, missing as 0, in 2 columns per marker. NOTE: not yet
#'    implented, use PLINK's --recodeA to convert this format to "raw".}
#'    \item{raw: }{Header row with SNP names, 6 descriptive columns,
#'     genotypes are coded as 0, 1, 2, missing as NA, in 1 column per marker.}
#'    \item{col: }{No header row, 1 descriptive column, genotypes are coded as
#'    numeric values, missing as 0, in 2 columns per marker.}
#'    \item{seq: }{No header row, 1 descriptive column genotypes are coded as
#'    0, 1, 2, missing as -9, in 1 column per marker.}
#'  }
#'
#' @param InFile character string with name of genotype file to be converted
#' @param InFormat One of "raw", "ped", "col" or "seq", see Details.
#' @param OutFile character string with name of converted file. If NA, return
#'   matrix with genotypes in console; if NULL, write to "GenoForSequoia.txt".
#' @param OutFormat as InFormat, default is "seq".
#' @param quiet suppress messages
#'
#' @return The converted genotype data is written to the specified file.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{sequoia}}
#'
#' @examples
#' \dontrun{
#' # Requires PLINK installed & in PATH:
#' system("cmd", input = "plink --file mydata --extract plink.prune.in
#'   --recodeA --out PlinkOUT")
#' #
#' GenoM <- GenoConvert(InFile = "PlinkOUT.raw", OutFile = NA)
#' }
#'
#' @export

GenoConvert <- function(InFile = NULL,
                        InFormat = "raw",
                        OutFile = NULL,
                        OutFormat = "seq",
                        quiet = FALSE) {
  if (OutFormat == "seq" & is.null(OutFile)) {
    OutFile <- "GenoForSequoia.txt"
  } else if (is.null(OutFile)) {
    stop("Please provide 'OutFile'")
  }

  if (interactive() & !quiet) {
    if (file.exists(OutFile)) {
      ANS <- readline(prompt = paste("WARNING: ", OutFile, " will be overwritten.",
                                     "Press <N> to abort, or any other key to continue."))
    } else {
      ANS <- readline(prompt = paste("Genotypes will be written to ", OutFile,
                                     " . Press <N> to abort, or any other key to continue."))
    }
    if (substr(ANS, 1, 1) %in% c("N", "n")) stop()
  }


  GenoTmp <- readLines(InFile)

  if (OutFormat == "seq") {  # TODO: shell:  sed -i 's/\bNA\b/-9/g' test.raw
    if (InFormat == "raw") {
      TmpL    <- strsplit(GenoTmp[-1], split = " ")  # skip row w marker names
      GenoOUT <- plyr::ldply(TmpL, function(x) x[-c(1, 3:6)])
      GenoTmp2 <- cbind(GenoOUT[, 1],
                        apply(GenoOUT[, -1], 2, function(x) gsub("NA", "-9", x)))
      if (!is.na(OutFile)) {
        utils::write.table(GenoTmp2, file = OutFile,
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      } else return(GenoTmp2)

    } else if (InFormat == "col") {
      TmpL    <- strsplit(GenoTmp, split = " ")
      GC <- plyr::ldply(TmpL)
      IDs_geno <- GC[, 1]
      GC <- as.matrix(GC[, -1])
      GCM <- matrix(NA, nrow(GC), ncol(GC))
      for (i in 1:ncol(GC)) {
        GCM[,i] <- as.numeric(as.factor(GC[, i]))-1
      }
      GCA <- array(dim=c(2, nrow(GC), ncol(GC)/2))
      GCA[1,,] <- GCM[, seq(1,ncol(GC)-1,2)]
      GCA[2,,] <- GCM[, seq(2,ncol(GC),2)]
      GS <- apply(GCA, 2:3, sum)
      GS[is.na(GS)] <- -9
      row.names(GS) <- IDs_geno
      if (!is.na(OutFile)) {
        utils::write.table(GS, OutFile, quote = FALSE, col.names = FALSE)
      } else return(GS)

    } else if (InFormat == "ped") {
      stop("not yet implemented")
    } else {
      stop("not implemented")
    }
  } else if (InFormat == "seq") {
      if (OutFormat == "col") {
        dc <- list("0" = c(1,1), "1" = c(1,2), "2" = c(2,2), "-9" = c(0,0))
        TmpL <- strsplit(GenoTmp, split=" ")
        Geno <- plyr::ldply(TmpL)
        IDs_geno <- Geno[,1]
        Geno <- as.matrix(Geno[, -1])
        GenoA <- array(dim=c(nrow(Geno), 2, ncol(Geno)))
        for (i in 1:nrow(Geno)) {
          GenoA[i,,] <- sapply(Geno[i,], function(z) dc[[z]])
        }
        GenoM <- matrix(GenoA, nrow(Geno))
        row.names(GenoM) <- IDs_geno
        if (!is.na(OutFile)) {
          utils::write.table(GenoM, OutFile, quote=FALSE, col.names=FALSE)
        } else return(GenoM)

      } else {
        stop("not yet implemented")
      }
  } else {
    stop("not implemented")
  }
}
