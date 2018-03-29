#' @title Convert genotype file
#'
#' @description Convert a genotype file from PLINK's .raw, or Colony's
#' 2-column-per-marker format, to sequoia's 1-column-per-marker format.
#'
#' @details The following formats can be specified by 'InFormat' and
#' 'OutFormat':
#'  \itemize{
#'   \item{col: }{No header row, 1 descriptive column, genotypes are coded as
#'    numeric values, missing as 0, in 2 columns per marker.}
#'    \item{ped: }{No header row, 6 descriptive columns, genotypes are coded as
#'    A, C, T, G, missing as 0, in 2 columns per marker. NOTE: not yet
#'    implented, use PLINK's --recodeA to convert this format to "raw".}
#'    \item{raw: }{Header row with SNP names, 6 descriptive columns,
#'     genotypes are coded as 0, 1, 2, missing as NA, in 1 column per marker.}
#'    \item{seq: }{No header row, 1 descriptive column genotypes are coded as
#'    0, 1, 2, missing as -9, in 1 column per marker.}
#'  }
#'
#' @param InFile character string with name of genotype file to be converted
#' @param InFormat One of "raw", "col" or "seq", see Details.
#' @param OutFile character string with name of converted file. If NA, return
#'   matrix with genotypes in console; if NULL, write to "GenoForSequoia.txt".
#' @param OutFormat as InFormat. Currently raw -> seq, col -> seq and seq ->
#' col are implemented.
#' @param UseFID Use the family ID column in the PLINK file. The
#'   resulting ids (rownames of GenoM) will be in the form FID__IID.
#' @param FIDsep characters inbetween FID and IID in composite-ID. By default a
#'  double underscore is used, to avoid problems when some IIDs contain an
#'   underscore. Only used when UseFID=TRUE.
#' @param quiet suppress messages
#'
#' @return A genotype matrix in the specified output format. If 'OutFile' is
#'   specified, the matrix is written to this file and nothing is returned
#'   inside R.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{LHConvert}, \link{PedStripFID}}
#'
#' @examples
#' \dontrun{
#' # Requires PLINK installed & in system PATH:
#'
#' # tinker with window size, window overlap and VIF to get a set of
#' # 400 - 800 markers (100-200 enough for just parentage):
#' system("cmd", input = "plink --file mydata --indep 50 5 2")
#'
#' system("cmd", input = "plink --file mydata --extract plink.prune.in
#'   --recodeA --out PlinkOUT")
#'
#' GenoM <- GenoConvert(InFile = "PlinkOUT.raw")
#' }
#'
#' @export

GenoConvert <- function(InFile = NULL,
                        InFormat = "raw",
                        OutFile = NA,
                        OutFormat = "seq",
                        UseFID = FALSE,
                        FIDsep = "__",
                        quiet = FALSE) {
  if (OutFormat == "seq" & is.null(OutFile)) {
    OutFile <- "GenoForSequoia.txt"
  } else if (is.null(OutFile)) {
    stop("please provide 'OutFile'")
  }
  if (is.null(InFile)) stop("please provide 'InFile'")
  if (!file.exists(InFile)) stop("cannot find 'InFile'")
  if (UseFID & FIDsep %in% c("", " ", "\t", "\n")) stop("sep can not be whitespace")

  if (interactive() & !quiet & !is.na(OutFile)) {
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

  if (OutFormat == "seq") {
    if (InFormat == "raw") {
      TmpL    <- strsplit(GenoTmp[-1], split = " ")  # skip row w marker names
      GenoOUT <- plyr::ldply(TmpL, function(x) x[-c(1, 3:6)])
      GenoTmp2 <- apply(GenoOUT[, -1], 2,
                        function(x) as.numeric(gsub("NA", "-9", x)))
      if (UseFID) {
        FID <- sapply(TmpL, function(x) x[1])
        rownames(GenoTmp2) <- paste(FID, GenoOUT[, 1], sep=FIDsep)
      } else {
        rownames(GenoTmp2) <- GenoOUT[, 1]
      }
      if (!is.na(OutFile)) {
        utils::write.table(GenoTmp2, file = OutFile,
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      } else {
       return(GenoTmp2)
      }

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




#######################################################################
#######################################################################

#' @title Extract sex and birthyear from PLINK file
#'
#' @description Convert the first six columns of a PLINK .fam, .ped or
#'  .raw file into a three-column lifehistory file for sequoia. Optionally
#'   FID and IID are combined.
#'
#' @details The first 6 columns of PLINK .fam, .ped and .raw files are by
#' default FID - IID - father ID (ignored) - mother ID (ignored) - sex -
#' phenotype.
#'
#' Use with caution, as not extensively tested yet.
#'
#' @param InFile character string with name of genotype file to be converted
#' @param UseFID Use the family ID column. The resulting ids (rownames of
#'  GenoM) will be in the form FID__IID
#' @param SwapSex change the coding from PLINK default (1=male, 2=female) to
#'  sequoia default (1=female, 2=male); any other numbers are set to NA
#' @param FIDsep characters inbetween FID and IID in composite-ID. By default a
#'  double underscore is used, to avoid problems when some IIDs contain an
#'   underscore. Only used when UseFID=TRUE.
#' @param LHIN  dataframe with additional sex and birth year info. In case of
#'   conflicts, LHIN takes priority, with a warning. If UseFID=TRUE, ids are
#'   assumed to be as FID__IID.
#'
#' @return a dataframe with id, sex and birth year, which can be used as input
#'  for \code{\link{sequoia}}
#'
#' @seealso \code{\link{GenoConvert}}, \code{\link{PedStripFID}} to reverse
#'  \code{UseFID}
#'
#' @export

LHConvert <- function(InFile = NULL, UseFID = FALSE,
                      SwapSex = TRUE, FIDsep="__", LHIN=NULL)
{
  if (is.null(InFile)) stop("please provide 'InFile'")
  if (!file.exists(InFile)) stop("cannot find 'InFile'")
  if (UseFID & FIDsep %in% c("", " ", "\t", "\n")) stop("sep can not be whitespace")

  ncol <- length(scan(InFile, nlines=1, what="real", quiet=TRUE))
  TMP <- scan(InFile, skip=1, what=as.list(c(rep("character", 2), rep("numeric", 4),
                             rep("NULL", ncol-6))), quiet=TRUE)

  LH <- data.frame(id = TMP[[2]],
                    Sex = TMP[[5]],
                    BY = TMP[[6]],
                    stringsAsFactors=FALSE)
  if (SwapSex) {
    LH$Sex <- ifelse(LH$Sex==1, 2,
                   ifelse(LH$Sex==2, 1,
                        NA))
  }

  if (UseFID) {
    IDX <- data.frame(id.old = TMP[[2]],
                        id.new = paste(TMP[[1]], TMP[[2]], sep=FIDsep),
                        stringsAsFactors=FALSE)
    LH <- merge(LH, IDX, by.x="id", by.y="id.old", all.x=TRUE)
    LH$id <- ifelse(!is.na(LH$id.new), LH$id.new, LH$id)
    LH <- LH[, c("id", "Sex", "BY")]
  }

  if (!is.null(LHIN)) {
    names(LHIN) <- c("id", "Sex", "BY")
    LH$Sex[!LH$Sex %in% c(1,2)] <- NA
    LHIN$Sex[!LHIN$Sex %in% c(1,2)] <- NA
    LH$BY[LH$BY < 0] <- NA
    LHIN$BY[LHIN$BY < 0] <- NA

    chk <- merge(LH, LHIN, by="id")
    if (any(!is.na(chk$Sex.x) & !is.na(chk$Sex.y) & chk$Sex.x != chk$Sex.y)) {
      warning(paste("There are", sum(chk$Sex.x != chk$Sex.y, na.rm=T), "sex mismatches"))
    } else if (any(!is.na(chk$BY.x) & !is.na(chk$BY.y) & chk$BY.x != chk$BY.y)) {
      warning(paste("There are", sum(chk$BY.x != chk$BY.y, na.rm=T),
                    "birth year mismatches"))
    }
    LH <- Merge(LH, LHIN, by="id", overwrite=TRUE, all=TRUE)
  }

  LH
}

