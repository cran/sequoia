#' @title Convert genotype data
#'
#' @description Convert genotype data in various formats to sequoia's
#'   1-column-per-marker format or Colony's 2-column-per-marker format.
#'
#' @section Input formats:
#' The following formats can be specified by \code{InFormat}:
#' \describe{
#'   \item{single}{1 column per marker, otherwise unspecified}
#'   \item{double}{2 columns per marker, otherwise unspecified}
#'   \item{col}{(Colony) genotypes are coded as numeric values, missing as 0, in
#'   2 columns per marker. Column 1 contains IDs.}
#'   \item{ped}{(PLINK) genotypes are coded as A, C, T, G, missing as 0, in 2
#'   columns per marker. The first 6 columns are descriptive (1:FID, 2:IID, 3 to
#'   6 ignored). }
#'   \item{raw}{(PLINK) genotypes are coded as 0, 1, 2, missing as NA, in 1
#'   column per marker. The first 6 columns are descriptive (1:FID, 2:IID, 3 to
#'   6 ignored), and there is a header row.}
#'   \item{seq}{(sequoia) genotypes are coded as 0, 1, 2, missing as \eqn{-9},
#'   in 1 column per marker. Column 1 contains IDs, there is no header row.}
#'   }
#'
#'  For each InFormat, its default values for \code{Missing, header, IDcol,
#'  FIDcol}, and \code{dropcol} can be overruled by specifying the corresponding
#'  input parameters.
#'
#' @param InFile character string with name of genotype file to be converted
#' @param InData dataframe or matrix with genotypes to be converted
#' @param InFormat One of 'single', 'double', 'col', 'ped', 'raw', or 'seq', see
#'   Details.
#' @param OutFile character string with name of converted file. If NA, return
#'   matrix with genotypes in console (default); if NULL, write to
#'   'GenoForSequoia.txt' in current working directory.
#' @param OutFormat as InFormat, currently only 'seq' and 'col' are implemented.
#' @param Missing vector with symbols interpreted as missing data.
#' @param sep vector with field separator strings that will be tried on
#'   \code{InFile}. The OutFile separator uses the write.table default, i.e. one
#'   blank space
#' @param header a logical value indicating whether the file contains the names
#'   of the variables as its first line. If NA (default), set to TRUE for 'raw',
#'   and FALSE otherwise.
#' @param IDcol single number giving the column which contains the individual
#'   IDs; 0 indicates the rownames (for InData only). If NA (default), set to 2
#'   for InFormat 'raw' and 'ped', and otherwise to 1 for InFile and 0
#'   (rownames) for InData, except when InData has a column labeled 'ID'.
#' @param FIDcol  column which contains the individual IDs, if any are wished to
#'   be used. This is column 1 for InFormat 'raw' and 'seq', but those are by
#'   default not used.
#' @param FIDsep string used to paste FID and IID together into a composite-ID
#'   (value passed to \code{paste}'s \code{collapse}). This joining can be
#'    reversed using \code{\link{PedStripFID}}.
#' @param dropcol  columns to exclude from the output data, on top of IDcol and
#'   FIDcol (which become rownames). When NA, defaults to columns 3-6 for
#'   InFormat 'raw' and 'seq'. Can also be used to drop some SNPs, see example
#'   below on how to do this for the 2-columns-per-SNP input formats.
#' @param quiet suppress messages and warnings
#'
#' @return A genotype matrix in the specified output format. If 'OutFile' is
#'   specified, the matrix is written to this file and nothing is returned
#'   inside R. When converting to 0/1/2 format, 2 is the homozygote for the
#'   minor allele, and 0 the homozygote for the major allele.
#'
#' @section Error messages: An occasional error when reading in a file with
#'   GenoConvert is that 'rows have unequal length'. GenoConvert makes use of
#'   \code{\link{readLines}} and \code{\link{strsplit}}, which is much faster
#'   than \code{\link{read.table}} for large datafiles, but also more sensitive
#'   to unusual line endings, unusual end-of-file characters, or invisible
#'   characters (spaces or tabs) after the end of some lines. In these cases,
#'   try to read the data from file using read.table or read.csv, and then use
#'   GenoConvert on the matrix, see example.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{CheckGeno}, \link{SnpStats}, \link{LHConvert}}
#'
#' @examples
#' \dontrun{
#' # Requires PLINK installed & in system PATH:
#'
#' # tinker with window size, window overlap and VIF to get a set of
#' # 400 - 800 markers (100-200 enough for just parentage):
#' system("cmd", input = "plink --file mydata --indep 50 5 2")
#' system("cmd", input = "plink --file mydata --extract plink.prune.in
#'   --recodeA --out PlinkOUT")
#'
#' GenoM <- GenoConvert(InFile = "PlinkOUT.raw")
#'
#' # save time on file conversion next time:
#' write.table(GenoM, file="Geno_for_sequoia.txt", quote=FALSE,
#'   col.names=FALSE)
#' GenoM <- read.table("Geno_for_sequoia.txt", row.names=1, header=FALSE)
#'
#' # drop some SNPs, e.g. after a warning of >2 alleles:
#' dropSNP <- c(5,68,101,128)
#' GenoM <- GenoConvert(ColonyFile, InFormat = "col",
#'                      dropcol = 1 + c(2*dropSNP-1, 2*dropSNP) )
#'
#' # circumvent a 'rows have unequal length' error:
#' GenoTmp <- as.matrix(read.table("mydata.txt", header=TRUE, row.names=1))
#' GenoM <- GenoConvert(InData=GenoTmp, InFormat="single", IDcol=0)
#' }
#'
#' @importFrom stats na.exclude
#'
#' @export

GenoConvert <- function(InFile = NULL,
                        InFormat = "raw",
                        OutFile = NA,
                        OutFormat = "seq",
                        InData = NULL,
                        Missing = c("-9", "??", "?", "NA", "NULL",
                                    c("0")[InFormat %in% c("col","ped")]),
                        sep = c(" ", "\t", ",", ";"),
                        header = NA,
												IDcol = NA, #
                        FIDcol = NA,
                        FIDsep = "__",
                        dropcol = NA,
                        quiet = FALSE) {

  if (is.null(InFile) & is.null(InData)) {
    stop("please provide 'InFile' or 'InData'")

  } else if (!is.null(InFile) & !is.null(InData)) {
    stop("please provide either 'InFile' or 'InData', not both")

  } else if (!is.null(InFile)) {
    if (is.character(InFile)) {
      if (!file.exists(InFile)) stop("cannot find 'InFile'")
    } else if (is.matrix(InFile) | is.data.frame(InFile)) {
      InData <- InFile
      InFile <- NULL
    } else {
      stop("'InFile' in unknown format, should be character string.")
    }
    if (is.na(header)) {
      header <- ifelse(InFormat == "raw", TRUE, FALSE)
    }
  }

  if (!InFormat %in% c("raw", "seq", "col", "ped", "single", "double")) {
    stop("invalid InFormat")
  }
  if (OutFormat %in% c("raw", "ped", "single")) {
    stop("OutFormat not (yet) implemented")
  } else if (!OutFormat %in% c("seq", "col")) {
    stop("invalid OutFormat")
  }

  if (!is.na(FIDcol) & FIDsep %in% c("", " ", "\t", "\n")) stop("sep can not be whitespace")
  UseFID <- ifelse(!is.na(FIDcol), TRUE, FALSE)

  if (is.na(IDcol)) {
    IDcol <- ifelse(InFormat %in% c("raw", "ped"),
                2,
                ifelse(!is.null(InFile),
                       1,
                       ifelse("ID" %in% colnames(InData),
                          which(colnames(InData) == "ID"),
                          0)))    # rownames
  }
  if (InFormat %in% c("raw", "ped")) {
    if(is.na(FIDcol))   FIDcol <- 1
    if(is.na(dropcol))  dropcol <- c(3:6)
  }

  if (OutFormat == "seq" & is.null(OutFile)) {
    OutFile <- "GenoForSequoia.txt"

  } else if (is.null(OutFile)) {
    stop("please provide 'OutFile'")

  }
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

  #~~~~~~~~~~~~~~~~~~~~~~~~

  if (!is.null(InData)) {
    GenoTmp <- as.matrix(InData)
    rm(InData)

  } else if (!is.null(InFile)) {
    GenoTmp <- readLines(InFile, warn=FALSE)
    if (header)  GenoTmp <- GenoTmp[-1]

    TmpL <- strsplit(GenoTmp, split = sep[1])
    if (length(TmpL[[1]])==1) {
      for (s in sep[-1]) {
        TmpL <- strsplit(GenoTmp, split = s)
        if (length(TmpL[[1]]) > 1) break
      }
    }
    if (length(TmpL[[1]])==1) {
      stop("unknown column separator, expecting ' ' (space), '\t' (tab), ',' or ';'")
    }

    if (length(table(sapply(TmpL, length))) > 1) {
      stop("rows have unequal length")
    }
    GenoTmp <- plyr::ldply(TmpL)
    rm(TmpL)
  }
  if (nrow(GenoTmp)<2)  stop("Genotype matrix must have at least 2 individuals")
  if (ncol(GenoTmp)<2)  stop("Genotype matrix must have at least 2 SNPs")

  #~~~~~~~~~

  if (IDcol==0) {
    IDs_geno <- rownames(GenoTmp)
  } else {
    IDs_geno <- GenoTmp[, IDcol]
  }
  if (!is.na(FIDcol))  FID <- GenoTmp[, FIDcol]
  dropcol <- na.exclude(c(FIDcol, IDcol, dropcol))
  if (any(dropcol!=0))   GenoTmp <- GenoTmp[, -dropcol]
  GenoTmp <- as.matrix(GenoTmp)

  if (UseFID) {
    IDs_geno <- paste(FID, IDs_geno, sep=FIDsep)
  }

  for (misX in Missing) {
    if (any(GenoTmp == misX, na.rm=TRUE)) {
      GenoTmp <- array(gsub(misX, NA, GenoTmp, fixed=TRUE),
                        dim=dim(GenoTmp))
    }
  }

  if (InFormat %in% c("raw", "seq")) {
    if (any(!GenoTmp %in% c(0,1,2, NA))) {
      stop(paste("Unexpected value! When InFormat=", InFormat, ", genotypes should be coded as 0/1/2.\n",
           "Choose InFormat='single' and/or different Missing, or check data."))
    }
  } else if (InFormat =="single") {
     if(!all(nchar(GenoTmp)==2, na.rm=T)) {
       stop(paste("Unexpected value! When InFormat='single', genotypes should be coded as 2 digits or characters.\n",
                  "Choose InFormat='col' for 2-columns-per-SNP format, or check data."))
     }
  }

  if (InFormat %in% c("col", "ped", "single", "double")) {  # A/C/T/G -> 0/1/2
    GCA <- array(dim=c(2, nrow(GenoTmp), ncol(GenoTmp)/2))
    if (InFormat %in% c("ped", "col", "double")) {
      GCA[1,,] <- GenoTmp[, seq(1,ncol(GenoTmp)-1,2)]
      GCA[2,,] <- GenoTmp[, seq(2,ncol(GenoTmp),2)]
    } else {
      GCA[1,,] <- substr(GenoTmp,1,1)
      GCA[2,,] <- substr(GenoTmp,2,2)
    }
    Alleles <- apply(GCA, 3, function(x) table(factor(x,
                                          levels=sort(na.exclude(unique(c(GenoTmp)))))))
    if (is.matrix(Alleles)) {
      NumAlleles <- apply(Alleles, 2, function(x) sum(x>0))
      minorAllele <- apply(Alleles, 2, function(x) names(sort(x[x>0]))[1])
    } else {  # list
      NumAlleles <- sapply(Alleles, length)
      minorAllele <- sapply(Alleles, function(x)  names(sort(x))[1])
    }
    if (any(NumAlleles > 2)) {
      n.problem <- sum(NumAlleles > 2)
      stop(paste("There are", n.problem, "SNPs with >2 alleles  ",
                 ifelse(n.problem<=10, paste(which(sapply(Alleles, length) > 2), collapse="-"),"")))
    }
    if (any(NumAlleles ==1) & !quiet) {
      warning(paste("There are", sum((NumAlleles ==1)), "monomorphic SNPs"))
    }

    GenoTmp2 <- sapply(1:dim(GCA)[3], function(i) {
      apply(GCA[,,i], 2, function(x) ifelse(is.na(x[1]) | is.na(x[2]), NA,
                                         ifelse(x[1] != x[2], 1,  # heterozygote
                                            ifelse(x[1] == minorAllele[i], 2, 0))) ) } )

  } else {
    AllHom0 <- apply(GenoTmp, 2, function(x) all(na.exclude(x) == 0))
    AllHom2 <- apply(GenoTmp, 2, function(x) all(na.exclude(x) == 2))
    if ((any(AllHom0) | any(AllHom2)) & !quiet) {
      warning(paste("There are", sum(AllHom0)+sum(AllHom2), "monomorphic SNPs"))
    }

    GenoTmp2 <- matrix(as.numeric(GenoTmp), nrow(GenoTmp))
  }
  rownames(GenoTmp2) <- IDs_geno
  rm(GenoTmp)

  if (OutFormat == "seq") {
    GenoOUT <- GenoTmp2
    GenoOUT[is.na(GenoOUT)] <- -9
    CheckGeno(GenoOUT, quiet=quiet)

  } else if (OutFormat == "col") {
    dc <- list("0" = c(1,1), "1" = c(1,2), "2" = c(2,2), "-9" = c(0,0))
    GenoTmp2[is.na(GenoTmp2)] <- -9
    GenoA <- array(dim=c(nrow(GenoTmp2), 2, ncol(GenoTmp2)))
    for (i in 1:nrow(GenoTmp2)) {
      GenoA[i,,] <- sapply(GenoTmp2[i,], function(z) dc[[z]])
    }
    GenoOUT <- matrix(GenoA, nrow(GenoTmp2))
    row.names(GenoOUT) <- IDs_geno

#  } else {
#    stop("OutFormat not implemented")  # caught above
  }

  if (!is.na(OutFile)) {
    utils::write.table(GenoOUT, file = OutFile,
              row.names = TRUE, col.names = FALSE, quote = FALSE)
  } else {
   return(GenoOUT)
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
#' When additionally a
#'
#' @param PlinkFile character string with name of genotype file to be converted
#' @param UseFID Use the family ID column. The resulting ids (rownames of GenoM)
#'   will be in the form FID__IID
#' @param SwapSex change the coding from PLINK default (1=male, 2=female) to
#'   sequoia default (1=female, 2=male); any other numbers are set to NA
#' @param FIDsep characters inbetween FID and IID in composite-ID. By default a
#'   double underscore is used, to avoid problems when some IIDs contain an
#'   underscore. Only used when UseFID=TRUE.
#' @param LifeHistData  dataframe with additional sex and birth year info. In
#'   case of conflicts, LifeHistData takes priority, with a warning. If
#'   UseFID=TRUE, IDs in LifeHistData are assumed to be already as FID__IID.
#'
#' @return a dataframe with id, sex and birth year, which can be used as input
#'  for \code{\link{sequoia}}
#'
#' @seealso \code{\link{GenoConvert}}, \code{\link{PedStripFID}} to reverse
#'  \code{UseFID}
#'
#' @examples
#' \dontrun{
#' # combine FID and IID in dataframe with additional sex & birth years
#' ExtraLH$FID_IID <- paste(ExtraLH$FID, ExtraLH$IID, sep = "__")
#' LH.new <- LHConvert(PlinkFile, UseFID = TRUE, FIDsep = "__",
#'                     LifeHistData = ExtraLH)
#' }
#'
#' @export

LHConvert <- function(PlinkFile = NULL, UseFID = FALSE,
                      SwapSex = TRUE, FIDsep="__", LifeHistData=NULL)
{
  if (is.null(PlinkFile)) stop("please provide 'InFile'")
  if (!file.exists(PlinkFile)) stop("cannot find 'PlinkFile'")
  if (UseFID & FIDsep %in% c("", " ", "\t", "\n")) stop("sep can not be whitespace")
  if (!is.null(LifeHistData)) {
    LHIN <- CheckLH(LifeHistData)
  }

  ncol <- length(scan(PlinkFile, nlines=1, what="real", quiet=TRUE))
  TMP <- scan(PlinkFile, skip=1, what=as.list(c(rep("character", 2), rep("numeric", 4),
                             rep("NULL", ncol-6))), quiet=TRUE)

  LH <- data.frame(id = TMP[[2]],
                   Sex = TMP[[5]],
                   BirthYear = TMP[[6]],
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
    LH <- LH[, c("id", "Sex", "BirthYear")]
  }

  if (!is.null(LHIN)) {
    names(LHIN) <- c("id", "Sex", "BirthYear")
    LH$Sex[!LH$Sex %in% c(1,2)] <- NA
    LHIN$Sex[!LHIN$Sex %in% c(1,2)] <- NA
    LH$BirthYear[LH$BirthYear < 0] <- NA
    LHIN$BirthYear[LHIN$BirthYear < 0] <- NA

    chk <- merge(LH, LHIN, by="id")
    n.sexmismatch <- sum(chk$Sex.x != chk$Sex.y, na.rm=T)
    n.BYmismatch <- sum(chk$BirthYear.x != chk$BirthYear.y, na.rm=T)
    if (n.sexmismatch > 0 & n.sexmismatch <= 10) {
      these <- with(chk, id[which(!is.na(Sex.x) & !is.na(Sex.y) & Sex.x!=Sex.y)])
      warning(paste("There are", n.sexmismatch, "sex mismatches: ",
                    paste(these, collapse=", ")))
    } else if (n.sexmismatch>10) {
      warning(paste("There are", n.sexmismatch, "sex mismatches"))
    }
    if (n.BYmismatch > 0 & n.BYmismatch <= 10) {
      these <- with(chk, id[which(!is.na(BirthYear.x) & !is.na(BirthYear.y) & BirthYear.x!=BirthYear.y)])
      warning(paste("There are", n.BYmismatch, "birth year mismatches: ",
                    paste(these, collapse=", ")))
    } else if (n.BYmismatch>10) {
      warning(paste("There are", n.BYmismatch, "BY mismatches"))
    }

    LH <- MergeFill(LH, LHIN, by="id", overwrite=TRUE, all=TRUE)
  }

  LH
}

