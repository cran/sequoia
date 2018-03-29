#' @title write sequoia output to excel or text files
#'
#' @description The various list elements returned by \code{sequoia} are each
#' written to
#'  text files in the specified folder, or to separate sheets in a single
#'  excel file (requires library \pkg{xlsx}).
#'
#' @details The text files can be used as input for the stand-alone Fortran
#' version of #'   sequoia, e.g. when the genotype data is too large for R. See
#'   \code{vignette('sequoia')} for further details.
#'
#' @param SeqList the list returned by \code{\link{sequoia}}, to be written out.
#' @param GenoM  the matrix with genetic data (optional). Ignored if
#'  OutFormat='xls', as the resulting file could become too large for excel.
#' @param PedComp a list with results from \code{\link{PedCompare}} (optional).
#'  \code{SeqList$DummyIDs} is combined with \code{PedComp$DummyMatch} if both
#'  are provided.
#' @param OutFormat 'xls' or 'txt'.
#' @param folder the directory where the text files will be written; will be
#'   created if it does not already exists. Relative to the current working
#'   directory, or NULL for current working directory. Ignored if
#'   OutFormat='xls'.
#' @param file the name of the excel file to write to, ignored if
#'  OutFormat='txt'.
#' @param quiet suppress messages.
#'
#'
#' @examples
#' \dontrun{
#' writeSeq(SeqList, OutFormat="xls", file="MyFile.xlsx")
#'
#' # add additional sheets to the excel file:
#' library(xlsx)
#' write.xlsx(MyData, file = "MyFile.xlsx", sheetName="ExtraData",
#'       col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
#' }
#'
#' @export

writeSeq <- function(SeqList, GenoM=NULL, PedComp=NULL,
                      OutFormat="txt",
                      folder="Sequoia-OUT",
                      file="Sequoia-OUT.xlsx", quiet=FALSE) {
  if (OutFormat == "xlsx") OutFormat <- "xls"
  if (OutFormat == "xls") {
    if (!requireNamespace("xlsx", quietly = TRUE)) {
      if (interactive() & !quiet) {
        ANS <- readline(prompt = paste("library 'xlsx' not found. Install Y/N? "))
        if (!substr(ANS, 1, 1) %in% c("Y", "y", "J", "j", "")) stop()
      }
      utils::install.packages("xlsx")
   }
   write.seq.xls(SeqList, file=file, PedComp=PedComp, quiet=quiet)

  } else if (OutFormat == "txt") {

  curdir <- getwd()
  if (is.null(folder)) folder = curdir
  dir.create(folder, showWarnings = FALSE)
  setwd(folder)
  if (any(file.exists("Geno.txt","Specs.txt", "AgePriors.txt", "PedigreePar.txt")) &
      interactive() & !quiet) {
    ANS <- readline(prompt = paste("Writing data to '", folder,
                                   "' will overwrite existing files. Continue Y/N? "))
    if (!substr(ANS, 1, 1) %in% c("Y", "y", "J", "j", "")) stop()
  }
  write(paste("These files were created by R package sequoia on ", date()),
        file="README.txt")

  if (!is.null(GenoM)) {
    if(nrow(GenoM)!= SeqList$Specs$NumberIndivGenotyped) {
      ANS <- readline(prompt = paste("Number of individuals according to Specs differs from number of rows in GenoM (", SeqList$Specs$NumberIndivGenotyped, "/", nrow(GenoM), ").\n Press Y to continue and fix manually in `SequoiaSpecs.txt' "))
      if (!substr(ANS, 1, 1) %in% c("Y", "y")) stop()
    }
    if(ncol(GenoM)!= SeqList$Specs$NumberSnps) {
      stop(paste("Number of SNPs according to Specs differs from number of rows in GenoM (", SeqList$Specs$NNumberSnps, "vs", ncol(GenoM), ")"))
    }
  }

  SpecsOUT <- cbind("Genofile" = "Geno.txt",
                    "LHfile" = "LifeHist.txt",
                    "nIndLH" = nrow(SeqList$LifeHist),
                    SeqList$Specs)
  SpecsOUT$Complexity <- c("full"=2, "simp"=1, "mono"=0)[SpecsOUT$Complexity]
  SpecsOUT$UseAge <- c("extra"=2, "yes"=1, "no"=0)[SpecsOUT$UseAge]
  SpecsOUT$FindMaybeRel <- as.numeric(SpecsOUT$FindMaybeRel)
  SpecsOUT$CalcLLR <- as.numeric(SpecsOUT$CalcLLR)

  OPT <- options()
  options(scipen = 10)
  if (!is.null(GenoM)) utils::write.table(GenoM, file="Geno.txt",
              quote=FALSE, row.names=TRUE, col.names=FALSE)
  utils::write.table(as.data.frame(t(SpecsOUT)), file="SequoiaSpecs.txt",
              sep = "\t,\t", quote = FALSE, col.names = FALSE)
  writeColumns(SeqList$AgePrior, "AgePriors.txt", row.names=FALSE)
  writeColumns(SeqList$LifeHist, "LifeHist.txt", row.names=FALSE)
  if ("PedigreePar" %in% names(SeqList)) {
    Par <- SeqList$PedigreePar
    for (x in c("dam", "sire")) Par[is.na(Par[, x]), x] <- "NA"
    for (x in c("id", "dam", "sire")) {
      Par[, paste0("row", x)] <- as.numeric(factor(Par[,x], levels=rownames(GenoM)))
      Par[is.na(Par[, paste0("row", x)]), paste0("row", x)] <- 0
    }
    for (x in c("LLRdam", "LLRsire", "LLRpair", "OHdam", "OHsire")) {
      Par[is.na(Par[, x]), x] <- -9
    }
    if (any(is.na(Par$id))) stop("Some id's in PedigreePar do not occur in GenoM!")
    writeColumns(Par, "Parents.txt", row.names=FALSE)
  }

  if ("Pedigree" %in% names(SeqList)) {
    writeColumns(SeqList$Pedigree, "Pedigree.txt", row.names=FALSE)
  }
  if ("DummyIDs" %in% names(SeqList)) {
    if ("DummyMatch" %in% names(PedComp)) {
      Dummies <- merge(PedComp$DummyMatch, SeqList$DummyIDs)
    } else {
      Dummies <- SeqList$DummyIDs
    }
    writeColumns(Dummies, "DummyIDs.txt", row.names=FALSE)
  }

  if (!is.null(PedComp)) {
    Counts.out <- cbind(Parent = rep(c("dam", "sire"), each=7),
                        Cat = rep(dimnames(PedComp$Counts)[[1]], 2))
    Counts.out <- cbind(Counts.out, rbind(PedComp$Counts[,,"dam"],
                                      PedComp$Counts[,,"sire"]))
    writeColumns(Counts.out, "PedCompare-Counts.txt", row.names=FALSE)

    writeColumns(PedComp$MergedPed, "PedCompare-MergedPed.txt",
                     row.names=FALSE)
  }

  options(OPT)
  setwd(curdir)
#  if(!quiet) message(paste("Output written to", folder))
  } else {
    stop("OutFormat not supported.")
  }
}




#==========================================================================


write.seq.xls <- function(SeqList, file, PedComp=NULL, quiet) {

  if (!requireNamespace("xlsx", quietly = TRUE)) {
      stop("library 'xlsx' not found")
   }

  # check if file exists
  if (file.exists(file) & interactive() & !quiet) {
    ANS <- readline(prompt = paste("Writing data to '", file,
                                   "' will overwrite existing file. Continue Y/N? "))
    if (!substr(ANS, 1, 1) %in% c("Y", "y", "J", "j", "")) stop()
  }

  if ("Pedigree" %in% names(SeqList)) {
    xlsx::write.xlsx(SeqList$Pedigree, file = file, sheetName="Pedigree",
             col.names=TRUE, row.names=FALSE, append=FALSE, showNA=FALSE)
    if ("PedigreePar" %in% names(SeqList)) {
      xlsx::write.xlsx(SeqList$PedigreePar, file = file, sheetName="Parents",
             col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
    }
  } else if ("PedigreePar" %in% names(SeqList)) {
    xlsx::write.xlsx(SeqList$PedigreePar, file = file, sheetName="Parents",
             col.names=TRUE, row.names=FALSE, append=FALSE, showNA=FALSE)
  } else {
    xlsx::write.xlsx("None", file = file, sheetName="Parents",
             col.names=TRUE, row.names=FALSE, append=FALSE, showNA=FALSE)
  }
  xlsx::write.xlsx(t(SeqList$Specs), file = file, sheetName="Run parameters",
             col.names=FALSE, row.names=TRUE, append=TRUE, showNA=FALSE)

  Names <- c("AgePriors", "LifeHist", "DupGenotype", "DupLifeHistID",
             "DummyIDs", "MaybeParent", "MaybeRel", "TotLikParents", "TotLikSib")
  for (x in Names) {
    if (!x %in% names(SeqList))  next
    if (is.null(dim(SeqList[[x]])))  next
    if (nrow(SeqList[[x]])==0)  next
    if (x == "DummyIDs") {
      if ("DummyMatch" %in% names(PedComp)) {
        SeqList$DummyIDs <- merge(PedComp$DummyMatch, SeqList$DummyIDs)
      }
    }
    xlsx::write.xlsx(SeqList[[x]], file = file, sheetName = x,
           col.names=TRUE, row.names=ifelse(x=="AgePriors", TRUE, FALSE),
           append=TRUE, showNA=FALSE)
  }

  if (!is.null(PedComp)) {  # pedcompare output
    Counts.out <- cbind(Parent = rep(c("dam", "sire"), each=7),
                        Cat = rep(dimnames(PedComp$Counts)[[1]], 2))
    Counts.out <- cbind(Counts.out, rbind(PedComp$Counts[,,"dam"],
                                      PedComp$Counts[,,"sire"]))
    rownames(Counts.out) <- 1:nrow(Counts.out)  # not used, but otherwise warning

    xlsx::write.xlsx(Counts.out, file = file, sheetName = "PedCompare-Counts",
                     col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
    xlsx::write.xlsx(PedComp$MergedPed, file=file, sheetName="PedCompare-MergedPed",
                     col.names=FALSE, row.names=FALSE, append=TRUE, showNA=FALSE)
  }
#  if(!quiet) message(paste("Output written to", file))
}




