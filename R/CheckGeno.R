#=======================================================================
#' @title check GenoM
#'
#' @description Check that the provided genotype matrix is in the correct format
#'
#' @param GenoM the genotype matrix
#' @param quiet suppress messages
#'
#' @return a list with low call rate individuals (<0.5, 'ExcludedInd') and low
#'   call rate SNPs (<0.1, 'ExcludedSnps'), if any, which will be excluded from
#'   pedigree reconstruction.
#'
#' @export

CheckGeno <- function(GenoM, quiet=FALSE) {
	if (is.null(GenoM)) stop("please provide 'GenoM'")
  if (!is.matrix(GenoM)) stop("'GenoM' should be a numeric matrix")
	if (!all(GenoM %in% c(0,1,2,-9))) {
	  UniqueValues <- unique(c(GenoM))
	  InvalidValues <- UniqueValues[!UniqueValues %in% c(0,1,2,-9)]
	  stop(paste0("'GenoM' includes invalid values: ", "'",
	             paste(InvalidValues, collapse="', '"), "'"))
	}

  if (is.null(rownames(GenoM))) stop("'GenoM' has no rownames, these should be the individual IDs")
	if (any(duplicated(rownames(GenoM))))  stop("'GenoM' has duplicate IDs. Please exclude or rename these samples, or run GenoConvert with UseFID=TRUE.")

  Excl <- list()
	Lscored <- apply(GenoM, 1, function(x) sum(x!=-9))
  if (any(Lscored < ncol(GenoM)/2)) {
    warning(paste("There are ", sum(Lscored < ncol(GenoM)/2)," individuals scored for <50% of SNPs, these will be excluded"),
            immediate.=TRUE)
		GenoM <<- GenoM[Lscored >= ncol(GenoM)/2, ]
		Excl[["ExcludedInd"]] <- rownames(GenoM)[which(Lscored < ncol(GenoM)/2)]
  }

  Nscored <- apply(GenoM, 2, function(x) sum(x!=-9))
  if (any(Nscored < nrow(GenoM)/10)) {
    warning(paste("There are ", sum(Nscored < nrow(GenoM)/10)," SNPs scored for <10% of individuals, these will be excluded"),
            immediate.=TRUE)
		GenoM <<- GenoM[, Nscored >= ncol(GenoM)/10]
		Excl[["ExcludedSnps"]] <- which(Nscored < nrow(GenoM)/10)
  }
  if (all(Lscored > ncol(GenoM)/2) & all(Nscored > nrow(GenoM)/10)) {
    if (!quiet)  message("'GenoM' looks OK!")
  }

  invisible( Excl )
}
