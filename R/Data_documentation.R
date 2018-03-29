#' @title Example pedigree
#'
#' @description This is Pedigree II in the paper.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. (2017) Pedigree reconstruction from SNP data: Parentage
#'   assignment, sibship clustering, and beyond. Molecular Ecology Resources
#'   17:1009--1024.
#'
#' @seealso \code{\link{LH_HSg5} \link{SimGeno_example} \link{sequoia}}
#'
#' @docType data
#' @keywords datasets, sequoia
#' @name Ped_HSg5
#' @usage data(Ped_HSg5)
#' @format A data frame with 1000 rows and 3 variables (id, dam, sire)
NULL



#' @title Example life history file
#'
#' @description This is the lifehistory file associated with Ped_HSg5, which
#'  is Pedigree II in the paper.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. (2017) Pedigree reconstruction from SNP data: Parentage
#'   assignment, sibship clustering, and beyond. Molecular Ecology Resources
#'   17:1009--1024.
#'
#' @seealso \code{\link{Ped_HSg5} \link{sequoia}}
#'
#' @docType data
#' @keywords datasets, sequoia
#' @name LH_HSg5
#' @usage data(LH_HSg5)
#' @format A data frame with 1000 rows and 3 variables: ID, Sex (1=female,
#' 2=male), and BY (birth year, here cohort)
NULL



#' @title Example genotype file
#'
#' @description Simulated genotype data for cohorts 1+2 in Pedigree Ped_HSg5
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{Ped_HSg5}, \link{SimGeno}}
#'
#' @docType data
#' @keywords datasets, sequoia
#' @name SimGeno_example
#' @usage data(SimGeno_example)
#' @format A data frame with 214 rows and 201 columns: id, followed by 1 column
#' per SNP coded as 0/1/2 or -9 for missing values.
NULL
