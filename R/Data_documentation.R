#' @title Example pedigree
#'
#' @description This is \strong{Pedigree II} in the paper, with discrete
#'   generations and considerable inbreeding
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. (2017) Pedigree reconstruction from SNP data:
#'   Parentage assignment, sibship clustering, and beyond. Molecular Ecology
#'   Resources 17:1009--1024.
#'
#' @seealso \code{\link{LH_HSg5} \link{SimGeno_example} \link{sequoia}}
#'
#' @docType data
#' @keywords datasets sequoia
#' @name Ped_HSg5
#' @usage data(Ped_HSg5)
#' @format A data frame with 1000 rows and 3 variables (id, dam, sire)
NULL


#===============================================================================
#' @title Example life history file
#'
#' @description This is the lifehistory file associated with
#'   \code{\link{Ped_HSg5}}, which is \strong{Pedigree II} in the paper.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. (2017) Pedigree reconstruction from SNP data:
#'   Parentage assignment, sibship clustering, and beyond. Molecular Ecology
#'   Resources 17:1009--1024.
#'
#' @seealso \code{\link{Ped_HSg5} \link{sequoia}}
#'
#' @docType data
#' @keywords datasets sequoia
#' @name LH_HSg5
#' @usage data(LH_HSg5)
#' @format A data frame with 1000 rows and 3 variables:
#' \describe{
#' \item{ID}{Female IDs start with 'a', males with 'b'; the next 2 numbers give
#' the generation number (00 -- 05), the last 3 numbers the individual ID number
#' (runs continuously across all generations)}
#' \item{Sex}{1 = female, 2 = male}
#' \item{BirthYear}{from 2000 (generation 0, founders) to 2005} }
#'
NULL


#===============================================================================
#' @title Example genotype file
#'
#' @description Simulated genotype data for cohorts 1+2 in Pedigree Ped_HSg5
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{Ped_HSg5}, \link{SimGeno}}
#'
#' @docType data
#' @keywords datasets sequoia
#' @name SimGeno_example
#' @usage data(SimGeno_example)
#' @format A genotype matrix with 214 rows (ids) and 200 columns (SNPs). Each
#'   SNP is coded as 0/1/2 copies of the reference allele, with -9 for missing
#'   values. Ids are stored as rownames.
NULL


#===============================================================================
#===============================================================================
#' @title Inheritance patterns
#'
#' @description Inheritance patterns used by SimGeno for non-autosomal SNPs,
#'   identical to those in Inherit.xlsx
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{SimGeno}}
#'
#' @docType data
#' @keywords datasets sequoia inherit
#' @name Inherit
#' @usage data(Inherit)
#' @format An array with the following dimensions:
#' \describe{
#'   \item{d1}{type: autosomal, x-chromosome, y-chromosome, or mtDNA}
#'   \item{d2}{offspring sex: female, male, or unknown}
#'   \item{d3}{offspring genotype: aa (0), aA (1), Aa (1), or AA (2)}
#'   \item{d4}{mother genotype}
#'   \item{d5}{father genotype}
#' }
NULL


#===============================================================================
#===============================================================================
#' @title Example pedigree: griffins
#'
#' @description Example pedigree used in the ageprior vignette, with overlapping
#'   generations.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso  \code{\link{LH_griffin}}; \code{\link{SeqOUT_griffin}} for a sequoia
#'   run on simulated genotype data based on this pedigree;
#'   \code{\link{Ped_HSg5}} for another pedigree, \code{\link{sequoia}}
#'
#' @docType data
#' @keywords datasets sequoia
#' @name Ped_griffin
#' @usage data(Ped_griffin)
#' @format A data frame with 200 rows and 4 variables (id, dam, sire, birthyear)
NULL


#===============================================================================
#' @title Example lifehistory data: griffins
#'
#' @description Example lifehistory data for griffin pedigree
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{Ped_griffin}}, \code{\link{SeqOUT_griffin}}
#'
#' @docType data
#' @keywords datasets sequoia
#' @name LH_griffin
#' @usage data(LH_griffin)
#' @format A data frame with 200 rows and 3 variables (ID, Sex, BirthYear)
NULL


#===============================================================================
#' @title Example sequoia output (griffins)
#'
#' @description Example output of a sequoia run including sibship clustering,
#' based on the griffin pedigree.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{Ped_griffin}, \link{sequoia}}
#'
#' @docType data
#' @keywords datasets sequoia
#' @name SeqOUT_griffin
#' @usage data(SeqOUT_griffin)
#' @format a list, see \code{\link{sequoia}}
#'
#' @examples
#' \dontrun{
#' GenoS <- SimGeno(Ped.griffin, nSnp=400, ParMis=0.4)
#' griffin.sex <- sapply(as.character(Ped_griffin$id), function(x)
#'                       substr(x, start=nchar(x), stop=nchar(x)))
#' LH.griffin <- data.frame(ID = Ped_griffin$id,
#'                          Sex = ifelse(griffin.sex=="F", 1, 2),
#'                          BirthYear = Ped_griffin$birthyear)
#' SeqOUT.GX <- sequoia(GenoS, LH.griffin,
#'                      Module = "ped",
#'                      args.AP = list(Smooth = FALSE))
#' }
NULL


#===============================================================================
#' @title Example field-observed mothers: griffins
#'
#' @description Example field pedigree used in vignette for
#'   \code{\link{PedCompare}} example. Non-genotyped females have IDs 'BlueRed',
#'   'YellowPink', etc.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{SeqOUT_griffin}} for a sequoia run on simulated genotype
#'   data, \code{\link{Ped_griffin}} for the 'true' pedigree.
#'
#' @docType data
#' @keywords datasets sequoia
#' @name FieldMums_griffin
#' @usage data(FieldMums_griffin)
#' @format A data frame with 144 rows and 2 variables (id, mum)
NULL


#===============================================================================
