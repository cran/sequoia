#' @title Example pedigree
#'
#' @description This is Pedigree II in the paper, with discrete generations and
#'   considerable inbreeding
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



#' @title Example life history file
#'
#' @description This is the lifehistory file associated with Ped_HSg5, which is
#'   Pedigree II in the paper.
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
#' @format A data frame with 1000 rows and 3 variables: ID, Sex (1=female,
#' 2=male), and BirthYear
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
#' @keywords datasets sequoia
#' @name SimGeno_example
#' @usage data(SimGeno_example)
#' @format A genotype matrix with 214 rows (ids) and 200 columns (SNPs). Each
#'   SNP is coded as 0/1/2 copies of the reference allele, with -9 for missing
#'   values. Ids are stored as rownames.
NULL


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


#' @title Example pedigree: griffins
#'
#' @description Example Pedigree used in the ageprior vignette, with overlapping
#'   generations.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{SeqOUT_griffin}} for a sequoia run on simulated genotype
#'   data based on this pedigree; \code{\link{Ped_HSg5}} for another pedigree,
#'   \code{\link{sequoia}}
#'
#' @docType data
#' @keywords datasets sequoia
#' @name Ped_griffin
#' @usage data(Ped_griffin)
#' @format A data frame with 200 rows and 4 variables (id, dam, sire, birthyear)
NULL


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
#' griffin.sex <- sapply(Ped.griffin$ID,
#'                     function(x) substr(x, start=nchar(x), stop=nchar(x)))
#' LH_griffin <- data.frame(ID = Ped_griffin$ID,
#'                          Sex = ifelse(griffin.sex=="F", 1, 2),
#'                          BirthYear = Ped_griffin$BY)
#' SeqOUT_griffin <- sequoia(GenoS, LH_griffin,
#'                   MaxSibIter = 10,
#'                   args.AP = list(Smooth = FALSE)) }
NULL
