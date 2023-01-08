#' @name Ped_HSg5
#' @docType data
#' @title Example pedigree: 'HSg5'
#'
#' @description A pedigree with five non-overlapping generations and considerable
#' inbreeding. Each female mated with two random males and each male with three
#' random females, producing four full-sib offspring per mating. This is
#' \strong{Pedigree II} in the paper.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @references Huisman, J. (2017) Pedigree reconstruction from SNP data:
#'   Parentage assignment, sibship clustering, and beyond. Molecular Ecology
#'   Resources 17:1009--1024.
#'
#' @seealso \code{\link{LH_HSg5} \link{SimGeno_example} \link{sequoia}}
#'
#' @keywords datasets sequoia
#' @usage data(Ped_HSg5)
#' @format A data frame with 1000 rows and 3 variables (id, dam, sire)
NULL


#===============================================================================
#' @name LH_HSg5
#' @docType data
#' @title Example life history file: 'HSg5'
#'
#' @description This is the life history file associated with
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
#' @keywords datasets sequoia
#' @usage data(LH_HSg5)
#' @format A data frame with 1000 rows and 3 variables:
#' \describe{
#' \item{ID}{Female IDs start with 'a', males with 'b'; the next 2 numbers give
#'   the generation number (00 -- 05), the last 3 numbers the individual ID
#'   number (runs continuously across all generations)}
#' \item{Sex}{1 = female, 2 = male}
#' \item{BirthYear}{from 2000 (generation 0, founders) to 2005} }
#'
NULL


#===============================================================================
#' @name SimGeno_example
#' @docType data
#' @title Example genotype file: 'HSg5'
#'
#' @description Simulated genotype data for cohorts 1+2 in Pedigree
#'   \code{\link{Ped_HSg5}}
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{LH_HSg5}, \link{SimGeno}}
#'
#' @keywords datasets sequoia
#' @usage data(SimGeno_example)
#' @format A genotype matrix with 214 rows (ids) and 200 columns (SNPs). Each
#'   SNP is coded as 0/1/2 copies of the reference allele, with -9 for missing
#'   values. Ids are stored as rownames.
NULL


#===============================================================================
#' @name Geno_HSg5
#' @docType data
#' @title Example genotype file: 'HSg5'
#'
#' @description Simulated genotype data for all* individuals in Pedigree
#'   \code{\link{Ped_HSg5}} (*: with 40% of parents are non-genotyped).
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{LH_HSg5}, \link{SimGeno}, \link{SeqOUT_HSg5}}
#'
#' @keywords datasets sequoia
#' @usage data(Geno_HSg5)
#' @format A genotype matrix with 920 rows (ids) and 200 columns (SNPs). Each
#'   SNP is coded as 0/1/2 copies of the reference allele, with -9 for missing
#'   values. Ids are stored as rownames.
#'
#' @examples
#' \dontrun{
#' # this output was created as follows:
#' Geno_HSg5 <- SimGeno(Ped = Ped_HSg5, nSnp = 200, ParMis=0.4,
#'                      CallRate = 0.9, SnpError = 0.005)
#' }
NULL


#===============================================================================
#' @name SeqOUT_HSg5
#' @docType data
#' @title Example output from pedigree inference: 'HSg5'
#'
#' @description Example output of a \code{\link{sequoia}} run including sibship
#'   clustering, based on Pedigree \code{\link{Geno_HSg5}}.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{Ped_HSg5}, \link{LH_HSg5}}
#'
#' @keywords datasets sequoia
#' @usage data(SeqOUT_HSg5)
#' @format a list, see \code{\link{sequoia}}
#'
#' @examples
#' \dontrun{
#' # this output was created as follows:
#' Geno <- SimGeno(Ped = Ped_HSg5, nSnp = 200)
#' SeqOUT_HSg5 <- sequoia(GenoM = Geno, LifeHistData = LH_HSg5, Module = "ped",
#'                        Err = 0.005)
#' }
#' # some ways to inspect the output; see vignette for more info:
#' names(SeqOUT_HSg5)
#' SeqOUT_HSg5$Specs
#' SummarySeq(SeqOUT_HSg5)
NULL


#===============================================================================
#===============================================================================
#' @name Ped_griffin
#' @docType data
#' @title Example pedigree: griffins
#'
#' @description Example pedigree with overlapping generations and polygamy.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso  \code{\link{LH_griffin}}; \code{\link{SeqOUT_griffin}} for a
#'   sequoia run on simulated genotype data based on this pedigree;
#'   \code{\link{Ped_HSg5}} for another pedigree; \code{\link{sequoia}}.
#'
#' @keywords datasets sequoia
#' @usage data(Ped_griffin)
#' @format A data frame with 200 rows and 4 variables (id, dam, sire, birthyear)
#'
#' @section Code:
#' The R code used to create this pedigree can be found in /data-raw.
NULL


#===============================================================================
#' @name LH_griffin
#' @docType data
#' @title Example life history data: griffins
#'
#' @description Example life history data associated with the griffin pedigree.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{Ped_griffin}}, \code{\link{SeqOUT_griffin}}
#'
#' @keywords datasets sequoia
#' @usage data(LH_griffin)
#' @format A data frame with 200 rows and 3 variables (ID, Sex, BirthYear)
#'
#' @examples
#' \dontrun{
#' BY <- rep(c(2001:2010), each=20)
#' Sex <- sample.int(n=2, size=200, replace=TRUE)
#' ID <- paste0("i", formatC(1:200, width=3, flag="0"), "_", BY, "_",
#'              ifelse(Sex==1, "F", "M"))
#' LH_griffin <- data.frame(ID, Sex, BirthYear = BY)
#' }
NULL


#===============================================================================
#' @name Geno_griffin
#' @docType data
#' @title Example genotype file: Griffins
#'
#' @description Simulated genotype data from Pedigree \code{\link{Ped_griffin}}
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{SimGeno}}
#'
#' @keywords datasets sequoia
#' @usage data(Geno_griffin)
#' @format A genotype matrix with 142 rows (individuals) and 200 columns (SNPs).
#'   Each SNP is coded as 0/1/2 copies of the reference allele, with -9 for
#'   missing values. Ids are stored as rownames.
NULL


#===============================================================================
#' @name SeqOUT_griffin
#' @docType data
#' @title Example output from pedigree inference: griffins
#'
#' @description Example output of a sequoia run including sibship clustering,
#'   with \code{\link{Geno_griffin}} as input (simulated from
#'   \code{\link{Ped_griffin}}).
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{sequoia}}
#'
#' @keywords datasets sequoia
#' @usage data(SeqOUT_griffin)
#' @format a list, see \code{\link{sequoia}}
#'
#' @examples
#' \dontrun{
#' SeqOUT_griffin <- sequoia(GenoM = Geno_griffin,
#'                           LifeHistData = LH_griffin,
#'                           Module = 'ped')
#' }
NULL


#===============================================================================
#' @name Conf_griffin
#' @docType data
#' @title Example output from estimating confidence probabilities: griffins
#'
#' @description Example output of \code{\link{EstConf}}, with the inferred
#'   pedigree in \code{\link{SeqOUT_griffin}} used as reference pedigree.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{Ped_griffin}}, \code{\link{Geno_griffin}},
#'
#' @keywords datasets sequoia
#' @usage data(Conf_griffin)
#' @format a list, see \code{\link{sequoia}}
#'
#' @examples
#' \dontrun{
#' Conf_griffin <- EstConf(Pedigree = SeqOUT_griffin$Pedigree,
#'                         LifeHistData = LH_griffin,
#'                         args.sim = list(nSnp = 400, SnpError = 0.001,
#'                                         ParMis=0.4),
#'                         args.seq = list(Module = 'ped', Err=0.001),
#'                         nSim = 20,
#'                         nCores = 5,
#'                         quiet = TRUE)
#' }
NULL


#===============================================================================
#' @name MaybeRel_griffin
#' @docType data
#' @title Example output from check for relatives: griffins
#'
#' @description Example output of a check for parent-offspring pairs and
#'   parent-parent-offspring trios with \code{\link{GetMaybeRel}}, with
#'   \code{\link{Geno_griffin}} as input (simulated from
#'   \code{\link{Ped_griffin}}).
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{SeqOUT_griffin}}
#'
#' @keywords datasets sequoia
#' @usage data(MaybeRel_griffin)
#' @format a list with 2 dataframes, 'MaybePar' and 'MaybeTrio'. See
#'   \code{\link{GetMaybeRel}} for further details.
#'
#' @examples
#' \dontrun{
#' MaybeRel_griffin <- GetMaybeRel(GenoM = Geno_griffin, Err=0.001,
#'                                 Module = 'par')
#' }
NULL


#===============================================================================
#' @name FieldMums_griffin
#' @docType data
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
#' @keywords datasets sequoia
#' @usage data(FieldMums_griffin)
#' @format A data frame with 144 rows and 2 variables (id, mum)
#' @examples
#' \dontrun{
#' PC_griffin <- PedCompare(Ped1 = cbind(FieldMums_griffin, sire=NA),
#'                          Ped2 = SeqOUT_griffin$Pedigree)
#' }
NULL



#===============================================================================
#===============================================================================
#' @name Inherit_patterns
#' @docType data
#' @title Inheritance patterns
#'
#' @description Inheritance patterns used by SimGeno for non-autosomal SNPs,
#'   identical to those in Inherit.xlsx
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @seealso \code{\link{SimGeno}}
#'
#' @keywords datasets sequoia inherit
#' @usage data(Inherit_patterns)
#' @format An array with the following dimensions:
#' \describe{
#'   \item{d1}{type: autosomal, x-chromosome, y-chromosome, or mtDNA}
#'   \item{d2}{offspring sex: female, male, or unknown}
#'   \item{d3}{offspring genotype: aa (0), aA (1), Aa (1), or AA (2)}
#'   \item{d4}{mother genotype}
#'   \item{d5}{father genotype}
#' }
NULL
