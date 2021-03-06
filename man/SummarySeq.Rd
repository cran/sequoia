% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeqListSummary.R
\name{SummarySeq}
\alias{SummarySeq}
\title{Summarise Sequoia Output or Pedigree}
\usage{
SummarySeq(
  SeqList = NULL,
  Pedigree = NULL,
  DumPrefix = c("F0", "M0"),
  SNPd = NULL,
  Plot = TRUE,
  Panels = "all"
)
}
\arguments{
\item{SeqList}{the list returned by \code{\link{sequoia}}. Only elements
'Pedigree' or 'PedigreePar' and 'AgePriors' are used.}

\item{Pedigree}{dataframe, pedigree with the first three columns being id - dam -
sire. Column names are ignored, as are additional columns.}

\item{DumPrefix}{character vector of length 2 with prefixes for dummy dams
(mothers) and sires (fathers). Will be read from \code{SeqList}'s 'Specs'
if provided. Used to distinguish between dummies and non-dummies. Length 3
in case of hermaphrodites.}

\item{SNPd}{character vector with ids of SNP genotyped individuals. Only when
\code{Pedigree} is provided instead of \code{SeqList}, then used to distinguish
between genetically assigned parents and 'observed' parents (e.g. observed
in the field, or assigned previously using microsatellites). If
\code{SeqList}'s 'PedigreePar' is provided, all ids in that dataframe will
be presumed genotyped.}

\item{Plot}{show barplots and histograms of the results, as well as of the
parental LLRs, Mendelian errors, and agepriors, if present.}

\item{Panels}{character vector with panel(s) to plot. Choose from 'all',
'G.parents' (parents of genotyped individuals), 'D.parents' (parents of
dummy individuals), 'sibships' (distribution of sibship sizes), 'LLR'
(log10-likelihood ratio parent/otherwise related), 'OH' (count of opposite
homozygote SNPs).}
}
\value{
A list with the following elements:
  \item{PedSummary}{a 2-column matrix with basic summary statistics, similar
  to what used to be returned by \pkg{Pedantics}' \code{pedStatSummary} (now
  archived on CRAN). First column refers to the complete pedigree, second
  column to SNP-genotyped individuals only. Maternal siblings sharing a dummy
  parent are counted in the 2nd column if both sibs are genotyped, but not if
  one of the sibs is a dummy individual.}
  \item{ParentCount}{a 2x3x2x4 array with the number of assigned parents,
  split by D1: genotyped vs dummy individuals; D2: female, male and
  unknown-sex individuals; D3: dams vs sires; D4: genotyped, dummy, observed
  vs no parent}
  \item{GPCount}{a 2x4x4 array with the number of assigned grandparents,
  split by D1: genotyped vs dummy individuals; D2 Maternal grandmother (MGM),
  maternal grandfather (MGF), paternal grandmother (PGM), paternal
  grandfather (PGF); D3: genotyped, dummy, observed vs no grandparent}
  \item{SibSize}{a list with as first element a table of maternal sibship
  sizes, and as second element a table of paternal sibship sizes. Each table
  is a matrix with a number of rows equal to the  maximum sibship size, and 3
  columns, splitting by the type of parent: genotyped, dummy, or observed.}
}
\description{
Number of assigned parents and grandparents and sibship sizes,
  split by genotyped, dummy, and 'observed'.
}
\examples{
\donttest{
data(SimGeno_example, LH_HSg5, package="sequoia")
SeqOUT <- sequoia(GenoM = SimGeno_example,
                   LifeHistData = LH_HSg5, Module="ped")
Ped_example <- SeqOUT[["Pedigree"]]
Ped_example$dam[1:20] <- paste0("Mum", 1:20)  # some fake field mums
PedSum1 <- SummarySeq(SeqOUT, Pedigree=Ped_example, Panels="G.parents")
summary(PedSum1)
PedSum1$PedSummary
}

}
\seealso{
\code{\link{sequoia}} for pedigree reconstruction and links to other
  functions.
}
