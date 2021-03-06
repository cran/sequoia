% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindFamilies.R
\name{FindFamilies}
\alias{FindFamilies}
\title{Assign Family IDs}
\usage{
FindFamilies(Ped = NULL, SeqList = NULL, UseMaybeRel = FALSE)
}
\arguments{
\item{Ped}{dataframe with columns id - parent1 - parent2; only the
first 3 columns will be used.}

\item{SeqList}{list as returned by \code{\link{sequoia}}. If 'Ped' is not
provided, the element 'Pedigree' from this list will be used if present,
and element 'Pedigreepar' otherwise.}

\item{UseMaybeRel}{use \code{SeqList$MaybeRel}, the dataframe with probable
but non-assigned relatives, to assign additional family IDs?}
}
\value{
A dataframe with the provided pedigree, with a column 'FID' added.
}
\description{
Add a column with family IDs (FIDs) to a pedigree, with each
 number denoting a cluster of connected individuals.
}
\details{
This function repeatedly finds all ancestors and all descendants of
 each individual in turn, and ensures they all have the same Family ID. Not
 all connected individuals are related, e.g. all grandparents of an
 individual will have the same FID, but will typically be unrelated.

When UseMaybeRel = TRUE, probable relatives are added to existing family
 clusters, or existing family clusters may be linked together. Currently no
 additional family clusters are created.
}
