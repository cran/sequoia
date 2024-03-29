% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimGeno.R
\name{DoErrors}
\alias{DoErrors}
\title{Fortran Simulate Genotyping Errors}
\usage{
DoErrors(SGeno, Act2Obs)
}
\arguments{
\item{SGeno}{matrix with genotype data, size nInd x nSnp.}

\item{Act2Obs}{array with conditional probability of observing genotype i
conditional on actual genotype j, size nSnp x 3 x 3.}
}
\value{
\code{SGeno} with errors.
}
\description{
Wrapper for Fortran function to simulate genotyping errors.
}
\keyword{internal}
