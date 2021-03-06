% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PedToNum.R
\name{getGDO}
\alias{getGDO}
\title{get individual category: genotyped, dummy, observed}
\usage{
getGDO(id, gID = NULL, DumPrefix = c("F0", "M0"))
}
\arguments{
\item{id}{character vector with ids, may include NA (as in parent columns in
pedigree)}

\item{gID}{character vector with ids of genotyped individuals}

\item{DumPrefix}{prefix that defines dummy individuals}
}
\value{
a factor with length equal to id, and levels Genotyped, Dummy,
  Observed, None
}
\description{
get individual category: genotyped, dummy, observed
}
\details{
In contrast to \code{\link{GetDummifiable}}, this function only
  identifies individuals who currently have a dummy code.
}
\keyword{internal}
