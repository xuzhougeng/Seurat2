% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_determination.R
\name{SaveClusters}
\alias{SaveClusters}
\title{Save cluster assignments to a TSV file}
\usage{
SaveClusters(object, file)
}
\arguments{
\item{object}{Seurat2 object with cluster assignments}

\item{file}{Path to file to write cluster assignments to}
}
\value{
No return value. Writes clusters assignments to specified file.
}
\description{
Save cluster assignments to a TSV file
}
\examples{
\dontrun{
pbmc_small
file.loc <- "~/Desktop/cluster_assignments.tsv"
SaveClusters(object = pbmc_small, file = file.loc)
}

}
