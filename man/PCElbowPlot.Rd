% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{PCElbowPlot}
\alias{PCElbowPlot}
\title{Quickly Pick Relevant PCs}
\usage{
PCElbowPlot(object, num.pc = 20)
}
\arguments{
\item{object}{Seurat2 object}

\item{num.pc}{Number of PCs to plot}
}
\value{
Returns ggplot object
}
\description{
Plots the standard deviations (or approximate singular values if running PCAFast)
of the principle components for easy identification of an elbow in the graph.
This elbow often corresponds well with the significant PCs and is much faster to run.
}
\examples{
PCElbowPlot(object = pbmc_small)

}
