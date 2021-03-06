% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interaction.R
\name{FetchData}
\alias{FetchData}
\title{Access cellular data}
\usage{
FetchData(object, vars.all = NULL, cells.use = NULL, use.imputed = FALSE,
  use.scaled = FALSE, use.raw = FALSE)
}
\arguments{
\item{object}{Seurat2 object}

\item{vars.all}{List of all variables to fetch}

\item{cells.use}{Cells to collect data for (default is all cells)}

\item{use.imputed}{For gene expression, use imputed values. Default is FALSE}

\item{use.scaled}{For gene expression, use scaled values. Default is FALSE}

\item{use.raw}{For gene expression, use raw values. Default is FALSE}
}
\value{
A data frame with cells as rows and cellular data as columns
}
\description{
Retreives data (gene expression, PCA scores, etc, metrics, etc.) for a set
of cells in a Seurat2 object
}
\examples{
pc1 <- FetchData(object = pbmc_small, vars.all = 'PC1')
head(x = pc1)

}
