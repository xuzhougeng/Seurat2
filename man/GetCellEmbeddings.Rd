% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimensional_reduction_utilities.R
\name{GetCellEmbeddings}
\alias{GetCellEmbeddings}
\title{Dimensional Reduction Cell Embeddings Accessor Function}
\usage{
GetCellEmbeddings(object, reduction.type = "pca", dims.use = NULL,
  cells.use = NULL)
}
\arguments{
\item{object}{Seurat2 object}

\item{reduction.type}{Type of dimensional reduction to fetch (default is PCA)}

\item{dims.use}{Dimensions to include (default is all stored dims)}

\item{cells.use}{Cells to include (default is all cells)}
}
\value{
Cell embedding matrix for given reduction, cells, and dimensions
}
\description{
Pull cell embeddings matrix for specified stored dimensional reduction
analysis
}
\examples{
pbmc_small
# Examine the head of the first 5 PC cell embeddings
head(GetCellEmbeddings(object = pbmc_small, reduction.type = "pca", dims.use = 1:5))

}
