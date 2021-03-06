% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimensional_reduction_utilities.R
\name{PrintDim}
\alias{PrintDim}
\title{Print the results of a dimensional reduction analysis}
\usage{
PrintDim(object, reduction.type = "pca", dims.print = 1:5,
  genes.print = 30, use.full = FALSE)
}
\arguments{
\item{object}{Seurat2 object}

\item{reduction.type}{Reduction technique to print results for}

\item{dims.print}{Number of dimensions to display}

\item{genes.print}{Number of genes to display}

\item{use.full}{Use full PCA (i.e. the projected PCA, by default FALSE)}
}
\value{
Set of genes defining the components
}
\description{
Prints a set of genes that most strongly define a set of components
}
\examples{
pbmc_small
PrintDim(object = pbmc_small, reduction.type = "pca")
# Options for how many dimensions and how many genes to print
PrintDim(object = pbmc_small, reduction.type = "pca", dims.print = 1:2, genes.print = 5)
# Can also print for the projected PCA
PrintDim(object = pbmc_small, reduction.type = "pca", use.full = TRUE)

}
