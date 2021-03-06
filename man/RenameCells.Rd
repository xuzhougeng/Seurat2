% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interaction.R
\name{RenameCells}
\alias{RenameCells}
\title{Rename cells}
\usage{
RenameCells(object, add.cell.id = NULL, new.names = NULL,
  for.merge = FALSE)
}
\arguments{
\item{object}{Seurat2 object}

\item{add.cell.id}{prefix to add cell names}

\item{new.names}{vector of new cell names}

\item{for.merge}{Only rename slots needed for merging Seurat2 objects.
Currently only renames the raw.data and meta.data slots.}
}
\value{
Seurat2 object with new cell names
}
\description{
Change the cell names in all the different parts of a Seurat2 object. Can
be useful before combining multiple objects.
}
\details{
If \code{add.cell.id} is set a prefix is added to existing cell names. If
\code{new.names} is set these will be used to replace existing names.
}
\examples{
head(pbmc_small@cell.names)
pbmc_small <- RenameCells(pbmc_small, add.cell.id = "Test")
head(pbmc_small@cell.names)

}
