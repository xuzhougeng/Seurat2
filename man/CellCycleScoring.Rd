% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoring.R
\name{CellCycleScoring}
\alias{CellCycleScoring}
\title{Score cell cycle phases}
\usage{
CellCycleScoring(object, g2m.genes, s.genes, set.ident = FALSE)
}
\arguments{
\item{object}{A Seurat2 object}

\item{g2m.genes}{A vector of genes associated with G2M phase}

\item{s.genes}{A vector of genes associated with S phases}

\item{set.ident}{If true, sets identity to phase assignments
Stashes old identities in 'old.ident'}
}
\value{
A Seurat2 object with the following columns added to object@meta.data: S.Score, G2M.Score, and Phase
}
\description{
Score cell cycle phases
}
\examples{
\dontrun{
# pbmc_small doesn't have any cell-cycle genes
# To run CellCycleScoring, please use a dataset with cell-cycle genes
# An example is available at http://satijalab.org/seurat2/cell_cycle_vignette.html
pbmc_small <- CellCycleScoring(
  object = pbmc_small,
  g2m.genes = cc.genes$g2m.genes,
  s.genes = cc.genes$s.genes
)
head(x = pbmc_small@meta.data)
}

}
\seealso{
\code{AddModuleScore}
}
