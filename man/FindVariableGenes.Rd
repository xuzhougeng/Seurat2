% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{FindVariableGenes}
\alias{FindVariableGenes}
\title{Identify variable genes}
\usage{
FindVariableGenes(object, mean.function = ExpMean,
  dispersion.function = LogVMR, do.plot = TRUE, set.var.genes = TRUE,
  x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1,
  y.high.cutoff = Inf, num.bin = 20, binning.method = "equal_width",
  selection.method = "mean.var.plot", top.genes = 1000, do.recalc = TRUE,
  sort.results = TRUE, do.cpp = TRUE, display.progress = TRUE, ...)
}
\arguments{
\item{object}{Seurat2 object}

\item{mean.function}{Function to compute x-axis value (average expression). Default
is to take the mean of the detected (i.e. non-zero) values}

\item{dispersion.function}{Function to compute y-axis value (dispersion). Default is to
take the standard deviation of all values/}

\item{do.plot}{Plot the average/dispersion relationship}

\item{set.var.genes}{Set object@var.genes to the identified variable genes
(default is TRUE)}

\item{x.low.cutoff}{Bottom cutoff on x-axis for identifying variable genes}

\item{x.high.cutoff}{Top cutoff on x-axis for identifying variable genes}

\item{y.cutoff}{Bottom cutoff on y-axis for identifying variable genes}

\item{y.high.cutoff}{Top cutoff on y-axis for identifying variable genes}

\item{num.bin}{Total number of bins to use in the scaled analysis (default
is 20)}

\item{binning.method}{Specifies how the bins should be computed. Available methods are:
\itemize{
\item{equal_width:}{ each bin is of equal width along the x-axis [default]}
\item{equal_frequency:}{ each bin contains an equal number of genes (can increase
statistical power to detect overdispersed genes at high expression values, at
the cost of reduced resolution along the x-axis)}
}}

\item{selection.method}{Specifies how to select the genes to store in @var.genes.
\itemize{
\item{mean.var.plot: }{Default method, placing cutoffs on the mean variablility plot}
\item{dispersion: }{Choose the top.genes with the highest dispersion}
}}

\item{top.genes}{Selects the genes with the highest value according to the
selection method.}

\item{do.recalc}{TRUE by default. If FALSE, plots and selects variable genes without recalculating statistics for each gene.}

\item{sort.results}{If TRUE (by default), sort results in object@hvg.info in decreasing order of dispersion}

\item{do.cpp}{Run c++ version of mean.function and dispersion.function if they
exist.}

\item{display.progress}{show progress bar for calculations}

\item{...}{Extra parameters to VariableGenePlot}
}
\value{
Returns a Seurat2 object, placing variable genes in object@var.genes.
The result of all analysis is stored in object@hvg.info
}
\description{
Identifies genes that are outliers on a 'mean variability plot'. First, uses
a function to calculate average expression (mean.function) and dispersion (dispersion.function)
for each gene. Next, divides genes into num.bin (deafult 20) bins based on
their average expression, and calculates z-scores for dispersion within each
bin. The purpose of this is to identify variable genes while controlling for
the strong relationship between variability and average expression.
}
\details{
Exact parameter settings may vary empirically from dataset to dataset, and
based on visual inspection of the plot.
Setting the y.cutoff parameter to 2 identifies genes that are more than two standard
deviations away from the average dispersion within a bin. The default X-axis function
is the mean expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
calculations are not performed in log-space, but the results are reported in log-space -
see relevant functions for exact details.
}
\examples{
pbmc_small <- FindVariableGenes(object = pbmc_small, do.plot = FALSE)
pbmc_small@var.genes

}
\seealso{
\code{VariableGenePlot}
}
