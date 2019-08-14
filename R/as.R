#' @include conversion.R
NULL

#' @rdname Convert
#' @export as.seurat2
#' @aliases as.seurat2
#'
as.seurat2 <- function(from) {
  UseMethod(generic = 'as.seurat2', object = from)
}

#' @rdname Convert
#' @export
#' @method as.seurat2 SingleCellExperiment
#'
as.seurat2.SingleCellExperiment <- function(from) {
  return(Convert(from = from, to = 'seurat2'))
}

#' @rdname Convert
#' @export as.SingleCellExperiment
#' @aliases as.SingleCellExperiment
#'
as.SingleCellExperiment <- function(from) {
  UseMethod(generic = 'as.SingleCellExperiment', object = from)
}

#' @rdname Convert
#' @export
#' @method as.SingleCellExperiment seurat2
#'
as.SingleCellExperiment.seurat2 <- function(from) {
  return(Convert(from = from, to = 'sce'))
}
