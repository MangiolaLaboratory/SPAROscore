
#' Compute SPAROscores
#'
#' Compute SPAROscores for one or more gene signatures from a variety of
#' supported data containers.
#' `sparoscore()` is an S4 generic with methods for matrix-like objects,
#' Seurat objects, and Bioconductor `SummarizedExperiment`-derived classes.
#' Depending on the input type, the function either returns a matrix of
#' SPAROscores or appends SPAROscores to the input object's metadata.
#'
#'
#' @param data Input data object.
#' Supported classes include:
#' \itemize{
#'   \item Dense matrices (`matrix`)
#'   \item Sparse matrices (`sparseMatrix`)
#'   \item Delayed matrices (`DelayedMatrix`)
#'   \item Data frames
#'   \item Seurat objects
#'   \item `SummarizedExperiment` (and derived objects)
#'   \item `SingleCellExperiment`
#'   \item `SpatialExperiment`
#'   \item `RangedSummarizedExperiment`
#' }
#'
#'
#' @param signatures Gene signature(s) to score.
#' Supported inputs include:
#' \itemize{
#'   \item A character vector representing a single gene signature.
#'   \item A named list of character vectors representing multiple signatures.
#'   \item A `GeneSet` object.
#'   \item A `GeneSetCollection` object.
#' }
#'
#'
#'
#' @param data_has_ranks Logical indicating whether the supplied data already
#' contains feature ranks. If `FALSE` (default), ranks are computed from
#' expression values prior to scoring.
#'
#'
#'
#' @param assay Character string specifying the name of
#' the assay containing the data to use.
#' When data_has_ranks = FALSE, this assay should contain expression counts.
#' When data_has_ranks = TRUE, it should contain pre-computed ranks.
#' Defaults to "RNA" for Seurat objects.
#' Defaults to "counts" for SummarizedExperiment-derived objects.
#' Ignored for matrix-like inputs.
#'
#'
#' @param layer Character string specifying the name
#' of the assay layer containing the data to use.
#' When data_has_ranks = FALSE, this layer should contain expression counts.
#' When data_has_ranks = TRUE, it should contain pre-computed ranks.
#' Defaults to "counts"
#' Used only for Seurat objects.
#'
#'
#' @param count_caps Optional numeric vector containing expression thresholds
#' values for each column in the counts data.
#' These values are used to determine the corresponding rank caps.
#' Only used when data_has_ranks = FALSE.
#' If NULL (default), values are computed using
#' \code{\link{compute_geometric_average}}.
#'
#'
#' @param rank_caps An optional named numeric vector containing the rank cap
#' associated with each column of ranks data used during
#' SPAROscore calculation. Typically obtained from the
#' rank_caps component returned by .
#' When data_has_ranks = TRUE and rank_cpas is NULL (default), rank caps are
#' obtained from \code{\link{get_ranks}} for matrix-like objects,
#' retrieved from the "rank_caps" column of metadata for Seurat and
#' SummarizedExperiment objects
#'
#' When both count_caps and rank_cpas are provided, rank_caps takes precedence.
#'
#'
#' @param handle_ties Character string specifying how tied expression values
#' should be ranked. Passed directly to
#' MatrixGenerics::colRanks(ties.method = ...).
#'
#' Supported values are:
#' \describe{
#' \item{"min"}{Assign the minimum rank to tied values (default).}
#' \item{"max"}{Assign the maximum rank to tied values.}
#' \item{"average"}{Assign the average rank to tied values.}
#' \item{"random"}{Break ties at random.}
#' }
#'
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from ranks should be handled.
#'
#' Supported options are:
#' \describe{
#' \item{"skip"}{
#' Exclude missing genes from score calculation (default).
#' }
#'
#' \item{"impute"}{
#' Include missing genes by assigning them the capped rank value.
#' }
#' }
#'
#' @param down_signatures Gene signature(s) to be considered for scoring
#'  the down-regulation effect. Supported inputs are same as signatures.
#'  If provided, names(down_signatures) must match names(signatures).
#'  Defaults to NULL.
#'
#' When down_signatures is not NULL,
#' Final Score <- Score(signatures) - Score(down_signatures)
#'
#'
#' @param prefix Character string to be appended before the headers of the
#' columns returning scores. Defaults to "".
#'
#'
#' @return
#' Method-dependent:
#' \itemize{
#'   \item Matrix-like inputs return a numeric matrix of SPAROscores.
#'   \item Seurat inputs return the Seurat object with
#'   \itemize{
#'              \item One or more metadata columns containing SPAROscores
#'              for the supplied signature(s).
#'              \item A metadata column named "rank_caps" containing the
#'              rank cap for each cell.
#'              \item A "ranks" layer added to the selected assay
#'              when ranks are calculated from expression data.
#'          }
#'   \item SummarizedExperiment-derived inputs return the same object with
#'   \itemize{
#'           \item One or more SPAROscore columns added to colData().
#'           \item A "rank_caps" column added to colData().
#'           \item A "ranks" assay added when ranks are computed from
#'           expression data.
#'          }
#' }
#'
#' @seealso
#' \code{\link{get_ranks}},
#' \code{\link{get_scores}},
#'
#' @export
#' @details
#' sparoscore() is an S4 generic that dispatches on the class of data.
#' The method selected determines how SPAROscores are computed and where the
#' results are stored.
#'
#' For matrix-like inputs (matrix, sparseMatrix, DelayedMatrix, and
#' data.frame), the function returns a numeric matrix of SPAROscores.
#'
#' For Seurat objects, SPAROscores and rank caps are appended to the object's
#' metadata. If ranks are computed from expression data, they are also stored
#' in a "ranks" assay layer for future reuse.
#'
#' For SummarizedExperiment-derived objects, SPAROscores and rank caps are
#' appended to colData(). If ranks are computed from expression data, they
#' are also stored in a "ranks" assay.
#'
#' When data_has_ranks = FALSE, feature ranks are first computed using
#' get_ranks(). When data_has_ranks = TRUE, the supplied ranks are used
#' directly and rank computation is skipped.
#'
#'
#' The supplied signature genes are first matched against the genes available
#' in ranks. Missing genes can either be excluded from scoring ("skip") or
#' assigned the capped rank value and included in score calculation
#' ("impute").
#'
#' For each sample, cell, or spatial domain, SPAROscores
#' evaluates the normalized Spearman footrule distance between observed
#' signature gene ranks and the column-specific rank cap.
#'
#' Rank caps typically correspond to the rank of the geometric mean expression
#' value estimated by \code{\link{get_ranks}}, although custom rank caps
#' may also be supplied.
#'
#'
#' Scores typically range from 0 to 1, where high scores indicate high
#' signature gene expression, and low scores indicate signature expression
#' closer to the cap value.
#'
#' @examples
#' # ------------------------------------------------------------------
#' # Matrix-like input
#' # ------------------------------------------------------------------
#'
#' counts <- matrix(
#'   sample(0:10, 500, replace = TRUE),
#'   nrow = 50,
#'   dimnames = list(
#'     paste0("gene", 1:50),
#'     paste0("cell", 1:10)
#'   )
#' )
#'
#' scores <- sparoscore(
#'   data = counts,
#'   signatures = c("gene1", "gene2", "gene3")
#' )
#'
#' # Multiple signatures
#' signatures <- list(
#'   SignatureA = c("gene1", "gene2", "gene3"),
#'   SignatureB = c("gene10", "gene11", "gene12")
#' )
#'
#' scores <- sparoscore(
#'   data = counts,
#'   signatures = signatures
#' )
#'
#' # ------------------------------------------------------------------
#' # Seurat object
#' # ------------------------------------------------------------------
#'
#' seurat_object <- Seurat::CreateSeuratObject(counts = counts)
#'
#' seurat_object <- sparoscore(
#'   data = seurat_object,
#'   signatures = c("gene1", "gene2", "gene3")
#' )
#'
#' # use custom count caps to measure distance from zero expression
#' zero_counts <- setNames(rep(0, ncol(seurat_object)), colnames(seurat_object))
#'
#' seurat_object <- sparoscore(
#'   data = seurat_object,
#'   signatures = c("gene1", "gene2", "gene3"),
#'   count_caps = zero_counts
#' )
#'
#' # use custom rank caps to consider top 5% genes
#' custom_ranks <- setNames(
#'   rep(0.05*nrow(seurat_object), ncol(seurat_object)),
#'   colnames(seurat_object)
#' )
#'
#' seurat_object <- sparoscore(
#'   data = seurat_object,
#'   signatures = c("gene1", "gene2", "gene3"),
#'   rank_caps = custom_ranks
#' )
#'
#' # Reuse previously computed ranks
#' seurat_object <- sparoscore(
#'   data = seurat_object,
#'   signatures = c("gene1", "gene2", "gene3"),
#'   data_has_ranks = TRUE,
#'   layer = "ranks"
#' )
#'
#' # ------------------------------------------------------------------
#' # SingleCellExperiment / SummarizedExperiment
#' # ------------------------------------------------------------------
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = counts)
#' )
#'
#' sce <- sparoscore(
#'   data = sce,
#'   signatures = c("gene1", "gene2", "gene3")
#' )
#'
#' # Reuse previously computed ranks
#' sce <- sparoscore(
#'   data = sce,
#'   signatures = c("gene1", "gene2", "gene3"),
#'   data_has_ranks = TRUE,
#'   assay = "ranks"
#' )



setGeneric("sparoscore",
           function(data,
                    signatures,
                    down_signatures = NULL,
                    data_has_ranks = FALSE,
                    assay = "RNA",
                    layer = "counts",
                    count_caps = NULL,
                    rank_caps = NULL,
                    handle_ties = "min",
                    handle_missing_genes = "skip",
                    prefix = "")
               standardGeneric("sparoscore"))


# set the method for sparoscore() where data is a matrix
#' Compute SPAROscores for matrix objects
#' @rdname sparoscore
#' @export
setMethod("sparoscore",
          signature(data = 'matrix'),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = NULL,
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip",
                   prefix = ""){

              #call the helper function
              sparoscores <- augment_sparoscores_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix)

              return(sparoscores)
          }
)


# set the method for sparoscore() where data is a sparseMatrix
#' Compute SPAROscores for sparseMatrix objects
#' @rdname sparoscore
#' @export
setMethod("sparoscore",
          signature(data = 'sparseMatrix'),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = NULL,
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip",
                   prefix = ""){

              #call the helper function
              sparoscores <- augment_sparoscores_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix)

              return(sparoscores)
          }
)

# set the method for sparoscore() where data is a DelayedMatrix
#' Compute SPARO scores for DelayedMatrix objects
#' @rdname sparoscore
#' @export
setMethod("sparoscore",
          signature(data = 'DelayedMatrix'),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = NULL,
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip",
                   prefix = ""){

              #call the helper function
              sparoscores <- augment_sparoscores_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix)

              return(sparoscores)
          }
)


# set the method for sparoscore() where data is a data.frame
#' Compute SPAROscores for data.frame objects
#' @rdname sparoscore
#' @export
setMethod("sparoscore",
          signature(data = 'data.frame'),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = NULL,
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip",
                   prefix = ""){

              #call the helper function
              sparoscores <- augment_sparoscores_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix)

              return(sparoscores)
          }
)


# set the method for sparoscore() where data is a Seurat object
#' Compute SPAROscores for Seurat objects
#' @rdname sparoscore
#' @export
setMethod("sparoscore",
          signature(data = 'Seurat'),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = "RNA",
                   layer = "counts",
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip",
                   prefix = ""){

              #call the helper function
              seurat_object <- augment_sparoscores_seurat(
                  seurat_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  assay = assay,
                  layer = layer,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix)

              return(seurat_object)
          }
)


# # set a union of class for bioconductor S4 objects
# setClassUnion("BioConductor_Objects",
#               c("SummarizedExperiment",
#                 "SingleCellExperiment",
#                 "SpatialExperiment",
#                 "RangedSummarizedExperiment"
#                 ))

# set the method for sparoscore() where data is one of the
# SummarizedExperiment-derived classes

#' Compute SPAROscores for SummarizedExperiment-like objects
#' @rdname sparoscore
#' @export
setMethod("sparoscore",
          signature(data = 'SummarizedExperiment'),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = "counts",
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip",
                   prefix = ""){

              #call the helper function
              sce_object <- augment_sparoscores_sce(
                  sce_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  assay = assay,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix)

              return(sce_object)
          }
)

