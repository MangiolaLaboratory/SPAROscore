#' Validate a gene signature against available genes
#'
#' Checks whether the genes in a user-supplied signature are present in a
#' reference set of genes with available pre-computed ranks. The function
#' returns the subset of signature genes that can be used for scoring and the
#' subset that are missing.
#'
#' A warning is issued if one or more signature genes are not found in
#' `all_genes`. An error is raised if none of the supplied signature genes are
#' present.
#'
#' @noRd
#'
#'
#' @param signature Character vector containing the genes that define the
#'   signature to be scored.
#'
#'
#' @param all_genes Character vector of genes for which pre-computed ranks are
#'   available.
#'
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{valid_genes}{Character vector of signature genes found in
#'   `all_genes` and therefore available for scoring.}
#'   \item{missing_genes}{Character vector of signature genes not found in
#'   `all_genes`.}
#' }
#'
#'
#' @details
#' The order of genes in `valid_genes` matches their order in the input
#' `signature`. Missing genes are reported in a warning message. If no
#' signature genes are found in `all_genes`, the function stops with an error.
#'
#'
#' Example of how to use the function
#' validate_signature(
#'   signature = c("CCR7", "CD62L"),
#'   all_genes = c("CCR7", "CD27", "CD28")
#' )
#'
#' validate_signature(
#'   signature = gene_signature,
#'   all_genes = rownames(ranks)
#' )
#'

validate_signature <- function(signature, all_genes){

    valid_signature <- signature[signature %in% all_genes]

    invalid_genes <- signature[!(signature %in% all_genes)]

    if(length(invalid_genes) > 0){
        warning(paste0("SPAROscore says: The following ", length(invalid_genes),
                       " signature genes are missing in the input dataset: ",
                       paste0(invalid_genes, collapse = ", ")), call. = FALSE)
    }

    if(length(valid_signature) == 0){
        stop("SPAROscore says:
             No signature genes are present in the input dataset")
    }

    return(list(valid_genes = valid_signature, missing_genes = invalid_genes))
}




#' Impute ranks for missing signature genes
#'
#' Replaces missing signature genes by assigning each of them a rank equal to
#' `rank_cap`. The imputed ranks are appended to the supplied vector of
#' available gene ranks, producing a complete rank vector for downstream
#' signature scoring.
#'
#' @noRd
#'
#'
#' @param incomplete_ranks Numeric vector containing the ranks of signature
#'   genes that are present in a cell, sample, or spatial location.
#'
#' @param rank_cap Numeric value used to impute ranks for missing genes.
#'   Typically represents the maximum rank considered during scoring.
#'
#'
#' @param missing_genes Character vector of signature genes that are absent
#'   from the input data and therefore require rank imputation.
#'
#'
#' @return A numeric vector containing the original ranks in
#'   `incomplete_ranks` followed by imputed ranks for each gene in
#'   `missing_genes`.
#'
#'
#' @details
#' One value equal to `rank_cap` is added for each gene listed in
#' `missing_genes`. If `missing_genes` is empty, the input vector is returned
#' unchanged.
#'
#'
#' Example of how to use the function
#' impute_missing_gene_ranks(
#'   incomplete_ranks = c(1, 2, 3, 4, 5, 6, 7, 8),
#'   rank_cap = 5,
#'   missing_genes = c("CD27", "CD28")
#' )
#'
#'
impute_missing_gene_ranks <- function(incomplete_ranks,rank_cap, missing_genes){

    # if no genes are missing, use incomplete ranks
    if (is.null(missing_genes) || length(missing_genes) == 0) {
        return(incomplete_ranks)
    }

    # Create the additional rows of imputed ranks as a separate matrix
    imputed_ranks <- rep(rank_cap, length(missing_genes))

    # Append the imputed matrix to in incomplete gene ranks
    full_gene_ranks <- c(incomplete_ranks, imputed_ranks)

    return(full_gene_ranks)
}


#' Append a numeric vector as a new row to a matrix-like object
#'
#' Appends a numeric vector to a matrix-like object while preserving the
#' original object class. Supported input types include base R matrices,
#' Matrix::sparseMatrix objects, and DelayedArray::DelayedMatrix objects.
#'
#' The supplied vector is converted to the appropriate one-row representation
#' before being combined with the input object.
#'
#' @noRd
#'
#'
#' @param matrix_like_object A matrix-like object to which a new row will be
#' appended. Supported classes are base matrices, sparse matrices from the
#' Matrix package, and delayed matrices from the DelayedArray package.
#' @param numeric_vector Numeric vector to append as a new row. Its length must
#' equal the number of columns in matrix_like_object.
#'
#'
#' @return An object of the same class as matrix_like_object, with
#' numeric_vector appended as an additional row.
#'
#'
#' @details
#' The function validates that the length of numeric_vector matches the
#' number of columns in matrix_like_object. An error is raised if the
#' dimensions are incompatible or if the supplied object class is not
#' supported.
#'
#'
#' Example of how to use the function
#' # Base matrix
#' mat <- matrix(1:6, nrow = 2)
#' append_to_matrix_like_object(mat, c(7, 8, 9))
#'
#' # Sparse matrix
#' sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
#' append_to_matrix_like_object(sparse_mat, c(7, 8, 9))
#'
#' # Delayed matrix
#' delayed_mat <- DelayedArray::DelayedArray(mat)
#' append_to_matrix_like_object(delayed_mat, c(7, 8, 9))
#'
#'
append_to_matrix_like_object <- function(matrix_like_object, numeric_vector){
    if(ncol(matrix_like_object) != length(numeric_vector)){
        stop("SPAROscore says:
             Cannot append cap values. Length not equal to columns")
    }

    if(inherits(matrix_like_object, "DelayedMatrix")){
        #convert numeric vector to delayedmatrix
        delayed_vector <- DelayedArray::DelayedArray(
            matrix(numeric_vector, nrow = 1))
        return(DelayedArray::rbind(matrix_like_object, delayed_vector))
    }
    else if(inherits(matrix_like_object, "sparseMatrix")){
        #convert numeric vector to sparsematrix
        sparse_vector <- Matrix::Matrix(numeric_vector, nrow = 1, sparse = TRUE)
        return(rbind(matrix_like_object, sparse_vector))
    }
    else if(is.matrix(matrix_like_object)){
        return(rbind(matrix_like_object, numeric_vector))
    }
    else{
        stop("SPAROscore says: Unsupported file type for counts/caps")
    }
}


#' Compute column-wise geometric mean expression values
#'
#' Computes the geometric mean expression for each sample, cell, or spatial
#' location in a count matrix. These values are used to derive rank caps for
#' downstream signature scoring.
#'
#' The geometric mean is calculated as:
#' \deqn{\exp(\mathrm{mean}(\log(1 + x)))}
#' where \eqn{x} represents the expression values in a column.
#'
#' @keywords internal
#'
#'
#' @param counts A matrix-like object containing expression counts, with
#' genes in rows and samples, cells, or spatial locations in columns.
#' Supported inputs include base matrices, sparse matrices
#' (Matrix::sparseMatrix), and delayed matrices
#' (DelayedArray::DelayedMatrix).
#'
#'
#' @return A numeric vector containing the column-wise geometric mean
#' expression values. The length of the vector equals the number of columns
#' in counts.
#'
#'
#' @details
#' A pseudocount of 1 is added to all expression values before log
#' transformation to avoid undefined values for zero counts. Missing values are
#' ignored when computing the mean log expression.
#'
#'
#' Example of how to use the function
#' counts <- matrix(
#' c(0, 5, 10,
#' 2, 0, 20,
#' 1, 3, 0),
#' nrow = 3
#' )
#'
#' compute_geometric_average(counts)
#'
#'
compute_geometric_average <-function(counts){
    message("SPAROscore says: Calculating column-wise geometric averages")
    # cap_values <- apply(counts, 2,
    #                     function(x){exp(mean(log(1+ x),
    #                                          na.rm = TRUE))})

    cap_values <- exp(MatrixGenerics::colMeans(log1p(counts), na.rm = TRUE))
    return(cap_values)
}


#' Generate ranks from expression count data
#'
#' Computes column-wise gene expression ranks from a count matrix and derives
#' rank caps for each sample, cell, or spatial location. Ranking is performed
#' using MatrixGenerics::colRanks(), allowing efficient processing of dense,
#' sparse, and delayed matrix representations.
#'
#' Rank caps are incorporated by temporarily appending a row of cap expression
#' values to the count matrix before ranking. The resulting rank of each cap
#' value is returned separately from the gene rank matrix.
#'
#' @noRd
#'
#'
#' @param counts A matrix-like object containing expression counts, with
#' genes in rows and samples, cells, or spatial locations in columns.
#' Supported inputs include base matrices, sparse matrices
#' (Matrix::sparseMatrix), and delayed matrices
#' (DelayedArray::DelayedMatrix).
#'
#'
#' @param count_caps Optional numeric vector containing the expression values
#' used as rank caps for each column. The length must equal the number of
#' columns in counts. By default, column-wise geometric mean
#' expression values computed by compute_geometric_average() are used.
#'
#'
#' @param handle_ties Character string specifying how tied expression values
#' are ranked. Passed directly to MatrixGenerics::colRanks() via the
#' ties.method argument. Supported options are:
#' \describe{
#' \item{"min"}{Assign the minimum rank to tied values (default).}
#' \item{"max"}{Assign the maximum rank to tied values.}
#' \item{"average"}{Assign the average rank to tied values.}
#' \item{"random"}{Break ties at random.}
#' }
#'
#'
#' @return A named list with two elements:
#' \describe{
#' \item{ranks}{A matrix of column-wise gene ranks. Lower rank values
#' correspond to higher expression levels. The matrix is typically of type
#' integer, except when handle_ties = "average", in which case it is
#' numeric.}
#'
#' \item{rank_caps}{A numeric vector containing the rank assigned to the cap
#' value in each column.}
#' }
#'
#'
#' @details
#' Gene expression values are ranked in descending order by applying
#' MatrixGenerics::colRanks() to the negated count matrix. Consequently, the
#' most highly expressed gene in a column receives rank 1.
#'
#' To determine cap ranks, the supplied cap values are appended as an
#' additional row prior to ranking. After ranking, gene ranks and cap ranks are
#' separated and returned independently.
#'
#' Informative messages are printed indicating whether ranking is being
#' performed on a base matrix, sparse matrix, or delayed matrix object.
#'
#'
#' Example of how to use the function
#' # Compute ranks using geometric mean expression as rank caps
#' rank_results <- get_ranks_from_counts(counts)
#'
#' # Specify custom cap values
#' rank_results <- get_ranks_from_counts(
#' counts,
#' count_caps = rep(10, ncol(counts))
#' )
#'
#'
get_ranks_from_counts <- function(counts,
                                count_caps =
                                    compute_geometric_average(counts),
                                handle_ties = "min"){



    # check validity of handle_ties input
    if(!(handle_ties %in%
         c("min", "max", "average", "random"))){
        stop("SPAROscore says: Invalid value provided for handle_ties
             Limit to using 'min', 'max', 'average' or 'random'")
    }

    # Ensure rownames exist
    if(is.null(rownames(counts))){
        stop("SPAROscore says:
             counts data do not have rownames. Add gene ids/names as rownames")
    }

    # Ensure colnames exist
    if(is.null(colnames(counts))){
        stop("SPAROscore says: counts data do not have colnmaes")
    }


    # Compute rank caps' expression values or use user input
    if(length(count_caps) == ncol(counts)){
        cap_values <- count_caps
    }
    else{
        stop("SPAROscore says: Length of user entered count_caps are
        not matching the column count of counts")
    }



    # get original rownames and column names fo count_data
    count_data_rnames <- rownames(counts)
    count_data_cnames <- colnames(counts)

    # append expression caps as the last row to count_data  before ranking
    counts <- append_to_matrix_like_object(counts, cap_values)


    # print message to user as to what object is being used
    if(inherits(counts, "DelayedMatrix")){
        message("SPAROscore says: Ranking a DelayedMatrix object")
    }
    else if(inherits(counts, "sparseMatrix")){
        message("SPAROscore says: Ranking a sparseMatrix object")
    }
    else if(is.matrix(counts)){
        message("SPAROscore says: Ranking a matrix object")
    }



    # Rank each column of gene expressions along with the expression caps
    full_rank_data <- MatrixGenerics::colRanks(-counts,
                                          ties.method = handle_ties,
                                          preserveShape = TRUE,
                                          useNames =  FALSE)

    # Force matrix representation
    full_rank_data <- as.matrix(full_rank_data)

    # port column names and rownames from count_data to full_rank_data
    rownames(full_rank_data) <- c(count_data_rnames, "rank_caps")
    colnames(full_rank_data) <- count_data_cnames

    # Extract gene ranks while preserving matrix dimensions
    gene_ranks <- full_rank_data[count_data_rnames, ,drop = FALSE]

    # Extract rank caps while converting to vector
    rank_caps <- full_rank_data["rank_caps", , drop = TRUE]


    # separate the ranks and rank caps and return

    return(list(ranks = gene_ranks,
                rank_caps = rank_caps))
}






#' Compute a SPAROscore for a single sample, cell, or spatial location
#'
#' Calculates the SPAROscore for a single expression profile using the ranks of
#' genes in a signature and a corresponding rank cap. The score is based on the
#' normalized Spearman footrule distance between the observed signature gene
#' ranks and the worst-case ranking in which all genes are assigned the rank
#' cap.
#'
#' Prior to scoring, all ranks greater than rank_cap are truncated to
#' rank_cap. Missing signature genes can either be ignored or imputed with
#' the capped rank.
#'
#' @noRd
#'
#'
#' @param signature_ranks_vector Numeric or integer vector containing the ranks
#' of signature genes for a single sample, cell, or spatial location.
#'
#'
#' @param rank_cap Numeric value used to cap gene ranks before score
#' calculation. Typically corresponds to the rank of the geometric mean
#' expression value for the sample.
#'
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from the input data should be handled:
#' \itemize{
#' \item "skip" (default): exclude missing genes from the score
#' calculation.
#' \item "impute": include missing genes by assigning them the capped
#' rank value (rank_cap).
#' }
#'
#'
#' @param missing_geneset Character vector containing signature genes that are
#' not present in the input dataset.
#'
#'
#' @return A single numeric value representing the SPAROscore for the specified
#' sample, cell, or spatial location.
#'
#'
#' @details
#' The SPAROscore is computed in four steps:
#' \enumerate{
#' \item All ranks greater than rank_cap are replaced with rank_cap.
#' \item Missing genes are optionally imputed with the capped rank.
#' \item The Spearman footrule distance between the observed ranks and the
#' worst-case ranking (all genes assigned rank_cap) is calculated.
#' \item The distance is normalized by the maximum theoretically possible
#' Spearman footrule distance for a signature of the same size.
#' }
#'
#' Scores range from 0 to 1, where larger values indicate that signature genes
#' tend to have higher expression ranks relative to the cap.
#'
#'
#' Example of how to use the function
#' # validate gene signatures
#' num_genes <- nrow(counts)
#' valid_gene_signature <- validate_signature(signatures,
#' rownames(counts))
#'
#' valid_genes <- valid_gene_signature$valid_genes
#' missing_genes <- valid_gene_signature$missing_genes
#'
#' # get ranks
#' ranks <- get_ranks(counts)
#'
#'
#'
#' compute_sparoscore_per_cell(
#' signature_ranks_vector = ranks$ranks[signatures, cell_id],
#' rank_cap = ranks$rank_caps[cell_id],
#' handle_missing_genes = "skip",
#' missing_geneset = missing_genes)
#'
#'
#'
compute_sparoscore_per_cell <- function(signature_ranks_vector, rank_cap,
                               handle_missing_genes = "skip", missing_geneset){


    # check validity of input
    if(!(handle_missing_genes %in% c("skip", "impute"))){
        stop("SPAROscore says: Invalid value provided for handle_missing_genes.
             Limit to using 'skip' or 'impute'")
    }


    # Cap ranks at rank_cap, if they are more than rank_cap
    signature_ranks_vector[signature_ranks_vector > rank_cap] <- rank_cap



    # Handle missing genes in gene_signature
    if(handle_missing_genes == "impute"){
        if(length(missing_geneset) > 0){
            signature_ranks_vector <- impute_missing_gene_ranks(
                signature_ranks_vector, rank_cap, missing_geneset)
        }
    }

    # worst case scenario of ranks for the geneset is all of them being rank_cap
    # Calculate Spearman's footrule distance from the worst case ranks
    spearman_footrule <- sum(abs(rank_cap - signature_ranks_vector))



    # Calculate theoretically possible max possible Spearman footrule distance
    max_spearman_footrule <- sum(abs(
        rank_cap - seq_along(signature_ranks_vector)))


    # Calculate normalised spearman footrule distance
    final_sparoscore <- spearman_footrule/max_spearman_footrule

    return(final_sparoscore)
}



#' Compute SPAROscores for a gene signature
#'
#' Computes SPAROscores for all samples, cells, or spatial locations using
#' pre-computed ranks and rank caps. Scores quantify the enrichment of a
#' gene signature based on the normalized Spearman footrule distance between
#' observed signature gene ranks and the corresponding rank cap.
#'
#' Signature genes that are not present in the ranked dataset can either be
#' excluded from the calculation or imputed using the capped rank value.
#'
#' @noRd
#'
#'
#' @param ranks A matrix of gene ranks produced by
#' get_ranks(). Rows correspond to genes and columns
#' correspond to samples, cells, or spatial locations.
#'
#'
#' @param rank_caps Numeric vector containing the rank cap for each column of
#' ranks. Typically obtained from the rank_caps element returned by
#' get_ranks().
#'
#'
#' @param signatures Character vector containing the genes that define the
#' signature of interest. Gene identifiers must match the row names of
#' ranks.
#'
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from ranks should be handled:
#' \itemize{
#' \item "skip" (default): exclude missing genes from score
#' calculation.
#' \item "impute": assign missing genes the capped rank value and include
#' them in the score calculation.
#' }
#'
#'
#' @return A named numeric vector of SPAROscores, with one score per column of
#' the ranks matrix.
#'
#'
#' @details
#' The supplied signature is first validated against the genes present in
#' ranks. Missing genes are reported and handled according to
#' handle_missing_genes.
#'
#' For each column, gene ranks are extracted for the signature genes and passed
#' to compute_sparoscore_per_cell(), which computes a normalized Spearman
#' footrule distance score. Higher scores indicate that signature genes tend to
#' occupy higher expression ranks relative to the column-specific rank cap.
#'
#'
#' Example of how to use the function
#' # Generate ranks and rank caps
#' rank_results <- get_ranks(counts)
#'
#' # Compute scores for a gene signature
#' sparoscores <- compute_sparoscores(
#' ranks = rank_results$ranks,
#' rank_caps = rank_results$rank_caps,
#' signatures = c("CCR7", "IL7R", "LTB")
#' )
#'
#' # Include missing signature genes by imputing capped ranks
#' sparoscores <- compute_sparoscores(
#' ranks = rank_results$ranks,
#' rank_caps = rank_results$rank_caps,
#' signatures = c("CCR7", "IL7R", "LTB"),
#' handle_missing_genes = "impute"
#' )
#'
#' @seealso
#' \code{\link{get_ranks}},
#' \code{\link{compute_sparoscore_per_cell}},
#' \code{\link{validate_signature}}
#'
#'
compute_sparoscores <- function(ranks,
                                rank_caps,
                                signatures,
                                handle_missing_genes = "skip"){


    #valide the ranks matrix
    if (anyNA(ranks)) {
        stop("SPAROscore says: ranks contain NA values")
    }
    if (any(!is.finite(ranks))) {
        stop("SPAROscore says: ranks must be finite numeric values")
    }

    # Get all available genes names from the ranks
    available_genes <- rownames(ranks)

    # validate the signatures
    valid_gene_signature <- validate_signature(signatures, available_genes)

    # subset the ranks for only the valid_genes from signatures and assign
    # missing genes
    if (handle_missing_genes == "skip") {
        signature_rank_matrix <- ranks[valid_gene_signature$valid_genes, ,
                                       drop = FALSE]
        missing_genes <- character(0)
    } else {
        signature_rank_matrix <- ranks[valid_gene_signature$valid_genes, ,
                                       drop = FALSE]
        missing_genes <- valid_gene_signature$missing_genes
    }


    # get the rank cap values
    if(length(rank_caps) == ncol(ranks)){
        rank_caps <- rank_caps
    }
    else{
        stop("SPAROscore says: Invalid number of rank caps provided")
    }

    #validate ranks caps' names
    if (is.null(names(rank_caps)) || any(names(rank_caps) == "")) {
        stop("SPAROscore says: rank_caps must be a named vector")
    }
    if (!all(colnames(ranks) %in% names(rank_caps))) {
        stop("SPAROscore says: rank_caps names must match input data columns")
    }

    # Calculate sparoscores for each column
    sparoscores <- stats::setNames(rep(0.0,
                                ncol(signature_rank_matrix)),
                            colnames(signature_rank_matrix))


    sparoscores <- vapply(
        X = names(sparoscores),
        FUN = function(cell) {
            compute_sparoscore_per_cell(
                signature_ranks_vector = signature_rank_matrix[, cell],
                rank_cap = rank_caps[cell],
                handle_missing_genes = handle_missing_genes,
                missing_geneset = valid_gene_signature$missing_genes
            )
        },
        FUN.VALUE = numeric(1)
    )
    return(sparoscores)
}



#' Compute SPAROscores from a matrix-like object
#'
#' Compute SPAROscores for one or more gene signatures from a matrix-like
#' object containing either expression values or pre-computed ranks.
#'
#' This function is implemented as an S4 generic and supports multiple
#' matrix-like input classes, including dense matrices, sparse matrices,
#' delayed matrices, and data frames.
#'
#' @noRd
#'
#'
#' @param matrix_object A matrix-like object with genes in rows and samples,
#' cells, or spatial locations in columns. The object should contain either
#' expression values or pre-computed ranks, depending on the value of
#' data_has_ranks.
#'
#'
#' @param signatures Gene signature(s) to score.
#'
#' Supported inputs include:
#' \itemize{
#' \item A character vector representing a single gene signature.
#' \item A named list of character vectors representing multiple gene
#' signatures.
#' \item A GeneSet object.
#' \item A GeneSetCollection object.
#' }
#'
#' Gene identifiers must match the row names of matrix_object.
#'
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
#' @param data_has_ranks Logical indicating whether matrix_object contains
#' pre-computed ranks. If TRUE, SPAROscores are computed directly from the
#' supplied ranks. If FALSE (default), ranks are first computed from the
#' supplied expression values.
#'
#'
#' @param count_caps Optional numeric vector of expression thresholds used to
#' determine rank caps during rank calculation. Only used when
#' data_has_ranks = FALSE. If NULL (default), values are computed using
#' compute_geometric_average().
#'
#'
#' @param rank_caps Optional numeric vector of rank cap values used during
#' SPAROscore calculation. If NULL (default), rank caps returned by get_ranks()
#' are used. When data_has_ranks = TRUE, this argument must be supplied.
#'
#'
#' @param handle_ties Character string specifying how tied expression values
#' are ranked. Passed to get_ranks(). Supported values are
#' "min", "max", "average", and "random".
#'
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from the dataset should be handled:
#' \itemize{
#' \item "skip" (default): exclude missing genes from score calculation.
#' \item "impute": include missing genes by assigning them the capped
#' rank value.
#' }
#'
#' @param prefix Character string to be appended before the headers of the
#' returned scores. Defaults to "".
#'
#'
#' @returns A numeric matrix of SPAROscores with one row per column of
#' matrix_object.
#'
#' \itemize{
#' \item For a single signature, the returned matrix contains one column
#' named "SPAROscore".
#' \item For multiple signatures, columns correspond to signatures and rows
#' correspond to samples, cells, or spatial locations.
#' }
#'
#'
#' @details
#' When data_has_ranks = FALSE, feature ranks are computed using
#' get_ranks() before SPAROscores are calculated. When
#' data_has_ranks = TRUE, the supplied ranks are used directly.
#'
#' SPAROscores are computed using get_scores().
#'
#' Example of how to use the function
#' # Compute scores directly from expression data
#' scores <- augment_sparoscores_matrix(
#' matrix_object = counts,
#' signatures = c("CCR7", "IL7R", "LTB")
#' )
#'
#' # Compute scores from pre-computed ranks
#' scores <- augment_sparoscores_matrix(
#' matrix_object = ranks,
#' signatures = c("CCR7", "IL7R", "LTB"),
#' data_has_ranks = TRUE,
#' rank_caps = rank_caps
#' )
#'
#' @seealso
#' \code{\link{get_ranks}},
#' \code{\link{get_scores}},
#' \code{\link{compute_geometric_average}}
#'
#'
augment_sparoscores_matrix <- function(matrix_object,
                                       signatures,
                                       down_signatures = NULL,
                                       data_has_ranks = FALSE,
                                       count_caps = NULL,
                                       rank_caps = NULL,
                                       handle_ties = "min",
                                       handle_missing_genes = "skip",
                                       prefix = ""){
    if(!is.logical(data_has_ranks)){
        stop("SPAROscore says: data_has_ranks should be a boolean")
    }

    # if ranks are provided, directly calculate scores
    if(data_has_ranks){
        if(is.null(rank_caps)){
            stop("SPAROscore says: rank_caps needed when data_has_ranks = TRUE")
        }

        sparoscores <- get_scores(ranks = matrix_object,
                                  rank_caps = rank_caps,
                                  signatures = signatures,
                                  down_signatures = down_signatures,
                                  handle_missing_genes = handle_missing_genes,
                                  prefix = prefix)
    }
    # if counts are provides, calculate ranks, rank_caps and then score
    else{
        # if user provides count_caps, then use it, if not, calculate default
        if(is.null(count_caps)){
            count_caps <- compute_geometric_average(matrix_object)
        }
        else{
            count_caps <- count_caps
        }

        # compute ranks and rank_caps
        ranks_output <- get_ranks(counts = matrix_object,
                                      count_caps = count_caps,
                                      handle_ties = handle_ties)

        #if user gives rank_caps, use them. If not, use default
        if(is.null(rank_caps)){
            rank_caps <- ranks_output$rank_caps
        }
        else{
            count_caps <- count_caps
        }

        sparoscores <- get_scores(ranks = ranks_output$ranks,
                                  rank_caps = rank_caps,
                                  signatures = signatures,
                                  down_signatures = down_signatures,
                                  handle_missing_genes = handle_missing_genes,
                                  prefix = prefix)
    }

    return(sparoscores)
}




#' Add SPAROscores to a Seurat object
#'
#' Computes SPAROscores for one or more gene signature and stores the resulting
#' scores in the metadata of a Seurat object.
#' Ranks can either be calculated from
#' an expression layer or retrieved from previously stored ranks.
#'
#' @noRd
#'
#'
#' @param seurat_object A Seurat object containing expression/rank data.
#'
#'
#' @param signatures Gene signature(s) to score.
#'
#' Supported inputs include:
#' \itemize{
#' \item A character vector representing a single gene signature.
#' \item A named list of character vectors representing multiple gene
#' signatures.
#' \item A GeneSet object.
#' \item A GeneSetCollection object.
#' }
#'
#'#' Gene identifiers must use the same naming convention as the row names of
#' ranks.
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
#' @param data_has_ranks Logical indicating whether previously
#' calculated ranks are stored in the assay/layer.
#' If TRUE, ranks are read directly from the assay layer.
#' If FALSE (default), ranks are computed from counts data in the assay layer.
#'
#'
#' @param assay Character string specifying the name of
#' the assay containing the data to use.
#' When data_has_ranks = FALSE, this assay should contain expression counts.
#' When data_has_ranks = TRUE, it should contain pre-computed ranks.
#' Defaults to "RNA"
#'
#'
#' @param layer Character string specifying the name
#' of the assay layer containing the data to use.
#' When data_has_ranks = FALSE, this layer should contain expression counts.
#' When data_has_ranks = TRUE, it should contain pre-computed ranks.
#' Defaults to "counts"
#'
#'
#' @param count_caps Optional numeric vector of expression values used to
#' derive rank caps during rank calculation. If NULL, geometric mean
#' expression values are computed using compute_geometric_average().
#'
#'
#' @param rank_caps Optional numeric vector of rank cap values used during
#' SPAROscore calculation. If NULL, rank caps returned by
#' get_ranks() are used.
#'
#'
#' @param handle_ties Character string specifying how tied expression values
#' are ranked. Passed to get_ranks(). Supported values are
#' "min", "max", "average", and "random".
#'
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from the dataset should be handled:
#' \itemize{
#' \item "skip" (default): exclude missing genes from score calculation.
#' \item "impute": include missing genes by assigning them the capped
#' rank value.
#' }
#'
#' @param prefix Character string to be appended before the headers of the
#' columns in the metada with the returned scores and rank_caps. Defaults to "".
#'
#'
#' @return A Seurat object with:
#' \itemize{
#' \item One or more metadata columns containing SPAROscores for the supplied
#' signature(s).
#' \item A metadata column named "rank_caps" containing the rank cap for
#' each cell.
#' \item A "ranks" layer added to the selected assay when ranks are
#' calculated from expression data.
#' }
#'
#' @details
#' When data_has_ranks = FALSE, the function:
#' \enumerate{
#' \item Retrieves expression data from the specified assay and layer.
#' \item Computes feature ranks using get_ranks().
#' \item Calculates SPAROscores using get_scores().
#' \item Stores the computed ranks in a new assay layer named "ranks".
#' }
#'
#' When data_has_ranks = TRUE, the function:
#' \enumerate{
#' \item Retrieves pre-computed ranks from the specified assay and layer.
#' \item Uses the supplied rank_caps or retrieves them from the
#' "rank_caps" metadata column.
#' \item Calculates SPAROscores using get_scores().
#' }
#'
#' Computed scores and rank caps are appended to the Seurat object's metadata
#' and can be accessed with seurat_object[[]].
#'
#' Example of how to use the function
#' # Compute SPAROscores from RNA counts
#' seurat_object <- augment_sparoscores_seurat(
#' seurat_object = seurat_object,
#' signatures = c("CCR7", "IL7R", "LTB")
#' )
#'
#' # Compute scores using pre-computed ranks stored in a ranks layer
#' seurat_object <- augment_sparoscores_seurat(
#' seurat_object = seurat_object,
#' signatures = c("CCR7", "IL7R", "LTB"),
#' data_has_ranks = TRUE,
#' layer = "ranks"
#' )
#'
#' @seealso
#' \code{\link{get_ranks}},
#' \code{\link{get_scores}},
#' \code{\link{compute_geometric_average}}
#'
#'
augment_sparoscores_seurat <- function(seurat_object,
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


    if(!is.logical(data_has_ranks)){
        stop("SPAROscore says: data_has_ranks should be a boolean")
    }

    if(!(assay %in% Seurat::Assays(seurat_object))){
        stop("SPAROscore says: assay passed is not in the seurat_object")
    }

    if(!(layer %in% seurat_object[[assay]][])){
        stop("SPAROscore says: layer passed is not in the seurat_object")
    }

    # if user ask to use pre-calculated ranks
    if(data_has_ranks){
        # extract precalculated ranks and rank_caps
        ranks <- Seurat::GetAssayData(seurat_object,
                                   assay = assay,
                                   layer = layer)

        # set rank_caps to pre-calculated values unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- as.numeric(unlist(seurat_object[["rank_caps"]]))
            names(rank_caps) <- rownames(seurat_object[[]])
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_scores(ranks =  ranks,
                                rank_caps = rank_caps,
                                signatures =  signatures,
                                down_signatures = down_signatures,
                                handle_missing_genes = handle_missing_genes,
                                prefix = prefix)
    }
    else{
        #extract count matrix from data
        counts <- Seurat::GetAssayData(seurat_object,
                                        assay = assay,
                                        layer = layer)

        # set count caps to geometric averages unless provided by the user
        if(is.null(count_caps)){
            count_caps <- compute_geometric_average(counts)
        }
        else{
            count_caps <- count_caps
        }


        # get outputs for ranking function
        get_ranks_output <- get_ranks(counts = counts,
                                      count_caps = count_caps,
                                      handle_ties = handle_ties)

        # set ranks for future input
        ranks <- get_ranks_output$ranks


        # set rank_caps to geometric averages unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- get_ranks_output$rank_caps
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_scores(ranks =  ranks,
                                rank_caps = rank_caps,
                                signatures =  signatures,
                                down_signatures = down_signatures,
                                handle_missing_genes = handle_missing_genes,
                                prefix = prefix)

        # store the ranks in new layer named ranks
        seurat_object <- Seurat::SetAssayData(object = seurat_object,
                                              assay = assay,
                                              layer = "ranks",
                                              new.data = ranks)
    }



    # append the rank_caps and sparoscores to metadata
    sparoscores_df <- data.frame(rank_caps)
    colnames(sparoscores_df) <- "rank_caps"

    sparoscores_df <- cbind(sparoscores_df, sparoscores)

    seurat_object <- Seurat::AddMetaData(object = seurat_object,
                                         metadata = sparoscores_df)

    return(seurat_object)
}


#' Add SPAROscores to a SummarizedExperiment-derived object
#'
#' Compute SPAROscores for a gene signature and store the resulting scores in
#' the colData of a SummarizedExperiment, SingleCellExperiment,
#' SpatialExperiment, or RangedSummarizedExperiment object.
#'
#' Depending on data_has_ranks, the function either computes feature ranks
#' from an expression assay or uses pre-computed ranks stored in an assay.
#' Rank caps are calculated automatically when not supplied.
#'
#' @noRd
#'
#'
#' @param sce_object A SummarizedExperiment-derived object containing
#' expression data or pre-computed ranks.
#'
#' @param signatures Gene signature(s) to score.
#'
#' Supported inputs include:
#' \itemize{
#' \item A character vector representing a single gene signature.
#' \item A named list of character vectors representing multiple gene
#' signatures.
#' \item A GeneSet object.
#' \item A GeneSetCollection object.
#' }
#'
#' Gene identifiers must use the same naming convention as the row names of
#' ranks.
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
#' @param data_has_ranks Logical indicating whether the specified assay already
#' contains feature ranks. If TRUE, ranks are read directly from the assay.
#' If FALSE (default), ranks are computed from expression data.
#'
#' @param assay Name of the assay containing the data to use. When
#' data_has_ranks = FALSE, this assay should contain expression counts.
#' When data_has_ranks = TRUE, it should contain pre-computed ranks.
#' Defaults to "counts".
#'
#' @param count_caps Optional numeric vector of expression thresholds used to
#' determine rank caps during rank calculation. If NULL, values are computed
#' using compute_geometric_average().
#'
#' @param rank_caps Optional numeric vector of rank cap values used during
#' SPAROscore calculation. If NULL, rank caps are obtained from
#' get_ranks() or retrieved from the "rank_caps" column of colData
#' when data_has_ranks = TRUE.
#'
#' @param handle_ties Method used to rank tied expression values. Passed to
#' get_ranks(). Supported values are "min", "max", "average",
#' and "random".
#'
#' @param handle_missing_genes Strategy for handling signature genes that are
#' absent from the dataset:
#' \itemize{
#' \item "skip" (default): exclude missing genes from score calculation.
#' \item "impute": assign missing genes the capped rank value and include
#' them in score calculation.
#' }
#'
#' @param prefix Character string to be appended before the headers of the
#' columns in colData with returned scores. Defaults to "".
#'
#'
#' @return The input object with:
#' \itemize{
#' \item One or more SPAROscore columns added to colData().
#' \item A "rank_caps" column added to colData().
#' \item A "ranks" assay added when ranks are computed from expression data.
#' }
#'
#' @details
#' When data_has_ranks = FALSE, the function:
#' \enumerate{
#' \item Retrieves expression values from the specified assay.
#' \item Computes feature ranks using get_ranks().
#' \item Calculates SPAROscores using get_scores().
#' \item Stores the computed ranks in a new assay named "ranks".
#' }
#'
#' When data_has_ranks = TRUE, the function:
#' \enumerate{
#' \item Retrieves pre-computed ranks from the specified assay.
#' \item Uses the supplied rank_caps or retrieves them from the
#' "rank_caps" column of colData().
#' \item Calculates SPAROscores using get_scores().
#' }
#'
#' Computed scores and rank caps are appended to colData() and can be
#' accessed using SummarizedExperiment::colData().
#'
#' Example of how to use the function
#' library(SingleCellExperiment)
#'
#' sce <- augment_sparoscores_sce(
#' sce_object = sce,
#' signatures = c("CCR7", "IL7R", "LTB")
#' )
#'
#' # Reuse previously computed ranks
#' sce <- augment_sparoscores_sce(
#' sce_object = sce,
#' signatures = c("CCR7", "IL7R", "LTB"),
#' data_has_ranks = TRUE,
#' assay = "ranks"
#' )
#'
#' @seealso
#' \code{\link{get_ranks}},
#' \code{\link{get_scores}},
#' \code{\link{compute_geometric_average}}
#'
#'
augment_sparoscores_sce <- function(sce_object,
                                    signatures,
                                    down_signatures = NULL,
                                    data_has_ranks = FALSE,
                                    assay = "counts",
                                    count_caps = NULL,
                                    rank_caps = NULL,
                                    handle_ties = "min",
                                    handle_missing_genes = "skip",
                                    prefix = ""){


    if(!is.logical(data_has_ranks)){
        stop("SPAROscore says: data_has_ranks should be a boolean")
    }

    if(!(assay %in% SummarizedExperiment::assayNames(sce_object))){
        stop("SPAROscore says: assay passed is not in the data")
    }

    # if user ask to use pre-calculated ranks
    if(data_has_ranks){
        # extract precalculated ranks and rank_caps
        ranks <- SummarizedExperiment::assay(sce_object, assay)

        # set rank_caps to pre-calculated values unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- as.numeric(
                                    unlist(
                                        SummarizedExperiment::colData(
                                            sce_object)[["rank_caps"]]))

            names(rank_caps) <- rownames(
                SummarizedExperiment::colData(sce_object))
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_scores(ranks =  ranks,
                                  rank_caps = rank_caps,
                                  signatures =  signatures,
                                  down_signatures = down_signatures,
                                  handle_missing_genes = handle_missing_genes,
                                  prefix = prefix)
    }
    else{
        #extract count matrix from data
        counts <- SummarizedExperiment::assay(sce_object, assay)

        # set count caps to geometric averages unless provided by the user
        if(is.null(count_caps)){
            count_caps <- compute_geometric_average(counts)
        }
        else{
            count_caps <- count_caps
        }


        # get outputs for ranking function
        get_ranks_output <- get_ranks(counts = counts,
                                      count_caps = count_caps,
                                      handle_ties = handle_ties)

        # set ranks for future input
        ranks <- get_ranks_output$ranks


        # set rank_caps to geometric averages unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- get_ranks_output$rank_caps
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_scores(ranks =  ranks,
                                  rank_caps = rank_caps,
                                  signatures =  signatures,
                                  down_signatures = down_signatures,
                                  handle_missing_genes = handle_missing_genes,
                                  prefix = prefix)

        # store the ranks in new assay named ranks
        SummarizedExperiment::assay(sce_object, "ranks") <- ranks
    }



    # append the rank_caps and sparoscores to metadata
    sparoscores_df <- data.frame(rank_caps)
    colnames(sparoscores_df) <- "rank_caps"

    sparoscores_df <- cbind(sparoscores_df, sparoscores)

    SummarizedExperiment::colData(sce_object) <- cbind(
        SummarizedExperiment::colData(sce_object), sparoscores_df
    )

    return(sce_object)
}


