
# ---------------------------------------------------------------------------
# Exponential tilt f(x) proportional to exp(theta * x) on the unit interval.
#
# Positions are expressed as fractions of the ranking, 0 at the bottom and 1 at
# the top, so a score is directly comparable across cells with different
# sparsity. Z_{0:1}(theta) = (e^theta - 1)/theta is the normalising area and
# Z_{0:q}(theta) = (e^{theta q} - 1)/theta the same area stopped at q.
#
# One unconstrained parameter: theta > 0 is an UP signature, theta < 0 a DOWN
# one, and the family passes smoothly through the uniform at theta = 0. A Beta
# parameterisation kinks there and destroys the curvature estimate.
# ---------------------------------------------------------------------------

# log Z_{0:1}(theta). Both expm1(theta) and theta flip sign together, so the
# ratio is positive on either side and can be logged on magnitudes.
.tilt_logZ <- function(theta) {
    if (abs(theta) < 1e-8) {
        log1p(theta / 2 + theta^2 / 6)
    } else {
        log(abs(expm1(theta))) - log(abs(theta))
    }
}

# Mean position implied by theta, which is the reported score. This is
# d/dtheta log Z_{0:1}(theta), pre-solved so nothing is integrated at runtime.
.tilt_mean <- function(theta) {
    if (abs(theta) < 1e-8) {
        0.5 + theta / 12
    } else {
        1 / (1 - exp(-theta)) - 1 / theta
    }
}

# log Z_{0:q}(theta) - log Z_{0:1}(theta), the contribution of a single
# censored gene. The 1/theta inside each area cancels between them.
#
# The negative-theta branch uses log1p on quantities that decay to zero rather
# than expm1 on quantities that saturate at -1. Written as
# expm1(theta * q)/expm1(theta) the ratio rounds to exactly 1 once
# abs(theta) * q exceeds log(4 / .Machine$double.eps), about 37.4, which flattens
# the likelihood and leaves a fully censored fit resting on an arbitrary point.
# This form stays strictly monotone until exp() itself underflows near -745.
.tilt_log_ratio <- function(q, theta) {
    if (q <= 0) return(-Inf)
    if (q >= 1) return(0)
    if (abs(theta) < 1e-8) return(log(q))
    if (theta > 0) {
        log(expm1(theta * q)) - log(expm1(theta))
    } else {
        log1p(-exp(theta * q)) - log1p(-exp(theta))
    }
}

# E[x | x <= q], the position the fit implies for every censored gene. Equal to
# q * .tilt_mean(theta * q), because truncating the tilt to [0, q] gives back a
# tilt on [0, 1] with parameter theta * q.
.tilt_cond_mean <- function(q, theta) {
    if (abs(theta) < 1e-8) {
        q / 2
    } else {
        q / (1 - exp(-theta * q)) - 1 / theta
    }
}


#' Fit the censored tilt to one cell
#'
#' Solves the censored maximum likelihood problem for a single column and
#' returns the fitted mean position together with its standard error.
#'
#' @noRd
#'
#' @param positions Numeric vector of detected gene positions on the unit
#' interval, all strictly above `cap_quantile`.
#'
#' @param censored_size Number of signature genes known only to fall below the
#' cap.
#'
#' @param cap_quantile Fraction of the ranking that sits at or below the cap.
#'
#' @param theta_bound Half-width of the search interval for `theta`.
#'
#' @return A numeric vector of length three: the mean position, its standard
#' error, and the fitted `theta`.
#'
#' @details
#' The log-likelihood gives each detected gene the log density at its known
#' position and each censored gene the log probability of the whole region below
#' the cap:
#'
#' \deqn{\ell(\theta) = \theta \sum_{det} x_i - k \log Z_{0:1}(\theta)
#'       + m \left( \log Z_{0:c}(\theta) - \log Z_{0:1}(\theta) \right)}
#'
#' Its derivative is the moment condition
#'
#' \deqn{\sum_{det} x_i + m \cdot E[x \mid x \le c] - S \mu(\theta)}
#'
#' which is strictly decreasing in `theta`, so the fit is found by root finding
#' rather than by a general optimiser. The censored genes are neither dropped
#' nor assigned a made-up rank; they enter only through their count and the cap.
.fit_tilt_per_cell <- function(positions, censored_size, cap_quantile,
                               theta_bound = 200) {

    detected_size <- length(positions)
    total_size <- detected_size + censored_size

    # With nothing detected the data place no upper constraint on the fit, so
    # the likelihood is monotone and any reported position would be an artefact
    # of where the search stopped. Say so instead.
    if (detected_size == 0L) {
        return(c(score = NA_real_, se = Inf, theta = NA_real_))
    }

    position_sum <- sum(positions)

    negative_loglik <- function(theta) {
        censored_term <- if (censored_size > 0L) {
            censored_size * .tilt_log_ratio(cap_quantile, theta)
        } else {
            0
        }
        -(theta * position_sum - detected_size * .tilt_logZ(theta) +
              censored_term)
    }

    # Derivative of the log-likelihood: observed positions plus the imputed
    # censored mass, against the mean the candidate theta implies.
    score_equation <- function(theta) {
        censored_mass <- if (censored_size > 0L) {
            censored_size * .tilt_cond_mean(cap_quantile, theta)
        } else {
            0
        }
        position_sum + censored_mass - total_size * .tilt_mean(theta)
    }

    lower <- score_equation(-theta_bound)
    upper <- score_equation(theta_bound)

    # A signature packed against either end of the ranking puts the root outside
    # any finite interval, so the bound itself is the answer.
    theta <- if (lower <= 0) {
        -theta_bound
    } else if (upper >= 0) {
        theta_bound
    } else {
        stats::uniroot(score_equation, c(-theta_bound, theta_bound),
                       f.lower = lower, f.upper = upper, tol = 1e-9)$root
    }

    # Observed information by central difference, carried onto the 0-1 mean
    # scale by the delta method. As censoring takes over, the peak flattens and
    # this grows without bound, which is the estimate declaring itself unusable.
    step <- 1e-3
    curvature <- (negative_loglik(theta + step) - 2 * negative_loglik(theta) +
                      negative_loglik(theta - step)) / step^2
    slope <- (.tilt_mean(theta + step) - .tilt_mean(theta - step)) / (2 * step)

    c(score = .tilt_mean(theta),
      se = abs(slope) / sqrt(max(curvature, 0)),
      theta = theta)
}


#' Score one signature across every column
#'
#' @noRd
#'
#' @param ranks Numeric matrix of gene ranks, genes in rows and columns as
#' returned by [get_ranks()]. Rank 1 is the most highly expressed gene.
#'
#' @param rank_caps Named numeric vector with one rank cap per column.
#'
#' @param signature Character vector of signature genes.
#'
#' @param handle_missing_genes Either "skip" or "impute".
#'
#' @param n_items Size of the ranking universe used to convert ranks into
#' positions on the unit interval.
#'
#' @param theta_bound Half-width of the search interval for `theta`.
#'
#' @return A two-column matrix with one row per column of `ranks`, holding the
#' score and its standard error.
.compute_sparoscores_mle <- function(ranks,
                                     rank_caps,
                                     signature,
                                     handle_missing_genes = "skip",
                                     n_items = NULL,
                                     theta_bound = 200) {

    if (!(handle_missing_genes %in% c("skip", "impute"))) {
        stop("SPAROscore says: Invalid value provided for handle_missing_genes.
             Limit to using 'skip' or 'impute'")
    }
    if (anyNA(ranks)) {
        stop("SPAROscore says: ranks contain NA values")
    }
    if (!all(is.finite(as.matrix(ranks)))) {
        stop("SPAROscore says: ranks must be finite numeric values")
    }
    if (length(rank_caps) != ncol(ranks)) {
        stop("SPAROscore says: Invalid number of rank caps provided")
    }
    if (is.null(names(rank_caps)) || any(names(rank_caps) == "")) {
        stop("SPAROscore says: rank_caps must be a named vector")
    }
    if (!all(colnames(ranks) %in% names(rank_caps))) {
        stop("SPAROscore says: rank_caps names must match input data columns")
    }

    valid_gene_signature <- validate_signature(signature, rownames(ranks))
    signature_rank_matrix <- as.matrix(
        ranks[valid_gene_signature$valid_genes, , drop = FALSE])

    # get_ranks() ranks the genes alongside an appended row of cap values, so
    # the universe the ranks were drawn from is one larger than the gene count.
    if (is.null(n_items)) {
        n_items <- nrow(ranks) + 1L
    }

    # Genes absent from the data are unobserved in exactly the sense the model
    # already handles, so imputing them means counting them as censored rather
    # than inventing a rank.
    extra_censored <- if (handle_missing_genes == "impute") {
        length(valid_gene_signature$missing_genes)
    } else {
        0L
    }

    fitted <- vapply(
        X = colnames(signature_rank_matrix),
        FUN = function(cell) {
            gene_ranks <- signature_rank_matrix[, cell]
            cap <- rank_caps[[cell]]

            detected <- gene_ranks < cap
            positions <- (n_items + 0.5 - gene_ranks[detected]) / n_items
            cap_quantile <- min(max((n_items + 1 - cap) / n_items, 0), 1)

            .fit_tilt_per_cell(
                positions = positions,
                censored_size = sum(!detected) + extra_censored,
                cap_quantile = cap_quantile,
                theta_bound = theta_bound)[c("score", "se")]
        },
        FUN.VALUE = numeric(2))

    t(fitted)
}


# Normalise every supported signature container down to a named list of
# character vectors, so the scoring path itself never branches on input class.
.as_signature_list <- function(signatures, default_name) {
    if (inherits(signatures, "GeneSetCollection")) {
        return(GSEABase::geneIds(signatures))
    }
    if (inherits(signatures, "GeneSet")) {
        return(stats::setNames(list(GSEABase::geneIds(signatures)),
                               GSEABase::setName(signatures)))
    }
    if (is.list(signatures)) {
        if (is.null(names(signatures)) || any(names(signatures) == "")) {
            stop("SPAROscore says: signatures supplied as a list must be named")
        }
        return(signatures)
    }
    stats::setNames(list(signatures), default_name)
}


# Score every signature in a named list and bind the results into matching
# score and standard-error matrices.
.score_signature_list <- function(ranks, rank_caps, signature_list,
                                  handle_missing_genes, n_items, theta_bound) {

    per_signature <- lapply(signature_list, function(genes) {
        .compute_sparoscores_mle(
            ranks = ranks,
            rank_caps = rank_caps,
            signature = genes,
            handle_missing_genes = handle_missing_genes,
            n_items = n_items,
            theta_bound = theta_bound)
    })

    scores <- vapply(per_signature, function(x) x[, "score"],
                     FUN.VALUE = numeric(ncol(ranks)))
    errors <- vapply(per_signature, function(x) x[, "se"],
                     FUN.VALUE = numeric(ncol(ranks)))

    scores <- matrix(scores, nrow = ncol(ranks),
                     dimnames = list(colnames(ranks), names(signature_list)))
    errors <- matrix(errors, nrow = ncol(ranks),
                     dimnames = list(colnames(ranks), names(signature_list)))

    list(score = scores, se = errors)
}


# Shared body for every matrix-like method: turn counts into ranks when needed,
# then score. Mirrors augment_sparoscores_matrix() for the shipped score.
.augment_sparoscores_mle_matrix <- function(matrix_object,
                                            signatures,
                                            down_signatures,
                                            data_has_ranks,
                                            count_caps,
                                            rank_caps,
                                            handle_ties,
                                            handle_missing_genes,
                                            prefix,
                                            theta_bound,
                                            return_se) {

    if (!is.logical(data_has_ranks)) {
        stop("SPAROscore says: data_has_ranks should be a boolean")
    }

    if (data_has_ranks) {
        if (is.null(rank_caps)) {
            stop("SPAROscore says: rank_caps needed when data_has_ranks = TRUE")
        }
        ranks <- matrix_object
    } else {
        if (is.null(count_caps)) {
            count_caps <- compute_geometric_average(matrix_object)
        }
        ranks_output <- get_ranks(counts = matrix_object,
                                  count_caps = count_caps,
                                  handle_ties = handle_ties)
        ranks <- ranks_output$ranks
        if (is.null(rank_caps)) {
            rank_caps <- ranks_output$rank_caps
        }
    }

    signature_list <- .as_signature_list(signatures, "SPAROscore_MLE")
    fitted <- .score_signature_list(
        ranks = ranks, rank_caps = rank_caps, signature_list = signature_list,
        handle_missing_genes = handle_missing_genes, n_items = NULL,
        theta_bound = theta_bound)

    # Both scores are means on the same 0-1 scale, so the contrast is their
    # difference and the errors add in quadrature.
    if (!is.null(down_signatures)) {
        down_list <- .as_signature_list(down_signatures, "SPAROscore_MLE")
        if (!identical(names(signature_list), names(down_list))) {
            stop("SPAROscore says: The names of gene sets in
                      down_signatures should match that of signatures")
        }
        down_fitted <- .score_signature_list(
            ranks = ranks, rank_caps = rank_caps, signature_list = down_list,
            handle_missing_genes = handle_missing_genes, n_items = NULL,
            theta_bound = theta_bound)
        fitted$score <- fitted$score - down_fitted$score
        fitted$se <- sqrt(fitted$se^2 + down_fitted$se^2)
    }

    colnames(fitted$score) <- paste0(prefix, colnames(fitted$score))
    colnames(fitted$se) <- paste0(prefix, colnames(fitted$se))

    if (return_se) {
        return(fitted)
    }

    scores <- fitted$score
    attr(scores, "standard_error") <- fitted$se
    scores
}


#' Compute SPAROscores by censored maximum likelihood
#'
#' Compute signature scores for one or more gene signatures by fitting a
#' censored maximum likelihood model to the observed ranks, rather than by
#' normalising a Spearman footrule distance as [sparoscore()] does.
#'
#' `sparoscore_mle()` is an S4 generic with methods for matrix-like objects. It
#' takes the same arguments as [sparoscore()] and returns a matrix of the same
#' shape, so calls can be swapped over without changing surrounding code. In
#' addition it reports a per-cell standard error, which makes the cases where
#' sparsity has destroyed the information visible instead of silent.
#'
#' The two scores are not numerically interchangeable. Both run from 0 to 1,
#' but this one is a mean position in the ranking, so a signature spread like
#' the background scores about 0.5, whereas the footrule score of [sparoscore()]
#' approaches 0 in that case. Compare across cells within one scoring method,
#' not across methods.
#'
#' @param data Input data object. Supported classes include:
#' \itemize{
#'   \item Dense matrices (`matrix`)
#'   \item Sparse matrices (`sparseMatrix`)
#'   \item Delayed matrices (`DelayedMatrix`)
#'   \item Data frames
#' }
#'
#' @param signatures Gene signature(s) to score. Supported inputs include:
#' \itemize{
#'   \item A character vector representing a single gene signature.
#'   \item A named list of character vectors representing multiple signatures.
#'   \item A `GeneSet` object.
#'   \item A `GeneSetCollection` object.
#' }
#'
#' @param down_signatures Gene signature(s) to be considered for scoring the
#' down-regulation effect. Supported inputs are the same as `signatures`. If
#' provided, `names(down_signatures)` must match `names(signatures)`. Defaults
#' to `NULL`.
#'
#' When `down_signatures` is not `NULL`,
#' `Final Score <- Score(signatures) - Score(down_signatures)`, and the
#' standard errors of the two fits are combined in quadrature.
#'
#' @param data_has_ranks Logical indicating whether the supplied data already
#' contains feature ranks. If `FALSE` (default), ranks are computed from
#' expression values prior to scoring.
#'
#' @param assay Ignored for matrix-like inputs. Present so that the argument
#' list matches [sparoscore()].
#'
#' @param layer Ignored for matrix-like inputs. Present so that the argument
#' list matches [sparoscore()].
#'
#' @param count_caps Optional numeric vector containing expression threshold
#' values for each column in the counts data. These values are used to
#' determine the corresponding rank caps. Only used when
#' `data_has_ranks = FALSE`. If `NULL` (default), values are computed using
#' [compute_geometric_average()].
#'
#' @param rank_caps An optional named numeric vector containing the rank cap
#' associated with each column of ranks data. Required when
#' `data_has_ranks = TRUE`; otherwise obtained from [get_ranks()].
#'
#' @param handle_ties Character string specifying how tied expression values
#' should be ranked. Passed to `MatrixGenerics::colRanks(ties.method = ...)`.
#' Supported values are `"min"` (default), `"max"`, `"average"` and
#' `"random"`.
#'
#' Note that `"average"` is a better choice here than the package default.
#' Sparse data produce large blocks of tied zero counts, and `"min"` gives
#' every gene in a tie the most favourable rank in the block, which biases
#' scores upward as sparsity increases.
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from the ranks should be handled. Supported options are:
#' \describe{
#' \item{"skip"}{Exclude missing genes from the score calculation (default).}
#' \item{"impute"}{Count missing genes as censored, that is, as known to fall
#' below the cap without assigning them a rank.}
#' }
#'
#' @param prefix Character string prepended to the names of the returned score
#' columns. Defaults to `""`.
#'
#' @param theta_bound Half-width of the interval searched for the tilt
#' parameter. The reportable score range follows from it, roughly
#' `1/theta_bound` to `1 - 1/theta_bound`, so the default of 200 admits scores
#' between about 0.005 and 0.995. Raise it only for signatures packed against
#' the very ends of the ranking.
#'
#' @param return_se Logical. If `FALSE` (default), return a matrix of scores
#' with the standard errors attached as the `"standard_error"` attribute, so
#' the result is shaped exactly like the output of [sparoscore()]. If `TRUE`,
#' return a list with elements `score` and `se`.
#'
#' @return
#' A numeric matrix of scores with one row per sample, cell, or spatial domain
#' and one column per signature, carrying the standard errors in its
#' `"standard_error"` attribute. When `return_se = TRUE`, a list of two such
#' matrices named `score` and `se` instead.
#'
#' Scores are mean positions in the ranking. A value of 0.5 means the signature
#' sits where an arbitrary set of genes would, above 0.5 means it is
#' concentrated towards the highly expressed end, and below 0.5 towards the
#' lowly expressed end.
#'
#' A score of `NA` with an infinite standard error marks a column in which no
#' signature gene was detected above the cap. The data genuinely cannot place
#' the signature in that case, so no number is reported.
#'
#' @details
#' Ranks below the cap are not observed: sparsity has collapsed them into ties
#' at or near zero, and their order carries no information. The shipped
#' [sparoscore()] handles this by truncating those ranks at the cap and
#' normalising the resulting footrule, which leaves the score depending on where
#' the cap happens to fall.
#'
#' This function instead treats the cap as a censoring point and estimates the
#' position of the signature in the full ranking. Each detected gene
#' contributes the log density at its known position; each censored gene
#' contributes only the log probability of falling below the cap, so its unknown
#' rank never enters. Maximising that likelihood over a one-parameter tilt
#' family, `f(x)` proportional to `exp(theta * x)`, is equivalent to matching
#' the mean position of all signature genes, with the censored ones filled in at
#' their conditional mean. The reported score is the fitted mean position.
#'
#' Because the cap enters the likelihood rather than the normaliser, the score
#' is close to invariant to the cap, whereas footrule-based normalisation drifts
#' with it. The cost is that a heavily censored signature is estimated from
#' little information, and the standard error is what reports this: it grows by
#' orders of magnitude as the detected count falls. Downstream analyses should
#' weight by `1 / se^2` rather than test any boundary condition.
#'
#' @seealso [sparoscore()], [get_ranks()], [get_scores()]
#'
#' @examples
#' counts <- matrix(
#'   sample(0:10, 500, replace = TRUE),
#'   nrow = 50,
#'   dimnames = list(
#'     paste0("gene", 1:50),
#'     paste0("cell", 1:10)
#'   )
#' )
#'
#' scores <- sparoscore_mle(
#'   data = counts,
#'   signatures = c("gene1", "gene2", "gene3")
#' )
#'
#' # the uncertainty travels with the score
#' attr(scores, "standard_error")
#'
#' # or ask for both explicitly
#' fitted <- sparoscore_mle(
#'   data = counts,
#'   signatures = c("gene1", "gene2", "gene3"),
#'   return_se = TRUE
#' )
#'
#' # multiple signatures, and average ranks for tied counts
#' scores <- sparoscore_mle(
#'   data = counts,
#'   signatures = list(
#'     SignatureA = c("gene1", "gene2", "gene3"),
#'     SignatureB = c("gene10", "gene11", "gene12")
#'   ),
#'   handle_ties = "average"
#' )
#'
#' @export
setGeneric("sparoscore_mle",
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
                    prefix = "",
                    theta_bound = 200,
                    return_se = FALSE)
               standardGeneric("sparoscore_mle"))


#' Compute censored maximum likelihood SPAROscores for matrix objects
#' @rdname sparoscore_mle
#' @export
setMethod("sparoscore_mle",
          signature(data = "matrix"),
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
                   prefix = "",
                   theta_bound = 200,
                   return_se = FALSE) {

              .augment_sparoscores_mle_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix,
                  theta_bound = theta_bound,
                  return_se = return_se)
          }
)


#' Compute censored maximum likelihood SPAROscores for sparseMatrix objects
#' @rdname sparoscore_mle
#' @export
setMethod("sparoscore_mle",
          signature(data = "sparseMatrix"),
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
                   prefix = "",
                   theta_bound = 200,
                   return_se = FALSE) {

              .augment_sparoscores_mle_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix,
                  theta_bound = theta_bound,
                  return_se = return_se)
          }
)


#' Compute censored maximum likelihood SPAROscores for DelayedMatrix objects
#' @rdname sparoscore_mle
#' @export
setMethod("sparoscore_mle",
          signature(data = "DelayedMatrix"),
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
                   prefix = "",
                   theta_bound = 200,
                   return_se = FALSE) {

              .augment_sparoscores_mle_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix,
                  theta_bound = theta_bound,
                  return_se = return_se)
          }
)


#' Compute censored maximum likelihood SPAROscores for data.frame objects
#' @rdname sparoscore_mle
#' @export
setMethod("sparoscore_mle",
          signature(data = "data.frame"),
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
                   prefix = "",
                   theta_bound = 200,
                   return_se = FALSE) {

              .augment_sparoscores_mle_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix,
                  theta_bound = theta_bound,
                  return_se = return_se)
          }
)
