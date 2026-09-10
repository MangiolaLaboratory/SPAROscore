# ---------------------------------------------------------------------------
# Discrete exponential tilt with ties as grouped ranks.
#
#   P(R = r | lambda) = exp(lambda * r / G) / Z(lambda),   r = 1, ..., G
#
# Rank 1 is the bottom of the ranking and G the top, so lambda > 0 is UP.
# A tie is not a point mass at the average rank. If gene i is only known to
# lie in [a_i, b_i], its contribution is the marginal
#
#   L_i(lambda) = sum_{r = a_i}^{b_i} P(R = r | lambda)
#
# and the log-likelihood factorises across genes:
#
#   ell(lambda) = sum_i log(sum_{r = a_i}^{b_i} e^{lambda r / G})
#               - S log Z(lambda)
#
# Untied genes have a_i = b_i = r_i and recover the exact-rank likelihood.
# The latent mean E[R_i | tie, lambda] is computed afterwards, if at all; it
# does not enter the fit.
# ---------------------------------------------------------------------------

# log(e^x - 1) for x > 0, without overflowing.
.geom_log_expm1 <- function(x) {
    ifelse(x > 18, x, log(expm1(x)))
}

# log sum_{r=lo}^{hi} exp(rate * r). Vectorised over lo, hi; rate is scalar.
.geom_logsum <- function(lo, hi, rate) {

    t <- hi - lo + 1
    out <- numeric(length(lo))
    tiny <- abs(rate) < 1e-10
    out[tiny] <- log(t[tiny]) + rate * (lo[tiny] + hi[tiny]) / 2
    pos <- !tiny & rate > 0
    neg <- !tiny & rate < 0
    if (any(pos)) {
        out[pos] <- rate * lo[pos] + .geom_log_expm1(rate * t[pos]) -
            .geom_log_expm1(rate)
    }
    if (any(neg)) {
        out[neg] <- rate * hi[neg] + .geom_log_expm1(-rate * t[neg]) -
            .geom_log_expm1(-rate)
    }
    out
}

# E[K] for K = 0, ..., t-1 with P(K) proportional to exp(rate * K).
.geom_ek <- function(t, rate) {

    if (t == 1L) return(0)
    if (abs(rate) < 1e-10) return((t - 1) / 2)
    if (rate > 0) return((t - 1) - .geom_ek(t, -rate))
    q <- exp(rate)
    if (rate * t < -37) return(q / (1 - q))
    q / (1 - q) - t * q^t / (1 - q^t)
}

# E[R | lo <= R <= hi] under P(R) proportional to exp(rate * r).
.geom_mean <- function(lo, hi, rate) {

    t <- hi - lo + 1
    out <- (lo + hi) / 2
    one <- t == 1
    out[one] <- lo[one]
    rest <- which(!one)
    if (length(rest) == 0L) return(out)
    out[rest] <- lo[rest] + vapply(rest, function(i) {
        .geom_ek(t[i], rate)
    }, numeric(1))
    out
}

# Var(R | lo <= R <= hi). Used only for the observed information.
.geom_var <- function(lo, hi, rate) {

    t <- hi - lo + 1
    out <- (t^2 - 1) / 12
    one <- t == 1
    out[one] <- 0
    rest <- which(!one & abs(rate) >= 1e-10)
    if (length(rest) == 0L) return(out)
    out[rest] <- vapply(rest, function(i) {
        ti <- t[i]
        # Var(K) = q/(1-q)^2 - t^2 q^t/(1-q^t)^2  for rate < 0, else reverse.
        r <- rate
        if (r > 0) r <- -rate
        q <- exp(r)
        vk <- if (r * ti < -37) {
            q / (1 - q)^2
        } else {
            q / (1 - q)^2 - ti^2 * q^ti / (1 - q^ti)^2
        }
        vk
    }, numeric(1))
    out
}


#' Tie intervals in ascending ranks, 1 at the bottom and G at the top.
#'
#' Recovers each gene's compatible integer interval from a column of ranks
#' (package descending ranks converted internally). Genes that share a rank
#' occupy consecutive integers covering 1..G, so `"average"`, `"min"` and
#' `"max"` all yield the same intervals; `"random"` has already broken the ties.
#'
#' @noRd
.tie_intervals <- function(gene_ranks, G) {

    asc <- unname((G + 1) - gene_ranks)
    blocks <- rle(sort(asc))
    t <- blocks$lengths
    hi <- cumsum(t)
    lo <- hi - t + 1
    j <- match(asc, blocks$values)
    list(lo = lo[j], hi = hi[j], G = G)
}


#' Fit the grouped-rank tilt to one signature in one cell
#'
#' Maximises the discrete log-likelihood that marginalises each gene over the
#' ranks inside its tie. No average rank is plugged in, and there is no cap.
#'
#' @noRd
.fit_grouped_tilt <- function(lo, hi, G, theta_bound = 200) {

    S <- length(lo)
    if (S == 0L) {
        return(c(score = NA_real_, se = Inf, theta = NA_real_, se_theta = Inf))
    }

    # Score of
    #   ell(theta) = sum log(sum_{r=a}^b exp(theta r / G)) - S log Z.
    # d ell / d theta = (sum E[R|block] - S mu) / G, strictly decreasing.
    score_equation <- function(theta) {
        rate <- theta / G
        (sum(.geom_mean(lo, hi, rate)) - S * .geom_mean(1, G, rate)) / G
    }

    lower <- score_equation(-theta_bound)
    upper <- score_equation(theta_bound)

    theta <- if (lower <= 0) {
        -theta_bound
    } else if (upper >= 0) {
        theta_bound
    } else {
        stats::uniroot(score_equation, c(-theta_bound, theta_bound),
                       f.lower = lower, f.upper = upper, tol = 1e-9)$root
    }

    rate <- theta / G
    var_full <- .geom_var(1, G, rate)
    var_hidden <- sum(.geom_var(lo, hi, rate))
    # information for theta: (S Var(R) - sum Var(R|block)) / G^2
    information <- (S * var_full - var_hidden) / G^2

    mu_rank <- .geom_mean(1, G, rate)
    slope <- var_full / G^2
    usable <- information > 0

    c(score = (mu_rank - 0.5) / G,
      se = if (usable) slope / sqrt(information) else Inf,
      theta = theta,
      se_theta = if (usable) 1 / sqrt(information) else Inf)
}



#' Score every signature against every column
#'
#' @noRd
#'
#' @param ranks Numeric matrix of gene ranks, genes in rows and columns.
#' Rank 1 is the top (highest expression), matching [get_ranks()] but without
#' an extra cap row.
#'
#' @param rank_caps Named numeric vector with one rank cap per column.
#'
#' @param signature_list Named list of character vectors.
#'
#' @param handle_missing_genes Either "skip" or "impute".
#'
#' @param n_items Size of the ranking universe, `G`. Defaults to the number of
#' genes. Ties are intervals inside `1:G`; there is no extra cap row.
#'
#' @param theta_bound Half-width of the search interval for `theta`.
#'
#' @param max_atoms Ignored. Kept so the argument list matches older callers.
#'
#' @return A list of matrices `score`, `se`, `theta` and `se_theta`, each with
#' one row per column of `ranks` and one column per signature, plus
#' `background`, the untilted mean position (always 0.5 on a complete ranking).
#'
#' @details
#' Intervals are a property of the column, so the loop ranks each column once
#' and then scores every signature against those intervals.
.compute_sparoscores_empirical_tilt <- function(ranks,
                                                rank_caps,
                                                signature_list,
                                                handle_missing_genes = "skip",
                                                n_items = NULL,
                                                theta_bound = 200,
                                                max_atoms = 1024L) {

    if (!(handle_missing_genes %in% c("skip", "impute"))) {
        stop("SPAROscore says: Invalid value provided for handle_missing_genes.
             Limit to using 'skip' or 'impute'")
    }
    if (anyNA(ranks)) {
        stop("SPAROscore says: ranks contain NA values")
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

    # Ties live in the column, so every gene is ranked once and then each
    # signature is scored against those intervals.
    ranks <- as.matrix(ranks)
    if (!all(is.finite(ranks))) {
        stop("SPAROscore says: ranks must be finite numeric values")
    }

    if (is.null(n_items)) {
        n_items <- nrow(ranks)
    }

    validated <- lapply(signature_list, validate_signature, rownames(ranks))
    signature_rows <- lapply(validated, function(x) {
        match(x$valid_genes, rownames(ranks))
    })

    extra_bottom <- vapply(validated, function(x) {
        if (handle_missing_genes == "impute") length(x$missing_genes) else 0L
    }, integer(1))

    cells <- colnames(ranks)
    signature_count <- length(signature_list)
    G <- as.integer(n_items)

    fitted <- vapply(cells, function(cell) {
        column <- ranks[, cell]
        intervals <- .tie_intervals(column, G)

        per_signature <- vapply(seq_len(signature_count), function(index) {
            rows <- signature_rows[[index]]
            lo <- intervals$lo[rows]
            hi <- intervals$hi[rows]
            extra <- extra_bottom[[index]]
            if (extra > 0L) {
                # genes absent from the matrix contribute the whole ranking
                # and cancel in the likelihood; they only dilute the information.
                lo <- c(lo, rep(1, extra))
                hi <- c(hi, rep(G, extra))
            }
            .fit_grouped_tilt(lo, hi, G, theta_bound)
        }, numeric(4))

        c(as.vector(per_signature), 0.5)
    }, numeric(4 * signature_count + 1))

    # Rows arrive as (score, se, theta, se_theta) per signature followed by the
    # background, and each output wants cells down the rows.
    unstack <- function(offset) {
        out <- t(fitted[seq(offset, by = 4, length.out = signature_count), ,
                        drop = FALSE])
        dimnames(out) <- list(cells, names(signature_list))
        out
    }

    list(score = unstack(1L),
         se = unstack(2L),
         theta = unstack(3L),
         se_theta = unstack(4L),
         background = stats::setNames(fitted[4 * signature_count + 1, ], cells))
}


# Shared body for every matrix-like method: turn counts into ranks when needed,
# then score. Mirrors .augment_sparoscores_mle_matrix().
.augment_sparoscores_empirical_tilt_matrix <- function(matrix_object,
                                                       signatures,
                                                       down_signatures,
                                                       data_has_ranks,
                                                       count_caps,
                                                       rank_caps,
                                                       handle_ties,
                                                       handle_missing_genes,
                                                       prefix,
                                                       theta_bound,
                                                       max_atoms,
                                                       return_se) {

    invisible(count_caps)
    invisible(max_atoms)

    if (!is.logical(data_has_ranks)) {
        stop("SPAROscore says: data_has_ranks should be a boolean")
    }

    if (data_has_ranks) {
        if (is.null(rank_caps)) {
            stop("SPAROscore says: rank_caps needed when data_has_ranks = TRUE")
        }
        ranks <- matrix_object
    } else {
        # The grouped-rank likelihood uses the full gene ranking. No cap row
        # is appended: zeros are a tie at the bottom, not censored.
        ranks <- as.matrix(MatrixGenerics::colRanks(
            -as.matrix(matrix_object),
            ties.method = handle_ties,
            preserveShape = TRUE,
            useNames = FALSE))
        rownames(ranks) <- rownames(matrix_object)
        colnames(ranks) <- colnames(matrix_object)
        if (is.null(rank_caps)) {
            rank_caps <- stats::setNames(
                rep(nrow(ranks) + 1L, ncol(ranks)),
                colnames(ranks))
        }
    }

    score_with <- function(sets) {
        .compute_sparoscores_empirical_tilt(
            ranks = ranks,
            rank_caps = rank_caps,
            signature_list = sets,
            handle_missing_genes = handle_missing_genes,
            n_items = NULL,
            theta_bound = theta_bound,
            max_atoms = max_atoms)
    }

    signature_list <- .as_signature_list(signatures, "SPAROscore_tilt")
    fitted <- score_with(signature_list)

    # Both scores are means on the same 0-1 scale, so the contrast is their
    # difference and the errors add in quadrature. The two tilts are reported
    # separately because a difference of tilts is not itself a tilt.
    if (!is.null(down_signatures)) {
        down_list <- .as_signature_list(down_signatures, "SPAROscore_tilt")
        if (!identical(names(signature_list), names(down_list))) {
            stop("SPAROscore says: The names of gene sets in
                      down_signatures should match that of signatures")
        }
        down_fitted <- score_with(down_list)
        fitted$score <- fitted$score - down_fitted$score
        fitted$se <- sqrt(fitted$se^2 + down_fitted$se^2)
        fitted$theta_down <- down_fitted$theta
        fitted$se_theta_down <- down_fitted$se_theta
    }

    for (element in c("score", "se", "theta", "se_theta", "theta_down",
                      "se_theta_down")) {
        if (!is.null(fitted[[element]])) {
            colnames(fitted[[element]]) <- paste0(prefix,
                                                  colnames(fitted[[element]]))
        }
    }

    if (return_se) {
        return(fitted)
    }

    scores <- fitted$score
    attr(scores, "standard_error") <- fitted$se
    attr(scores, "theta") <- fitted$theta
    attr(scores, "se_theta") <- fitted$se_theta
    attr(scores, "background") <- fitted$background
    scores
}


#' Compute SPAROscores by a discrete exponential tilt of ranks
#'
#' Compute signature scores by maximising a discrete exponential-tilt likelihood
#' in which each tied gene is marginalised over every latent rank compatible
#' with its observed tie. The average of a tie is never plugged in; the latent
#' mean inside a tie is computed afterwards, if at all, at the fitted tilt.
#'
#' `sparoscore_empirical_tilt()` is an S4 generic with methods for matrix-like
#' objects. It takes the same arguments as [sparoscore()] and
#' [sparoscore_mle()] and returns a matrix of the same shape, so calls can be
#' swapped over without changing surrounding code.
#'
#' It differs from [sparoscore_mle()] in the observation model. That function
#' tilts a continuous uniform and treats a rank as a point (or as censored
#' below a cap). This one tilts the discrete ranks `1, ..., G` and, for a gene
#' observed only as a tie spanning `[a_i, b_i]`, sums the tilted mass over that
#' interval. Untied genes have `a_i = b_i` and recover the exact-rank
#' likelihood. There is no rank cap: zeros are the bottom tie block.
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
#' expression values prior to scoring, with no cap row appended.
#'
#' @param assay Ignored for matrix-like inputs. Present so that the argument
#' list matches [sparoscore()].
#'
#' @param layer Ignored for matrix-like inputs. Present so that the argument
#' list matches [sparoscore()].
#'
#' @param count_caps Ignored. The grouped-rank likelihood does not censor at a
#' count threshold. Kept so the argument list matches [sparoscore()].
#'
#' @param rank_caps Required when `data_has_ranks = TRUE` so the argument list
#' matches [sparoscore()]; otherwise filled with a dummy past the bottom of the
#' ranking. The fit itself does not censor at these values.
#'
#' @param handle_ties Character string specifying how tied expression values
#' should be ranked. Passed to `MatrixGenerics::colRanks(ties.method = ...)`.
#' Supported values are `"average"` (default here), `"min"`, `"max"` and
#' `"random"`.
#'
#' `"average"`, `"min"` and `"max"` all recover the same integer intervals: a
#' block of `t` genes occupies `t` consecutive ranks, and the likelihood sums
#' over those ranks. `"random"` has already broken the ties, so each gene is
#' treated as an exact rank.
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from the ranks should be handled. Supported options are:
#' \describe{
#' \item{"skip"}{Exclude missing genes from the score calculation (default).}
#' \item{"impute"}{Treat missing genes as known only to lie somewhere in
#' `1:G`. That contribution cancels in the likelihood, so the point estimate
#' is unchanged and the information is diluted.}
#' }
#'
#' @param prefix Character string prepended to the names of the returned score
#' columns. Defaults to `""`.
#'
#' @param theta_bound Half-width of the interval searched for the tilt
#' parameter. The reportable score range follows from it, roughly
#' `1/theta_bound` to `1 - 1/theta_bound`.
#'
#' @param max_atoms Ignored. Kept so the argument list matches older callers.
#'
#' @param return_se Logical. If `FALSE` (default), return a matrix of scores
#' carrying the standard errors, tilts and background means as attributes, so
#' the result is shaped exactly like the output of [sparoscore()]. If `TRUE`,
#' return them as a named list instead.
#'
#' @return
#' A numeric matrix of scores with one row per sample, cell, or spatial domain
#' and one column per signature. Its `"standard_error"`, `"theta"` and
#' `"se_theta"` attributes are matrices of the same shape, and its
#' `"background"` attribute is one value per column of the input.
#'
#' When `return_se = TRUE`, a named list of those same pieces: `score`, `se`,
#' `theta`, `se_theta`, `background`, and `theta_down` and `se_theta_down` when
#' `down_signatures` is supplied.
#'
#' Scores are mean positions in the ranking, 1 being the top, on the same scale
#' as [sparoscore_mle()]. At the maximum they equal the average of the
#' conditional latent ranks `E[R_i | a_i \le R_i \le b_i, \hat\theta]`,
#' written as a position on the unit interval. `background` is the untilted
#' mean position, 0.5 on any complete ranking of `G` genes.
#'
#' `theta` is the fitted tilt: 0 for a signature sitting where an arbitrary gene
#' sits, positive for UP and negative for DOWN. Because the likelihood is the
#' discrete exponential family on `1:G`, two cells with different tie patterns
#' remain comparable on `theta`. Its standard error, `se_theta`, is
#' `1 / sqrt(I(theta))`.
#'
#' A score of `NA` with an infinite standard error marks a column in which no
#' signature gene was present to score.
#'
#' @details
#' The latent ranking is discrete,
#'
#' \deqn{P(R = r \mid \theta) = e^{\theta r / G} / Z(\theta),
#'       \qquad r = 1, \ldots, G.}
#'
#' Rank 1 is the bottom of the ranking and `G` the top, so `\theta > 0` is UP.
#' If gene `i` is observed in a tied block spanning ranks `[a_i, b_i]`, its
#' contribution is the marginal
#'
#' \deqn{L_i(\theta) = \sum_{r = a_i}^{b_i} P(R = r \mid \theta)
#'       = \frac{\sum_{r = a_i}^{b_i} e^{\theta r / G}}{Z(\theta)}.}
#'
#' The genes factorise given `\theta`, so the log-likelihood is
#'
#' \deqn{\ell(\theta) = \sum_{i = 1}^{S}
#'       \log\Bigl(\sum_{r = a_i}^{b_i} e^{\theta r / G}\Bigr)
#'       - S \log Z(\theta).}
#'
#' Combinations of latent ranks across genes are never enumerated. After
#' `\hat\theta` is obtained, the latent mean inside a tie is
#'
#' \deqn{E[R_i \mid a_i \le R_i \le b_i, \hat\theta]
#'       = \frac{\sum_{r = a_i}^{b_i} r\, e^{\hat\theta r / G}}
#'         {\sum_{r = a_i}^{b_i} e^{\hat\theta r / G}}.}
#'
#' That posterior mean is not an input to the fit. Plugging the average of a
#' tie into an exact-rank likelihood is the `\theta = 0` special case of this
#' formula, and is biased once `\theta \ne 0`.
#'
#' The score equation equates the average of those conditional means to the
#' model's mean rank. The observed information subtracts the residual
#' within-block variances, so a gene in a singleton contributes a full slot and
#' a gene in a giant zero block contributes almost nothing.
#'
#' @seealso [sparoscore_mle()], [sparoscore()], [get_ranks()]
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
#' scores <- sparoscore_empirical_tilt(
#'   data = counts,
#'   signatures = c("gene1", "gene2", "gene3")
#' )
#'
#' # the uncertainty, the tilt and the neutral point travel with the score
#' attr(scores, "standard_error")
#' attr(scores, "theta")
#' attr(scores, "se_theta")
#' attr(scores, "background")
#'
#' # or ask for all of them explicitly
#' fitted <- sparoscore_empirical_tilt(
#'   data = counts,
#'   signatures = c("gene1", "gene2", "gene3"),
#'   return_se = TRUE
#' )
#'
#' @export
setGeneric("sparoscore_empirical_tilt",
           function(data,
                    signatures,
                    down_signatures = NULL,
                    data_has_ranks = FALSE,
                    assay = "RNA",
                    layer = "counts",
                    count_caps = NULL,
                    rank_caps = NULL,
                    handle_ties = "average",
                    handle_missing_genes = "skip",
                    prefix = "",
                    theta_bound = 200,
                    max_atoms = 1024L,
                    return_se = FALSE)
               standardGeneric("sparoscore_empirical_tilt"))


#' Compute empirical tilt SPAROscores for matrix objects
#' @rdname sparoscore_empirical_tilt
#' @export
setMethod("sparoscore_empirical_tilt",
          signature(data = "matrix"),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = NULL,
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "average",
                   handle_missing_genes = "skip",
                   prefix = "",
                   theta_bound = 200,
                   max_atoms = 1024L,
                   return_se = FALSE) {

              .augment_sparoscores_empirical_tilt_matrix(
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
                  max_atoms = max_atoms,
                  return_se = return_se)
          }
)


#' Compute empirical tilt SPAROscores for sparseMatrix objects
#' @rdname sparoscore_empirical_tilt
#' @export
setMethod("sparoscore_empirical_tilt",
          signature(data = "sparseMatrix"),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = NULL,
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "average",
                   handle_missing_genes = "skip",
                   prefix = "",
                   theta_bound = 200,
                   max_atoms = 1024L,
                   return_se = FALSE) {

              .augment_sparoscores_empirical_tilt_matrix(
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
                  max_atoms = max_atoms,
                  return_se = return_se)
          }
)


#' Compute empirical tilt SPAROscores for DelayedMatrix objects
#' @rdname sparoscore_empirical_tilt
#' @export
setMethod("sparoscore_empirical_tilt",
          signature(data = "DelayedMatrix"),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = NULL,
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "average",
                   handle_missing_genes = "skip",
                   prefix = "",
                   theta_bound = 200,
                   max_atoms = 1024L,
                   return_se = FALSE) {

              .augment_sparoscores_empirical_tilt_matrix(
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
                  max_atoms = max_atoms,
                  return_se = return_se)
          }
)


#' Compute empirical tilt SPAROscores for data.frame objects
#' @rdname sparoscore_empirical_tilt
#' @export
setMethod("sparoscore_empirical_tilt",
          signature(data = "data.frame"),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   assay = NULL,
                   layer = NULL,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "average",
                   handle_missing_genes = "skip",
                   prefix = "",
                   theta_bound = 200,
                   max_atoms = 1024L,
                   return_se = FALSE) {

              .augment_sparoscores_empirical_tilt_matrix(
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
                  max_atoms = max_atoms,
                  return_se = return_se)
          }
)
