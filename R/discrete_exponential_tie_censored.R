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


#' Maximum tilt of an uncensored signature
#'
#' The largest `theta` an uncensored signature of size `signature_size` can
#' produce in a ranking of `n_genes` genes. With no ties the signature fills
#' the top of the ranking, the continuous interval from `G - S` to `G`, whose
#' average rank is `(G + (G - S)) / 2`. On the unit interval that average is
#'
#' \deqn{m = 1 - \frac{S}{2G}.}
#'
#' The continuous exponential on `[0, 1]` has mean
#' `1 / (1 - exp(-theta)) - 1 / theta`. At the top of the ranking
#' `exp(-theta)` is negligible, so the mean is `1 - 1 / theta`. Setting that
#' equal to `m` gives
#'
#' \deqn{\theta = \frac{2G}{S}.}
#'
#' The smallest tilt is the negative of this value. A signature that contains
#' every gene has mean `1 / 2` and tilt 0.
#'
#' @param n_genes Size of the ranking, `G`.
#' @param signature_size Number of signature genes, `S`. A single number, or
#' a vector of sizes sharing the same `n_genes`.
#'
#' @return The maximum theoretical `theta`, the same length as
#' `signature_size`.
#'
#' @export
max_theoretical_theta <- function(n_genes, signature_size) {

    if (length(n_genes) != 1L || !is.finite(n_genes) || n_genes < 1 ||
        n_genes != as.integer(n_genes)) {
        stop("SPAROscore says: n_genes must be a single positive integer")
    }
    if (length(signature_size) == 0L || any(!is.finite(signature_size)) ||
        any(signature_size < 1) || any(signature_size != as.integer(signature_size))) {
        stop("SPAROscore says: signature_size must be positive integers")
    }
    G <- as.integer(n_genes)
    if (any(signature_size > G)) {
        stop("SPAROscore says: signature_size cannot exceed n_genes")
    }

    S <- as.integer(signature_size)
    ifelse(S >= G, 0, 2 * G / S)
}


# Genes tied across the whole ranking do not move the maximum: their block
# mean is the mean of the tilt at every theta. A signature with none left
# has a flat likelihood, and symmetry fixes theta at 0.
.fit_limit <- function(lo, hi, G, theta_bound = NULL) {

    S_bound <- sum(lo > 1L | hi < G)
    if (S_bound == 0L) return(0)
    .theta_search_limit(G, S_bound, theta_bound)
}


# Search limit for one fit. The theoretical maximum is the limit; a supplied
# theta_bound can only make the interval smaller.
.theta_search_limit <- function(G, S, theta_bound = NULL) {

    theo <- max_theoretical_theta(G, S)
    if (is.null(theta_bound)) return(theo)
    min(theta_bound, theo)
}


#' Fit the grouped-rank tilt to one signature in one cell
#'
#' Maximises the discrete log-likelihood that marginalises each gene over the
#' ranks inside its tie. No average rank is plugged in, and there is no cap.
#'
#' @noRd
.fit_grouped_tilt <- function(lo, hi, G, theta_bound = NULL) {

    S <- length(lo)
    if (S == 0L) {
        return(c(score = NA_real_, se = Inf, theta_mode = NA_real_, se_theta = Inf))
    }

    theta_bound <- .fit_limit(lo, hi, G, theta_bound)
    if (theta_bound == 0) {
        information <- (S * .geom_var(1, G, 0) -
                            sum(.geom_var(lo, hi, 0))) / G^2
        se_theta <- if (information > 0) 1 / sqrt(information) else Inf
        return(.tilt_report(0, se_theta, lo, hi, G))
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

    se_theta <- if (information > 0) 1 / sqrt(information) else Inf
    .tilt_report(theta, se_theta, lo, hi, G)
}


# Map a fitted tilt and its standard error onto the length-4 result that every
# optimiser must return. The score is the mean rank of the tilt, written on
# the unit interval, and its standard error follows by the delta method.
.tilt_report <- function(theta, se_theta, lo, hi, G) {

    rate <- theta / G
    mu_rank <- .geom_mean(1, G, rate)
    slope <- .geom_var(1, G, rate) / G^2
    usable <- is.finite(se_theta) && se_theta > 0
    c(score = (mu_rank - 0.5) / G,
      se = if (usable) slope * se_theta else Inf,
      theta_mode = theta,
      se_theta = if (usable) se_theta else Inf)
}


# Central interval of a normal approximation to the posterior, optionally
# truncated to the parameter bounds Stan enforces. The median equals the mode
# when that truncation does not cut the distribution.
.posterior_interval <- function(center, sd, prob,
                                lower_bound = -Inf, upper_bound = Inf) {

    if (!is.finite(center) || !is.finite(sd) || sd <= 0) {
        return(c(median = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    center <- min(max(center, lower_bound), upper_bound)
    alpha <- stats::pnorm(lower_bound, center, sd)
    beta <- stats::pnorm(upper_bound, center, sd)
    width <- beta - alpha
    if (!is.finite(width) || width <= 0) {
        return(c(median = center, lower = center, upper = center))
    }
    quantile_at <- function(p) {
        q <- stats::qnorm(alpha + p * width, center, sd)
        min(max(q, lower_bound), upper_bound)
    }
    tail <- (1 - prob) / 2
    c(median = quantile_at(0.5),
      lower = quantile_at(tail),
      upper = quantile_at(1 - tail))
}


.with_interval <- function(reported, prob, lower_bound = -Inf, upper_bound = Inf) {

    interval <- .posterior_interval(reported[["theta_mode"]], reported[["se_theta"]],
                                    prob, lower_bound, upper_bound)
    c(reported,
      theta_median = unname(interval[["median"]]),
      theta_lower = unname(interval[["lower"]]),
      theta_upper = unname(interval[["upper"]]))
}


# Compiled grouped-tilt model, cached the way sccomp:::load_model caches its
# CmdStan executables: one binary per package version, reused on later calls.
.load_grouped_tilt_model <- function(cache_dir = NULL, force = FALSE) {

    if (!requireNamespace("instantiate", quietly = TRUE) ||
        !instantiate::stan_cmdstan_exists()) {
        stop("SPAROscore says: optimiser = \"stan\" needs the instantiate
             package and a working CmdStan installation.")
    }
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
        stop("SPAROscore says: optimiser = \"stan\" needs the cmdstanr package.")
    }

    if (is.null(cache_dir)) {
        cache_dir <- tools::R_user_dir("SPAROscore", which = "cache")
    }
    cache_dir <- file.path(cache_dir, as.character(utils::packageVersion("SPAROscore")))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    stan_file <- system.file("stan", "grouped_tilt.stan", package = "SPAROscore")
    if (!nzchar(stan_file) || !file.exists(stan_file)) {
        stop("SPAROscore says: grouped_tilt.stan was not found. Install
             SPAROscore, or load the source tree with pkgload.")
    }

    cache_rds <- file.path(cache_dir, "grouped_tilt.rds")
    stamp_file <- file.path(cache_dir, "grouped_tilt.stamp")
    stamp <- as.character(file.info(stan_file)$mtime)

    if (!force && file.exists(cache_rds) && file.exists(stamp_file) &&
        identical(readLines(stamp_file, warn = FALSE), stamp)) {
        mod <- readRDS(cache_rds)
        exe <- tryCatch(mod$exe_file(), error = function(e) "")
        if (nzchar(exe) && file.exists(exe)) {
            return(mod)
        }
    }

    message("Compiling the grouped-tilt Stan model...")
    mod <- cmdstanr::cmdstan_model(
        stan_file = stan_file,
        exe_file = file.path(cache_dir, "grouped_tilt"),
        compile = TRUE,
        force_recompile = TRUE
    )
    saveRDS(mod, cache_rds)
    writeLines(stamp, stamp_file)
    message("Model compiled and saved to cache.")
    mod
}


# MAP of the grouped-rank likelihood under theta ~ normal(0, prior_sd).
# CmdStan's L-BFGS finds the mode. The standard error is the Laplace
# curvature at that mode: observed information of the likelihood plus the
# prior precision 1/prior_sd^2. That is the Hessian Stan would use for the
# Laplace approximation, evaluated in closed form rather than from random
# draws of it.
.fit_grouped_tilt_stan <- function(lo, hi, G, theta_bound, prior_sd, prob, mod) {

    S <- length(lo)
    if (S == 0L) {
        return(c(score = NA_real_, se = Inf, theta_mode = NA_real_, se_theta = Inf,
                 theta_median = NA_real_, theta_lower = NA_real_,
                 theta_upper = NA_real_))
    }
    if (!is.finite(prior_sd) || prior_sd <= 0) {
        stop("SPAROscore says: prior_sd must be a positive finite number")
    }

    theta_bound <- .fit_limit(lo, hi, G, theta_bound)
    if (theta_bound == 0) {
        information <- (S * .geom_var(1, G, 0) -
                            sum(.geom_var(lo, hi, 0))) / G^2 + 1 / prior_sd^2
        # The likelihood does not depend on theta, so the posterior is the
        # prior and the interval is read from that normal.
        return(.with_interval(
            .tilt_report(0, 1 / sqrt(information), lo, hi, G), prob))
    }

    stan_data <- list(
        S = S,
        G = as.integer(G),
        lo = as.integer(lo),
        hi = as.integer(hi),
        prior_sd = prior_sd,
        theta_max = theta_bound
    )
    opt <- mod$optimize(
        data = stan_data,
        init = list(list(theta = 0)),
        jacobian = FALSE,
        seed = 1L,
        refresh = 0,
        show_messages = FALSE,
        show_exceptions = FALSE
    )
    theta <- unname(opt$mle()[["theta"]])
    if (!is.finite(theta)) {
        return(c(score = NA_real_, se = Inf, theta_mode = NA_real_, se_theta = Inf,
                 theta_median = NA_real_, theta_lower = NA_real_,
                 theta_upper = NA_real_))
    }

    rate <- theta / G
    lik_info <- (S * .geom_var(1, G, rate) - sum(.geom_var(lo, hi, rate))) / G^2
    if (lik_info < 0 && lik_info > -1e-8) {
        lik_info <- 0
    }
    information <- lik_info + 1 / prior_sd^2
    se_theta <- if (information > 0) 1 / sqrt(information) else Inf

    .with_interval(.tilt_report(theta, se_theta, lo, hi, G), prob,
                   lower_bound = -theta_bound, upper_bound = theta_bound)
}


# `optimiser` is "uniroot", "stan", or a function (lo, hi, G, theta_bound).
.resolve_tilt_optimiser <- function(optimiser, prior_sd, prob,
                                    cache_stan_model) {

    if (is.function(optimiser)) {
        return(optimiser)
    }
    if (!is.character(optimiser) || length(optimiser) != 1L ||
        !optimiser %in% c("uniroot", "stan")) {
        stop("SPAROscore says: optimiser must be \"uniroot\", \"stan\", or a function")
    }
    if (optimiser == "uniroot") {
        return(.fit_grouped_tilt)
    }

    mod <- .load_grouped_tilt_model(cache_stan_model)
    function(lo, hi, G, theta_bound) {
        .fit_grouped_tilt_stan(lo, hi, G, theta_bound, prior_sd, prob, mod)
    }
}



#' Score every signature against every column
#'
#' @noRd
#'
#' @param ranks Numeric matrix of gene ranks, genes in rows and columns.
#' Rank 1 is the top (highest expression), matching [get_ranks()] but without
#' an extra cap row.
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
#' @param optimiser `"uniroot"`, `"stan"`, or a function. `"uniroot"` solves
#' the score equation. `"stan"` maximises the same likelihood plus a wide
#' normal prior with CmdStan's L-BFGS. A function must accept `lo`, `hi`, `G` and
#' `theta_bound` and return `c(score, se, theta_mode, se_theta)`. `"stan"` also
#' returns `theta_median`, `theta_lower` and `theta_upper`.
#'
#' @param prior_sd Standard deviation of the `normal(0, prior_sd)` prior used
#' when `optimiser = "stan"`.
#'
#' @param prob Posterior probability of the central interval returned when
#' `optimiser = "stan"`. The default, 0.95, is the central 95% interval.
#'
#' @param cache_stan_model Directory in which the compiled CmdStan model is
#' saved and reused. `NULL` uses the per-user package cache.
#'
#' @return A list of matrices, each with one row per column of `ranks` and one
#' column per signature: `score`, `se`, `theta_mode` and `se_theta`. When
#' `optimiser = "stan"` the list also contains `theta_median`, `theta_lower`
#' and `theta_upper`. `background` is the untilted mean position (always 0.5
#' on a complete ranking).
#'
#' @details
#' Intervals are a property of the column, so the loop ranks each column once
#' and then scores every signature against those intervals.
.compute_sparoscores_empirical_tilt <- function(ranks,
                                                signature_list,
                                                handle_missing_genes = "skip",
                                                n_items = NULL,
                                                theta_bound = NULL,
                                                optimiser = "uniroot",
                                                prior_sd = 1000,
                                                prob = 0.95,
                                                cache_stan_model = NULL) {

    if (!(handle_missing_genes %in% c("skip", "impute"))) {
        stop("SPAROscore says: Invalid value provided for handle_missing_genes.
             Limit to using 'skip' or 'impute'")
    }
    if (length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob >= 1) {
        stop("SPAROscore says: prob must be a single number between 0 and 1")
    }
    if (anyNA(ranks)) {
        stop("SPAROscore says: ranks contain NA values")
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
    fit_one <- .resolve_tilt_optimiser(optimiser, prior_sd, prob,
                                       cache_stan_model)

    fits <- lapply(cells, function(cell) {
        column <- ranks[, cell]
        intervals <- .tie_intervals(column, G)

        per_signature <- lapply(seq_len(signature_count), function(index) {
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
            fit <- fit_one(lo, hi, G, theta_bound)
            stored_names <- names(fit)
            fit <- as.numeric(fit)
            if (length(fit) < 4L) {
                stop("SPAROscore says: optimiser must return score, se, theta_mode and se_theta")
            }
            if (is.null(stored_names) || any(stored_names == "")) {
                stored_names <- c("score", "se", "theta_mode", "se_theta")
                if (length(fit) > 4L) {
                    stored_names <- c(stored_names,
                                      paste0("V", seq(5, length(fit))))
                }
            }
            names(fit) <- stored_names
            names(fit)[names(fit) == "theta"] <- "theta_mode"
            fit
        })
        do.call(rbind, per_signature)
    })

    field_names <- colnames(fits[[1]])
    out <- lapply(field_names, function(field) {
        raw <- vapply(fits, function(one) one[, field],
                      numeric(signature_count))
        if (signature_count == 1L) {
            m <- matrix(raw, ncol = 1L)
        } else {
            m <- t(raw)
        }
        dimnames(m) <- list(cells, names(signature_list))
        m
    })
    names(out) <- field_names
    out$background <- stats::setNames(rep(0.5, length(cells)), cells)
    out
}


# Shared body for every matrix-like method: turn counts into ranks when needed,
# then score. Mirrors .augment_sparoscores_mle_matrix().
.augment_sparoscores_empirical_tilt_matrix <- function(matrix_object,
                                                       signatures,
                                                       down_signatures,
                                                       data_has_ranks,
                                                       handle_missing_genes,
                                                       prefix,
                                                       theta_bound,
                                                       optimiser = "uniroot",
                                                       prior_sd = 1000,
                                                       prob = 0.95,
                                                       cache_stan_model = NULL) {

    if (!is.logical(data_has_ranks)) {
        stop("SPAROscore says: data_has_ranks should be a boolean")
    }

    if (data_has_ranks) {
        ranks <- matrix_object
    } else {
        # The grouped-rank likelihood uses the full gene ranking. No cap row
        # is appended: zeros are a tie at the bottom, not censored. Average
        # ranks keep tied genes on one shared rank, which is the block the
        # interval likelihood sums over.
        ranks <- as.matrix(MatrixGenerics::colRanks(
            -as.matrix(matrix_object),
            ties.method = "average",
            preserveShape = TRUE,
            useNames = FALSE))
        rownames(ranks) <- rownames(matrix_object)
        colnames(ranks) <- colnames(matrix_object)
    }

    score_with <- function(sets) {
        .compute_sparoscores_empirical_tilt(
            ranks = ranks,
            signature_list = sets,
            handle_missing_genes = handle_missing_genes,
            n_items = NULL,
            theta_bound = theta_bound,
            optimiser = optimiser,
            prior_sd = prior_sd,
            prob = prob,
            cache_stan_model = cache_stan_model)
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
        fitted$theta_down_mode <- down_fitted$theta_mode
        fitted$se_theta_down <- down_fitted$se_theta
        for (piece in c("median", "lower", "upper")) {
            value <- down_fitted[[paste0("theta_", piece)]]
            if (!is.null(value)) {
                fitted[[paste0("theta_down_", piece)]] <- value
            }
        }
    }

    for (element in names(fitted)) {
        if (is.matrix(fitted[[element]]) && !is.null(colnames(fitted[[element]]))) {
            colnames(fitted[[element]]) <- paste0(prefix,
                                                  colnames(fitted[[element]]))
        }
    }

    fitted
}


#' Compute SPAROscores by a discrete exponential tilt of ranks
#'
#' Compute signature scores by maximising a discrete exponential-tilt likelihood
#' in which each tied gene is marginalised over every latent rank compatible
#' with its observed tie. The average of a tie is never plugged in; the latent
#' mean inside a tie is computed afterwards, if at all, at the fitted tilt.
#'
#' `discrete_exponential_tie_censored()` is an S4 generic with methods for matrix-like
#' objects. It returns the score, its standard error, and the fitted tilt.
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
#' expression values prior to scoring, with no cap row appended. Tied genes
#' are given their average rank, so a block of `t` genes shares one rank and
#' the likelihood sums over the `t` integers that block occupies.
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
#' @param theta_bound Optional cap on the absolute value of `theta`. `NULL`
#' (default) uses [max_theoretical_theta()], the largest tilt an uncensored
#' signature of this size can produce, and the negative of that value as the
#' lower limit. A supplied number can only shrink that interval.
#'
#' @param optimiser How to maximise the grouped-rank likelihood.
#' `"uniroot"` (default) solves the score equation on the interval from
#' [max_theoretical_theta()] down to its negative. `"stan"` maximises the same
#' likelihood on that same interval, plus a
#' `normal(0, prior_sd)` prior using CmdStan's L-BFGS. The standard error is
#' the curvature of that log posterior at the mode, and `prob` sets the
#' central posterior interval around it. The compiled model is
#' cached and reused. A function may be supplied instead: it receives `lo`,
#' `hi`, `G` and `theta_bound`, and returns `c(score, se, theta_mode, se_theta)`.
#'
#' @param prior_sd Standard deviation of the normal prior on `theta` used when
#' `optimiser = "stan"`. The default, 1000, is flat across the range where the
#' likelihood has curvature, and keeps the mode finite when the likelihood is
#' flat.
#'
#' @param prob Posterior probability covered by the central interval of
#' `theta`, used when `optimiser = "stan"`. The default is 0.95. The interval
#' is the same normal approximation whose standard deviation is `se_theta`,
#' truncated to the bounds on `theta`. `theta_mode` is the posterior mode.
#' `theta_median` equals that mode when the truncation does not cut the
#' distribution.
#'
#' @param cache_stan_model Directory for the compiled CmdStan executable.
#' `NULL` uses the per-user package cache, versioned like
#' `sccomp:::load_model()`.
#'
#' @return
#' A named list. `score`, `se`, `theta_mode` and `se_theta` are matrices with one
#' row per sample, cell, or spatial domain and one column per signature.
#' `background` is the untilted mean position, one value per column of the
#' input. When `down_signatures` is supplied, `theta_down_mode` and
#' `se_theta_down` are the mode of that signature and its standard error.
#'
#' When `optimiser = "stan"`, the list also contains `theta_median`,
#' `theta_lower` and `theta_upper`: the median and the central `prob` interval
#' of the posterior. The matching `theta_down_median`,
#' `theta_down_lower` and `theta_down_upper` are included when
#' `down_signatures` is supplied.
#'
#' Scores are mean positions in the ranking, 1 being the top, on the same scale
#' as [sparoscore_mle()]. At the maximum they equal the average of the
#' conditional latent ranks `E[R_i | a_i \le R_i \le b_i, \hat\theta]`,
#' written as a position on the unit interval. `background` is the untilted
#' mean position, 0.5 on any complete ranking of `G` genes.
#'
#' `theta_mode` is the fitted posterior mode: 0 for a signature sitting where
#' an arbitrary gene sits, positive for UP and negative for DOWN. Because the
#' likelihood is the discrete exponential family on `1:G`, two cells with
#' different tie patterns remain comparable on `theta_mode`. Its standard
#' error, `se_theta`, is `1 / sqrt(I(theta_mode))`.
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
#' scores <- discrete_exponential_tie_censored(
#'   data = counts,
#'   signatures = c("gene1", "gene2", "gene3")
#' )
#'
#' scores$score
#' scores$se
#' scores$theta_mode
#' scores$se_theta
#' scores$background
#'
#' @export
setGeneric("discrete_exponential_tie_censored",
           function(data,
                    signatures,
                    down_signatures = NULL,
                    data_has_ranks = FALSE,
                    handle_missing_genes = "skip",
                    prefix = "",
                    theta_bound = NULL,
                    optimiser = "uniroot",
                    prior_sd = 1000,
                    prob = 0.95,
                    cache_stan_model = NULL)
               standardGeneric("discrete_exponential_tie_censored"))


#' Compute empirical tilt SPAROscores for matrix objects
#' @rdname discrete_exponential_tie_censored
#' @export
setMethod("discrete_exponential_tie_censored",
          signature(data = "matrix"),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   handle_missing_genes = "skip",
                   prefix = "",
                   theta_bound = NULL,
                   optimiser = "uniroot",
                   prior_sd = 1000,
                   prob = 0.95,
                   cache_stan_model = NULL) {

              .augment_sparoscores_empirical_tilt_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix,
                  theta_bound = theta_bound,
                  optimiser = optimiser,
                  prior_sd = prior_sd,
                  prob = prob,
                  cache_stan_model = cache_stan_model)
          }
)


#' Compute empirical tilt SPAROscores for sparseMatrix objects
#' @rdname discrete_exponential_tie_censored
#' @export
setMethod("discrete_exponential_tie_censored",
          signature(data = "sparseMatrix"),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   handle_missing_genes = "skip",
                   prefix = "",
                   theta_bound = NULL,
                   optimiser = "uniroot",
                   prior_sd = 1000,
                   prob = 0.95,
                   cache_stan_model = NULL) {

              .augment_sparoscores_empirical_tilt_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix,
                  theta_bound = theta_bound,
                  optimiser = optimiser,
                  prior_sd = prior_sd,
                  prob = prob,
                  cache_stan_model = cache_stan_model)
          }
)


#' Compute empirical tilt SPAROscores for DelayedMatrix objects
#' @rdname discrete_exponential_tie_censored
#' @export
setMethod("discrete_exponential_tie_censored",
          signature(data = "DelayedMatrix"),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   handle_missing_genes = "skip",
                   prefix = "",
                   theta_bound = NULL,
                   optimiser = "uniroot",
                   prior_sd = 1000,
                   prob = 0.95,
                   cache_stan_model = NULL) {

              .augment_sparoscores_empirical_tilt_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix,
                  theta_bound = theta_bound,
                  optimiser = optimiser,
                  prior_sd = prior_sd,
                  prob = prob,
                  cache_stan_model = cache_stan_model)
          }
)


#' Compute empirical tilt SPAROscores for data.frame objects
#' @rdname discrete_exponential_tie_censored
#' @export
setMethod("discrete_exponential_tie_censored",
          signature(data = "data.frame"),
          function(data,
                   signatures,
                   down_signatures = NULL,
                   data_has_ranks = FALSE,
                   handle_missing_genes = "skip",
                   prefix = "",
                   theta_bound = NULL,
                   optimiser = "uniroot",
                   prior_sd = 1000,
                   prob = 0.95,
                   cache_stan_model = NULL) {

              .augment_sparoscores_empirical_tilt_matrix(
                  matrix_object = data,
                  signatures = signatures,
                  down_signatures = down_signatures,
                  data_has_ranks = data_has_ranks,
                  handle_missing_genes = handle_missing_genes,
                  prefix = prefix,
                  theta_bound = theta_bound,
                  optimiser = optimiser,
                  prior_sd = prior_sd,
                  prob = prob,
                  cache_stan_model = cache_stan_model)
          }
)
