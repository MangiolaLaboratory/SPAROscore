# helper testing function
make_tilt_counts <- function(genes = 200L, cells = 4L, seed = 11L) {
    set.seed(seed)
    matrix(
        stats::rpois(genes * cells, lambda = 3),
        nrow = genes,
        dimnames = list(
            paste0("gene", seq_len(genes)),
            paste0("cell", seq_len(cells))
        )
    )
}


test_that("discrete_exponential_tie_censored returns the fit and its uncertainty", {

    counts <- make_tilt_counts()

    res <- discrete_exponential_tie_censored(
        data = counts,
        signatures = paste0("gene", 1:20)
    )

    expect_named(res, c("score", "se", "theta", "se_theta", "background"))
    expect_true(is.matrix(res$score))
    expect_equal(nrow(res$score), ncol(counts))
    expect_equal(colnames(res$score), "SPAROscore_tilt")
    expect_true(all(res$score > 0 & res$score < 1))
    expect_equal(dim(res$se), dim(res$score))
    expect_equal(dim(res$theta), dim(res$score))
    expect_equal(dim(res$se_theta), dim(res$score))
    expect_equal(length(res$background), ncol(counts))
    expect_true(all(res$se > 0))
    expect_true(all(res$se < res$se_theta))
})


test_that("discrete_exponential_tie_censored agrees across matrix-like input classes", {

    skip_if_not_installed("Matrix")

    counts <- make_tilt_counts()
    signature <- paste0("gene", 1:20)

    dense <- discrete_exponential_tie_censored(counts, signature)
    sparse <- discrete_exponential_tie_censored(Matrix::Matrix(counts, sparse = TRUE),
                                       signature)
    frame <- discrete_exponential_tie_censored(as.data.frame(counts), signature)

    expect_equal(as.vector(dense$score), as.vector(sparse$score))
    expect_equal(as.vector(dense$theta), as.vector(sparse$theta))
    expect_equal(as.vector(dense$score), as.vector(frame$score))
})


test_that("discrete_exponential_tie_censored accepts pre-computed ranks", {

    counts <- make_tilt_counts()
    signature <- paste0("gene", 1:20)

    ranks <- as.matrix(MatrixGenerics::colRanks(
        -counts, ties.method = "average", preserveShape = TRUE,
        useNames = FALSE))
    dimnames(ranks) <- dimnames(counts)

    from_counts <- discrete_exponential_tie_censored(counts, signature)
    from_ranks <- discrete_exponential_tie_censored(
        data = ranks,
        signatures = signature,
        data_has_ranks = TRUE
    )

    expect_equal(as.vector(from_counts$score), as.vector(from_ranks$score))
})


test_that("discrete_exponential_tie_censored handles multiple signatures and a prefix", {

    counts <- make_tilt_counts()

    res <- discrete_exponential_tie_censored(
        data = counts,
        signatures = list(
            SignatureA = paste0("gene", 1:20),
            SignatureB = paste0("gene", 50:80)
        ),
        prefix = "tilt_"
    )

    expect_equal(colnames(res$score), c("tilt_SignatureA", "tilt_SignatureB"))
    expect_equal(colnames(res$se), colnames(res$score))

    # one signature at a time must give the same numbers as both at once
    alone <- discrete_exponential_tie_censored(counts, paste0("gene", 50:80))
    expect_equal(unname(res$score[, 2]), as.vector(alone$score))
    expect_error(
        discrete_exponential_tie_censored(counts, list(paste0("gene", 1:20))),
        "named"
    )
})


test_that("discrete_exponential_tie_censored subtracts down_signatures", {

    counts <- make_tilt_counts()
    up <- list(Signature = paste0("gene", 1:20))
    down <- list(Signature = paste0("gene", 50:80))

    combined <- discrete_exponential_tie_censored(counts, up, down_signatures = down)
    up_only <- discrete_exponential_tie_censored(counts, up)
    down_only <- discrete_exponential_tie_censored(counts, down)

    expect_equal(combined$score, up_only$score - down_only$score)
    expect_equal(combined$se, sqrt(up_only$se^2 + down_only$se^2))
    expect_equal(combined$theta_down, down_only$theta)
    expect_equal(combined$se_theta_down, down_only$se_theta)
    expect_error(
        discrete_exponential_tie_censored(counts, up,
                                  down_signatures = list(Other = "gene1")),
        "should match"
    )
})


test_that("discrete_exponential_tie_censored scores zeros instead of dropping them", {

    counts <- make_tilt_counts()
    signature <- paste0("gene", 1:20)

    # push every signature gene to the bottom of the first cell: they sit in
    # the zero block and are still part of the likelihood
    counts[signature, 1] <- 0L
    counts[setdiff(rownames(counts), signature), 1] <- 100L

    fitted <- discrete_exponential_tie_censored(counts, signature)

    expect_true(is.finite(fitted$score[1, 1]))
    expect_true(fitted$theta[1, 1] < 0)
    expect_true(all(is.finite(fitted$score)))
})


test_that("discrete_exponential_tie_censored imputes missing genes as uninformative", {

    counts <- make_tilt_counts()
    present <- paste0("gene", 1:20)
    padded <- c(present, "absent1", "absent2", "absent3")

    skipped <- suppressWarnings(
        discrete_exponential_tie_censored(counts, padded,
                                  handle_missing_genes = "skip"))
    imputed <- suppressWarnings(
        discrete_exponential_tie_censored(counts, padded,
                                  handle_missing_genes = "impute"))

    expect_equal(as.vector(skipped$score),
                 as.vector(discrete_exponential_tie_censored(counts, present)$score))
    expect_equal(as.vector(imputed$theta), as.vector(skipped$theta),
                 tolerance = 1e-8)
    expect_true(all(imputed$se_theta >= skipped$se_theta - 1e-12))
    expect_error(
        suppressWarnings(
            discrete_exponential_tie_censored(counts, padded,
                                      handle_missing_genes = "nonsense")),
        "handle_missing_genes"
    )
})


test_that("geometric sums match a direct calculation", {

    lo <- c(1, 10, 50)
    hi <- c(1, 12, 80)
    rate <- 0.02

    direct <- vapply(seq_along(lo), function(i) {
        log(sum(exp(rate * seq(lo[i], hi[i]))))
    }, numeric(1))

    expect_equal(SPAROscore:::.geom_logsum(lo, hi, rate), direct,
                 tolerance = 1e-10)
    expect_equal(
        SPAROscore:::.geom_mean(10, 12, rate),
        sum(seq(10, 12) * exp(rate * seq(10, 12))) /
            sum(exp(rate * seq(10, 12)))
    )

    # large |rate| must not overflow
    expect_true(is.finite(SPAROscore:::.geom_logsum(1, 10000, 1)))
    expect_true(is.finite(SPAROscore:::.geom_logsum(1, 10000, -1)))
})


test_that("an untied fit is exact moment matching", {

    set.seed(8)
    genes <- 300L
    cells <- 3L
    counts <- sapply(seq_len(cells), function(i) sample.int(genes))
    dimnames(counts) <- list(paste0("gene", seq_len(genes)),
                             paste0("cell", seq_len(cells)))
    signature <- paste0("gene", 1:40)

    ranks <- as.matrix(MatrixGenerics::colRanks(
        -counts, ties.method = "average", preserveShape = TRUE,
        useNames = FALSE))
    dimnames(ranks) <- dimnames(counts)

    fitted <- discrete_exponential_tie_censored(
        data = ranks,
        signatures = signature,
        data_has_ranks = TRUE
    )

    # unique ranks: the grouped likelihood reduces to the exact-rank one,
    # so the fitted mean is the sample mean of the signature's ranks
    asc <- (genes + 1) - ranks[signature, ]
    truth <- colMeans((asc - 0.5) / genes)
    expect_equal(as.vector(fitted$score), unname(truth), tolerance = 1e-7)
    expect_equal(unname(fitted$background), rep(0.5, cells), tolerance = 1e-12)
})


test_that("theta is zero for a signature sitting at the background", {

    counts <- make_tilt_counts(genes = 300L, cells = 3L)

    fitted <- discrete_exponential_tie_censored(
        data = counts,
        signatures = rownames(counts)
    )

    expect_equal(as.vector(fitted$theta), rep(0, ncol(counts)),
                 tolerance = 1e-6)
    expect_equal(as.vector(fitted$score), unname(fitted$background),
                 tolerance = 1e-9)
})


test_that("a tie is an interval, not an atom at the average rank", {

    # 300 genes at zero, 200 at one, the rest distinct. The zero block occupies
    # ascending ranks 1:300, the ones occupy 301:500.
    counts <- c(rep(0, 300), rep(1, 200), 2:501)
    n_items <- length(counts)
    ranks <- rank(-counts, ties.method = "average")

    intervals <- SPAROscore:::.tie_intervals(ranks, n_items)

    expect_equal(unique(intervals$lo[1:300]), 1)
    expect_equal(unique(intervals$hi[1:300]), 300)
    expect_equal(unique(intervals$lo[301:500]), 301)
    expect_equal(unique(intervals$hi[301:500]), 500)

    # min and max recover the same intervals
    intervals_min <- SPAROscore:::.tie_intervals(
        rank(-counts, ties.method = "min"), n_items)
    expect_equal(intervals$lo, intervals_min$lo)
    expect_equal(intervals$hi, intervals_min$hi)
})


test_that("the grouped likelihood is not the midpoint plug-in", {

    lo <- c(1, 8000)
    hi <- c(6800, 8000)
    G <- 10000
    theta <- 5
    rate <- theta / G
    S <- length(lo)

    grouped <- sum(SPAROscore:::.geom_logsum(lo, hi, rate)) -
        S * SPAROscore:::.geom_logsum(1, G, rate)
    mid <- (lo + hi) / 2
    plugged <- sum(rate * mid) - S * SPAROscore:::.geom_logsum(1, G, rate)

    expect_false(isTRUE(all.equal(grouped, plugged)))

    # after fitting an UP-ish signature, the latent mean in the zero block sits
    # above the midpoint
    fitted <- SPAROscore:::.fit_grouped_tilt(lo, hi, G)
    latent <- SPAROscore:::.geom_mean(lo, hi, fitted[["theta"]] / G)
    expect_gt(fitted[["theta"]], 0)
    expect_gt(latent[1], (lo[1] + hi[1]) / 2)
    expect_equal(latent[2], 8000)
})


test_that("average ranks recover the same intervals as min and max", {

    counts <- make_tilt_counts(genes = 400L, cells = 2L)
    signature <- paste0("gene", 1:50)
    G <- nrow(counts)

    intervals_of <- function(ties) {
        ranks <- rank(-counts[, 1], ties.method = ties)
        SPAROscore:::.tie_intervals(ranks, G)
    }

    average <- intervals_of("average")
    expect_equal(intervals_of("min"), average)
    expect_equal(intervals_of("max"), average)
    expect_equal(
        unname(discrete_exponential_tie_censored(counts, signature)$background),
        rep(0.5, ncol(counts))
    )
})


test_that("the grouped tilt does not depend on a rank cap", {

    counts <- make_tilt_counts(genes = 200L, cells = 4L)
    signature <- paste0("gene", 1:40)

    ranks <- as.matrix(MatrixGenerics::colRanks(
        -counts, ties.method = "average", preserveShape = TRUE,
        useNames = FALSE))
    dimnames(ranks) <- dimnames(counts)

    caps <- c(40, 80, 120, 160)
    capped <- vapply(caps, function(cap) {
        cap_vector <- stats::setNames(rep(cap, ncol(counts)),
                                      colnames(counts))
        mean(sparoscore(
            data = ranks,
            signatures = signature,
            data_has_ranks = TRUE,
            rank_caps = cap_vector
        ))
    }, numeric(1))

    tilt <- discrete_exponential_tie_censored(
        data = ranks,
        signatures = signature,
        data_has_ranks = TRUE
    )

    expect_gt(diff(range(capped)), 0)
    expect_equal(as.vector(tilt$score),
                 as.vector(discrete_exponential_tie_censored(counts, signature)$score))
})


test_that("a custom optimiser is called with the tie intervals", {

    counts <- make_tilt_counts(genes = 30L, cells = 2L)
    seen <- list()
    fake <- function(lo, hi, G, theta_bound) {
        seen[[length(seen) + 1L]] <<- list(lo = lo, hi = hi, G = G)
        c(score = 0.5, se = 1, theta = 0, se_theta = 1)
    }

    res <- discrete_exponential_tie_censored(
        counts,
        paste0("gene", 1:5),
        optimiser = fake
    )

    expect_equal(unname(res$theta[, 1]), c(0, 0))
    expect_equal(length(seen), ncol(counts))
    expect_true(all(vapply(seen, function(x) x$G, integer(1)) == 30L))
    expect_true(all(vapply(seen, function(x) {
        length(x$lo) == 5L && all(x$lo >= 1L & x$hi <= x$G & x$lo <= x$hi)
    }, logical(1))))
})


test_that("max_theoretical_theta is the continuous tilt of the top block", {

    # Top 200 of 1000 occupy [800, 1000], average 900, so m = 0.9 and
    # 1 - 1/theta = 0.9 gives theta = 10.
    expect_equal(max_theoretical_theta(1000, 200), 10)
    expect_equal(max_theoretical_theta(100, 2), 100)
    expect_equal(max_theoretical_theta(20000, 50), 800)
    expect_equal(max_theoretical_theta(10, 1), 20)
    expect_equal(max_theoretical_theta(10, 10), 0)
    expect_error(max_theoretical_theta(10, 11), "signature_size")
})


test_that("prob must lie between 0 and 1", {

    counts <- make_tilt_counts(genes = 20L, cells = 1L)
    expect_error(
        discrete_exponential_tie_censored(counts, paste0("gene", 1:3), prob = 95),
        "prob"
    )
})


test_that("a boundary posterior moves the median inside the support", {

    edge <- SPAROscore:::.posterior_interval(5, 1, 0.95,
                                             lower_bound = -5, upper_bound = 5)
    expect_lt(edge[["median"]], 5)
    expect_gte(edge[["lower"]], -5)
    expect_lte(edge[["upper"]], 5)
    expect_equal(edge[["lower"]],
                 stats::qnorm(stats::pnorm(-5, 5, 1) +
                              0.025 * (0.5 - stats::pnorm(-5, 5, 1)), 5, 1))
})


test_that("optimiser must be uniroot, stan, or a function", {

    counts <- make_tilt_counts(genes = 20L, cells = 1L)
    expect_error(
        discrete_exponential_tie_censored(counts, paste0("gene", 1:3), optimiser = "nlm"),
        "optimiser"
    )
})


test_that("stan optimiser matches an interior uniroot and stays finite on a flat likelihood", {

    skip_if_not_installed("instantiate")
    skip_if_not_installed("cmdstanr")
    skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan is not installed")

    counts <- make_tilt_counts(genes = 40L, cells = 2L)
    signature <- paste0("gene", 1:8)
    cache <- tempfile("sparoscore-stan-")

    by_root <- discrete_exponential_tie_censored(
        counts, signature, optimiser = "uniroot")
    by_stan <- discrete_exponential_tie_censored(
        counts, signature, optimiser = "stan",
        cache_stan_model = cache)

    expect_equal(by_stan$theta, by_root$theta, tolerance = 1e-3)
    expect_equal(by_stan$se_theta, by_root$se_theta, tolerance = 1e-2)
    expect_true(all(is.finite(by_stan$se_theta)))
    expect_true(all(by_stan$theta_lower < by_stan$theta_median))
    expect_true(all(by_stan$theta_median < by_stan$theta_upper))
    # Interior mode: truncation does not move the median off the mode, and
    # the 95% interval is the normal one with that standard error.
    expect_equal(by_stan$theta_median, by_stan$theta, tolerance = 1e-6)
    z <- stats::qnorm(0.975)
    expect_equal(as.vector(by_stan$theta_upper - by_stan$theta_lower),
                 as.vector(2 * z * by_stan$se_theta), tolerance = 1e-6)

    narrower <- discrete_exponential_tie_censored(
        counts, signature, optimiser = "stan", prob = 0.8,
        cache_stan_model = cache)
    expect_true(all(narrower$theta_upper - narrower$theta_lower <
                    by_stan$theta_upper - by_stan$theta_lower))

    # Every gene shares one count, so the likelihood does not depend on theta.
    # uniroot reports the edge of the search interval; the prior keeps Stan at 0.
    flat <- matrix(0, nrow = 10, ncol = 1,
                   dimnames = list(paste0("gene", 1:10), "cell"))
    root_flat <- discrete_exponential_tie_censored(
        flat, paste0("gene", 1:5), optimiser = "uniroot",
        theta_bound = 50)
    stan_flat <- discrete_exponential_tie_censored(
        flat, paste0("gene", 1:5), optimiser = "stan",
        prior_sd = 1000, cache_stan_model = cache)

    # Every gene shares the single block [1, G], so the likelihood is flat
    # and symmetry puts theta at 0 rather than on the edge of the bracket.
    expect_equal(unname(root_flat$theta[1, 1]), 0)
    expect_equal(unname(stan_flat$theta[1, 1]), 0, tolerance = 1e-6)
    # Flat likelihood: the posterior is the prior, so se_theta is prior_sd.
    expect_equal(unname(stan_flat$se_theta[1, 1]), 1000, tolerance = 1e-6)
    expect_equal(unname(stan_flat$theta_median[1, 1]), 0, tolerance = 1e-6)
    expect_equal(unname(stan_flat$theta_lower[1, 1]),
                 stats::qnorm(0.025) * 1000, tolerance = 1e-4)
    expect_equal(unname(stan_flat$theta_upper[1, 1]),
                 stats::qnorm(0.975) * 1000, tolerance = 1e-4)
})
