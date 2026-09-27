# helper testing function
make_mle_counts <- function(genes = 200L, cells = 4L, seed = 11L) {
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


test_that("sparoscore_mle matrix dispatch returns correct structure", {

    counts <- make_mle_counts()

    res <- sparoscore_mle(
        data = counts,
        signatures = paste0("gene", 1:20)
    )

    expect_true(is.matrix(res))
    expect_equal(nrow(res), ncol(counts))
    expect_equal(colnames(res), "SPAROscore_MLE")

    # the uncertainty travels with the score
    errors <- attr(res, "standard_error")
    expect_equal(dim(errors), dim(res))
    expect_true(all(errors > 0))
})


test_that("sparoscore_mle scores lie on the unit interval", {

    counts <- make_mle_counts()

    res <- sparoscore_mle(
        data = counts,
        signatures = paste0("gene", 1:20)
    )

    expect_true(all(res > 0 & res < 1))
})


test_that("sparoscore_mle return_se gives both matrices", {

    counts <- make_mle_counts()

    fitted <- sparoscore_mle(
        data = counts,
        signatures = paste0("gene", 1:20),
        return_se = TRUE
    )

    expect_named(fitted, c("score", "se"))
    expect_equal(dim(fitted$score), dim(fitted$se))
    expect_equal(
        as.vector(fitted$score),
        as.vector(sparoscore_mle(counts, paste0("gene", 1:20)))
    )
})


test_that("sparoscore_mle agrees across matrix-like input classes", {

    skip_if_not_installed("Matrix")

    counts <- make_mle_counts()
    signature <- paste0("gene", 1:20)

    dense <- sparoscore_mle(counts, signature)
    sparse <- sparoscore_mle(Matrix::Matrix(counts, sparse = TRUE), signature)
    frame <- sparoscore_mle(as.data.frame(counts), signature)

    expect_equal(as.vector(dense), as.vector(sparse))
    expect_equal(as.vector(dense), as.vector(frame))
})


test_that("sparoscore_mle accepts pre-computed ranks", {

    counts <- make_mle_counts()
    signature <- paste0("gene", 1:20)

    ranked <- get_ranks(counts)

    from_counts <- sparoscore_mle(counts, signature)
    from_ranks <- sparoscore_mle(
        data = ranked$ranks,
        signatures = signature,
        data_has_ranks = TRUE,
        rank_caps = ranked$rank_caps
    )

    expect_equal(as.vector(from_counts), as.vector(from_ranks))
    expect_error(
        sparoscore_mle(ranked$ranks, signature, data_has_ranks = TRUE),
        "rank_caps needed"
    )
})


test_that("sparoscore_mle handles multiple signatures and a prefix", {

    counts <- make_mle_counts()

    res <- sparoscore_mle(
        data = counts,
        signatures = list(
            SignatureA = paste0("gene", 1:20),
            SignatureB = paste0("gene", 50:80)
        ),
        prefix = "mle_"
    )

    expect_equal(colnames(res), c("mle_SignatureA", "mle_SignatureB"))
    expect_equal(colnames(attr(res, "standard_error")), colnames(res))
    expect_error(sparoscore_mle(counts, list(paste0("gene", 1:20))), "named")
})


test_that("sparoscore_mle subtracts down_signatures and pools their error", {

    counts <- make_mle_counts()
    up <- list(Signature = paste0("gene", 1:20))
    down <- list(Signature = paste0("gene", 50:80))

    combined <- sparoscore_mle(counts, up, down_signatures = down,
                               return_se = TRUE)
    up_only <- sparoscore_mle(counts, up, return_se = TRUE)
    down_only <- sparoscore_mle(counts, down, return_se = TRUE)

    expect_equal(combined$score, up_only$score - down_only$score)
    expect_equal(combined$se, sqrt(up_only$se^2 + down_only$se^2))
    expect_error(
        sparoscore_mle(counts, up, down_signatures = list(Other = "gene1")),
        "should match"
    )
})


test_that("tilt means stay finite at the theta = 0 singularity", {

    # the expanded score has 1/theta poles that cancel; the series is what
    # remains, and the closed form is already garbage at 1e-8
    expect_equal(SPAROscore:::.tilt_mean(0), 0.5)
    expect_equal(SPAROscore:::.tilt_cond_mean(0.4, 0), 0.2)
    expect_equal(
        SPAROscore:::.tilt_mean(1e-8),
        0.5 + 1e-8 / 12,
        tolerance = 1e-16
    )
    expect_equal(
        SPAROscore:::.tilt_cond_mean(0.4, 1e-8),
        0.2 + 1e-8 * 0.16 / 12,
        tolerance = 1e-16
    )
})


test_that("uniroot can sit on the uniform without returning NaN", {

    # R + Q * c/2 - G/2 = 2.1 + 2 * 0.2 - 2.5 = 0, so the root is theta = 0
    fitted <- SPAROscore:::.fit_tilt_per_cell(
        positions = c(0.5, 0.7, 0.9),
        censored_size = 2,
        cap_quantile = 0.4
    )

    expect_equal(unname(fitted[["theta"]]), 0, tolerance = 1e-8)
    expect_equal(unname(fitted[["score"]]), 0.5, tolerance = 1e-8)
    expect_true(is.finite(fitted[["se"]]))
})


test_that("sparoscore_mle reports NA when nothing is detected", {

    counts <- make_mle_counts()
    signature <- paste0("gene", 1:20)

    # push every signature gene to the bottom of the first cell and everything
    # else to the top, so no signature gene survives the cap
    counts[signature, 1] <- 0L
    counts[setdiff(rownames(counts), signature), 1] <- 100L

    fitted <- sparoscore_mle(counts, signature, return_se = TRUE)

    expect_true(is.na(fitted$score[1, 1]))
    expect_true(is.infinite(fitted$se[1, 1]))

    # the remaining cells are untouched and still scored
    expect_true(all(is.finite(fitted$score[-1, 1])))
})


test_that("sparoscore_mle imputes missing genes as censored", {

    counts <- make_mle_counts()
    present <- paste0("gene", 1:20)
    padded <- c(present, "absent1", "absent2", "absent3")

    skipped <- suppressWarnings(
        sparoscore_mle(counts, padded, handle_missing_genes = "skip"))
    imputed <- suppressWarnings(
        sparoscore_mle(counts, padded, handle_missing_genes = "impute"))

    expect_equal(as.vector(skipped),
                 as.vector(sparoscore_mle(counts, present)))
    expect_true(all(imputed < skipped))
    expect_error(
        suppressWarnings(
            sparoscore_mle(counts, padded, handle_missing_genes = "nonsense")),
        "handle_missing_genes"
    )
})


test_that("sparoscore_mle is far less cap dependent than the footrule score", {

    counts <- make_mle_counts(genes = 2000L, cells = 4L)
    signature <- paste0("gene", 1:200)

    ranked <- get_ranks(counts)
    caps <- c(400, 800, 1200, 1600)

    spread <- function(scorer) {
        means <- vapply(caps, function(cap) {
            cap_vector <- stats::setNames(rep(cap, ncol(counts)),
                                          colnames(counts))
            mean(scorer(
                data = ranked$ranks,
                signatures = signature,
                data_has_ranks = TRUE,
                rank_caps = cap_vector
            ))
        }, numeric(1))
        diff(range(means))
    }

    expect_lt(spread(sparoscore_mle), spread(sparoscore) / 5)
})
