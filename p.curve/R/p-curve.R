#' Draw samples
#'
#' Draws samples of size `n` from two normal distributions.
#' @param delta Real difference in means
#' @param n Sample size for each group
#' @param sigma Standard deviation, same for both groups
#' @export

draw_samples = function(delta, n, sigma = 1) {
    group1 = rnorm(n, mean = 0, sd = sigma)
    group2 = rnorm(n, mean = delta, sd = sigma)
    return(list(group1, group2))
}

#' Tidy version of a t-test
#'
#' Given two samples, do an unpaired, two-sided t-test and return the results in tidy format. This simulates doing 1 study.
#' @param samples List of length 2; the two samples
#' @return A one-row `tibble::tibble`
#' @export
t_test = function(samples) {
    test = t.test(samples[[1]], samples[[2]])
    return(broom::tidy(test))
}

#' Do several studies
#'
#' Do several studies using draw_samples() and t_test().  This simulates doing a meta-analysis of `N` studies.
#' @param N Number of studies
#' @param delta Real effect size
#' @param n Sample size (per group)
#' @param seed Seed for the RNG; if `NULL` (the default) no seed will be used
#' @return A `tibble::tibble` with one row per study. Column `samples` contains the samples for each study.
#' @export
draw_studies = function(N, delta, n, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    tibble::tibble(study_id = 1:N) %>%
        dplyr::mutate(delta = delta,
                      n = n) %>%
        dplyr::mutate(samples = purrr::map2(delta, n,
                                            draw_samples)) %>%
        dplyr::mutate(test = purrr::map(samples, t_test)) %>%
        tidyr::unnest(test) %>%
        dplyr::mutate(rank = rank(p.value),
                      p.uniform = rank / (max(rank) + 1))
}



## Young's p-value plot ----

#' Young's p-value plot
#'
#' S. Stanley Young's "p-value plot" or "p-curve." Equivalent to a QQ-plot against a uniform distribution with a rescaled x-axis.
#' @param studies A dataframe of studies to be plotted, as returned by `draw_studies()`. Minimally, a dataframe with columns `rank` and `p.value`.
#' @param alpha Threshold for statistical significance; default (.05)
#' @param draw_alpha Draw a horizontal line at `alpha`?
#' @param ... Other aesthetics, passed at the plot level
#' @return A `ggplot2` plot
#' @details Young's p-curve plots the p-value of each study against its rank, by ascending p-value (1 = smallest p-value). It is equivalent to a QQ-plot against a uniform distribution, with the x-axis running from 1 to the number of studies rather than the quantiles `(1:N)/N`.
#' @export
young_curve = function(studies,
                       alpha = .05,
                       draw_alpha = TRUE,
                       ...) {
    plot = ggplot2::ggplot(studies,
                           ggplot2::aes(rank, p.value, ...)) +
        ggplot2::geom_point()
    if (draw_alpha) {
        plot = plot + ggplot2::geom_hline(yintercept = alpha)
    }
    return(plot)
}

#' Slope of Young's p-curve
#'
#' Calculates the slope of Young's p-curve using a linear regression.
#' @param studies A dataframe of studies, as returned by `draw_studies()`. Minimally, a dataframe with columns `rank` and `p.value`.
#' @return The slope, a length-1 numeric
#' @export
young_slope = function(studies) {
    model = lm(p.value ~ rank, data = studies)
    slope = model$coefficients[['rank']]
    return(slope)
}


#' QQ-plot
#'
#' A QQ-plot of the p-values against a uniform distribution.  Equivalent to Young's p-curve with a rescaled x-axis.
#' @param studies A dataframe of studies to be plotted, as returned by `draw_studies()`. Minimally, a dataframe with columns `p.uniform` and `p.value`.
#' @return A `ggplot2` plot
#' @details `p.uniform` is calculated in `draw_studies()` as `rank/(max(rank)+1)`. This plot includes a "45-degree line," a line with slope 1 and intercept 0.
#' @export
qq_plot = function(studies) {
    ggplot2::ggplot(studies,
                    ggplot2::aes(p.uniform, p.value)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 1, intercept = 0)
}


#' Slope of a QQ-plot
#'
#' Calculates the slope of the QQ-plot against the uniform distribution.
#' @param studies A dataframe of studies, as returned by `draw_studies()`. Minimally, a dataframe with columns `p.uniform` and `p.value`.
#' @return The slope, a length-1 numeric
#' @export
qq_slope = function(studies) {
    model = lm(p.value ~ p.uniform, data = studies)
    slope = model$coefficients[['p.uniform']]
    return(slope)
}
# qq_slope(studies)

#' Test linearity of a QQ-plot
#'
#' Checks the linearity of a QQ-plot (against a uniform distribution) by comparing linear and quadratic regressions of `p.value` against `p.uniform`.
#' @param studies A dataframe of studies, as returned by `draw_studies()`. Minimally, a dataframe with columns `p.uniform` and `p.value`.
#' @param alpha Threshold for statistical significance; used by the automated F-test.
#' @details Two checks are used: an F-test for the ANOVA of the two regression models, and a comparison of AICs.
#' @return A one-row `tibble::tibble` with the following columns:
#'   \item{f_stat}{Observed value of the F statistic}
#'   \item{alpha}{Threshold for statistical significance used in the F-test}
#'   \item{f_p}{p-value of the F statistic}
#'   \item{f_comp}{Inference from the F-test; one of `significant` or `non-significant`}
#'   \item{aic_linear}{AIC statistic of the linear regression}
#'   \item{aic_quad}{AIC statistic of the quadratic regression}
#'   \item{aic_comp}{Inference from comparing AICs; one of `quadratic` or `linear`}
qq_linear = function(studies, alpha = .05) {
    ## Checks linearity of the QQ plot by comparing linear and quadratic fits, using (a) F-test and (b) AICs
    model_linear = lm(p.value ~ p.uniform, data = studies)
    model_quad = lm(p.value ~ p.uniform + I(p.uniform^2),
                    data = studies)

    f_stat = anova(model_linear, model_quad)[['F']][[2]]
    f_p = anova(model_linear, model_quad)[['Pr(>F)']][[2]]
    f_comp = dplyr::if_else(f_p < alpha,
                     'significant',
                     'non-significant')
    aic_linear = AIC(model_linear)
    aic_quad = AIC(model_quad)
    aic_comp = dplyr::if_else(aic_quad < aic_linear,
                       'quadratic',
                       'linear')
    return(tibble::tibble(f_stat, alpha, f_p, f_comp,
                  aic_linear, aic_quad,
                  aic_comp))
}


## Simonsohn et al.'s p-curve ----
#' Simonsohn et al.'s p-curve
#'
#' Simonsohn et al.'s p-curve is equivalent to a histogram with binwidth .01 on the subset of p-values below .05
#' @param studies A dataframe of studies to be plotted, as returned by `draw_studies()`.  Minimally, a dataframe with column `p.value`.
#' @param ... Other aesthetics, passed at the plot level
#' @return A `ggplot2` plot
#' @details Very very small p-values (below 1e-90) are trimmed to keep the plot from drawing a value at 0.
#' @export
simonsohn_curve = function(studies, ...) {
    studies %>%
        filter(p.value < .05) %>%
        ggplot(aes(p.value, ...)) +
        geom_freqpoly(breaks = c(0:5 / 100),
                      position = position_nudge(x = .005)) +
        stat_bin(geom = 'point',
                 breaks = c(0:5 / 100),
                 position = position_nudge(x = .005)) +
        geom_rug() +
        xlim(1e-90, .05)
}

## Schweder and Spjøtvoll's p-value plot ----
#' Schweder and Spjøtvoll's p-value plot
#'
#' Schweder and Spjøtvolls "p-value plot."
#' @param studies A dataframe of studies to be plotted, as returned by `draw_studies()`. Minimally, a dataframe with columns `rank` and `p.value`.
#' @param ... Other aesthetics, passed at the plot level
#' @return A `ggplot2` plot
#' @details Schweder and Spjøtvoll's p-value plot plots the descending rank of each study (1 = largest p-value) against 1-p.  It stands in a 1-1 relationship with Young's p-value plot, and so a QQ-plot; but swaps the axes.
#' @export
schsp_curve = function(studies, ...) {
    ggplot(studies, aes(1 - p.value,
                        max(rank) - rank,
                        ...)) +
        geom_point()
}

#' Slope of Schweder and Spjøtvoll's p-value plot
#'
#' Calculates the slope of Schweder and Spjøtvoll's p-value plot using a linear regression.
#' @param studies A dataframe of studies, as returned by `draw_studies()`. Minimally, a dataframe with columns `rank` and `p.value`.
#' @return The slope, a length-1 numeric
#' @export
schsp_slope = function(studies) {
    model = studies %>%
        mutate(rank.rev = max(rank) - rank,
               p.value.rev = 1 - p.value) %>%
        lm(rank.rev ~ p.value.rev, data = .)
    slope = model$coefficients[['p.value.rev']]
    return(slope)
}


## Do many meta-analyses ----
#' Run the simulation
#'
#' Top-level function for running the simulation
#' @param NN Number of times to run the simulation for each condition
#' @param N Number of studies to conduct in each simulation run (vector)
#' @param delta Real effect size (vector)
#' @param n Sample size per group for each simulated study (vector)
#' @details The condition parameters `N`, `delta`, and `n` can be vectors; the function will automatically construct all combinations.  `NN` must have length 1.
#' @return A dataframe containing one row per simulation run and the following columns:
#'   \item{meta_idx}{Index of the simulation run within the condition}
#'   \item{N, delta, n}{Simulation conditions}
#'   \item{studies}{Nested samples and t-test results for each study}
#'   \item{young_slope, schsp_slope, qq_slope}{Slopes for the Young, Schweder and Spjøtvoll, and QQ plots}
#'   \item{f_stat, alpha, f_p, f_comp}{F statistic, alpha threshold, p-value, and inference for the F test of linearity of the QQ plot}
#'   \item{aic_linear, aic_quad, aic_comp}{AIC values for the linear and quadratic regression, and inference, for the AIC test of linearity of the QQ plot}
many_metas = function(NN,
                      N, delta, n) {
    assertthat::assert_that(assertthat::are_equal(length(NN),
                                                  1L),
                            msg = 'NN must be length 1')

    meta_params = purrr::cross_df(tibble::lst(meta_idx = 1:NN,
                                              N, delta, n))
    # return(meta_params)

    meta_params %>%
        dplyr::mutate(studies = purrr::pmap(list(N, delta, n),
                                            draw_studies)) %>%
        dplyr::mutate(young_slope = purrr::map_dbl(studies,
                                                   young_slope),
               schsp_slope = purrr::map_dbl(studies, schsp_slope),
               qq_slope = purrr::map_dbl(studies, qq_slope),
               qq_linear = purrr::map(studies, qq_linear)) %>%
        tidyr::unnest(qq_linear)
}


