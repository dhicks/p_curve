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



## Plots ----

#' Young's p-curve
#'
#' S. Stanley Young's "p-value plot" or "p-curve."
#' @param studies A dataframe of studies to be plotted, as returned by `draw_studies()`. Minimally, a dataframe with columns `rank` and `p.value`.
#' @param alpha Threshold for statistical significance; default (.05)
#' @param draw_alpha Draw a horizontal line at `alpha`?
#' @return A `ggplot2` plot
#' @details Young's p-curve plots the p-value of each study against its rank, by ascending p-value (1 = smallest p-value). It is equivalent to a QQ-plot against a uniform distribution, with the x-axis running from 1 to the number of studies rather than the quantiles `(1:N)/N`.
#' @export
young_curve = function(studies,
                       alpha = .05,
                       draw_alpha = TRUE) {
    plot = ggplot2::ggplot(studies,
                           ggplot2::aes(rank, p.value)) +
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
#' @param f_stat Observed value of the F statistic
#' @param alpha Threshold for statistical significance used in the F-test
#' @param f_p p-value of the F statistic
#' @param f_comp Inference from the F-test; one of `significant` or `non-significant`
#' @param aic_linear AIC statistic of the linear regression
#' @param aic_quad AIC statistic of the quadratic regression
#' @param aic_comp Inference from comparing AICs; one of `quadratic` or `linear`
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
## TODO
# studies_med$studies[[2]] %>%
#     # filter(p.value < .05) %>%
#     ggplot(aes(p.value)) +
#     stat_bin(binwidth = .01, geom = 'point') +
#     stat_bin(binwidth = .01, geom = 'line')

#*[pick up documentation here]*
## Schweder and Spjøtvoll's p-curve ----
schsp_curve = function(studies) {
    ggplot(studies, aes(1 - p.value,
                        max(rank) - rank)) +
        geom_point()
}

## Schweder and Spjøtvoll slope
schsp_slope = function(studies) {
    model = studies %>%
        mutate(rank.rev = max(rank) - rank,
               p.value.rev = 1 - p.value) %>%
        lm(rank.rev ~ p.value.rev, data = .)
    slope = model$coefficients[['p.value.rev']]
    return(slope)
}


## Do many meta-analyses ----
many_metas = function(NN, ## how many meta-studies
                      N, delta, n) {
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


