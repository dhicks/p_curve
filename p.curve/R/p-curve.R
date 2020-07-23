library(magrittr)

#' Draw samples
#'
#' Draws samples of size `n` from two normal distributions.
#' @param delta Real difference in means
#' @param n Sample size for each group
#' @param sigma Standard deviation, same for both groups
#' @details If `delta` has length > 1, means will be recycled, silently creating subsamples within the second group.  This simulates a case where the intervention group is heterogeneous, with even samples from two latent populations.
#' @export

draw_samples = function(delta, n, sigma = 1) {
    group1 = rnorm(n, mean = 0, sd = sigma)
    group2 = rnorm(n, mean = delta, sd = sigma)
    return(list(group1, group2))
}

#' Tidy version of a t-test
#'
#' Given two samples, do an unpaired, two-sided t-test and return the results in tidy format. This simulates the analysis for a single study.
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
#' @param delta Real effect size, or multiple sizes for a mixed/heterogeneous case
#' @param n Sample size (per group)
#' @param seed Seed for the RNG; if `NULL` (the default) no seed will be used
#' @return A `tibble::tibble` with one row per study. Column `samples` contains the samples for each study.
#' @export
draw_studies = function(N, delta, n, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    tibble::tibble(study_id = 1:N) %>%
        dplyr::mutate(delta = list(delta),
                      n = n) %>%
        dplyr::mutate(samples = purrr::map2(delta, n,
                                            ~draw_samples(purrr::simplify(.x),
                                                          .y))) %>%
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
#' Checks the linearity of a QQ-plot (against a uniform distribution) by comparing linear and quadratic regressions of `p.value` against `p.uniform` and using a Kolmogorov-Smirnov test of uniformity.
#' @param studies A dataframe of studies, as returned by `draw_studies()`. Minimally, a dataframe with columns `p.uniform` and `p.value`.
#' @param alpha Threshold for statistical significance; used by F-test and Kolmogorov-Smirnof test
#' @details Two checks are used: an F-test for the ANOVA of the two regression models, and a comparison of AICs.
#' @return A one-row `tibble::tibble` with the following columns:
#'   \item{f_stat}{Observed value of the F statistic}
#'   \item{alpha}{Threshold for statistical significance used in the F-test and KS test}
#'   \item{f_p}{p-value of the F statistic}
#'   \item{f_comp}{Inference from the F-test; one of `non-linear` or `linear`}
#'   \item{aic_linear}{AIC statistic of the linear regression}
#'   \item{aic_quad}{AIC statistic of the quadratic regression}
#'   \item{aic_comp}{Inference from comparing AICs; one of `non-linear` or `linear`}
#'   \item{ks_stat}{Observed value of the KS statistic}
#'   \item{ks_p}{p-value of the KS statistic}
#'   \item{ks_comp}{Inference from the KS stat; one of `non-linear` or `linear`}
qq_linear = function(studies, alpha = .05) {
    model_linear = lm(p.value ~ p.uniform, data = studies)
    model_quad = lm(p.value ~ p.uniform + I(p.uniform^2),
                    data = studies)

    f_stat = anova(model_linear, model_quad)[['F']][[2]]
    f_p = anova(model_linear, model_quad)[['Pr(>F)']][[2]]
    f_comp = dplyr::if_else(f_p < alpha,
                            'non-linear',
                            'linear')
    aic_linear = AIC(model_linear)
    aic_quad = AIC(model_quad)
    aic_comp = dplyr::if_else(aic_quad < aic_linear,
                              'non-linear',
                              'linear')

    ks_test = ks.test(studies$p.value, studies$p.uniform)
    ks_stat = ks_test$statistic
    ks_p = ks_test$p.value
    ks_comp = dplyr::if_else(ks_p < alpha, 'non-linear', 'linear')


    return(tibble::tibble(f_stat, alpha, f_p, f_comp,
                          aic_linear, aic_quad, aic_comp,
                          ks_stat, ks_p, ks_comp))
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
        dplyr::mutate(rank.rev = max(rank) - rank,
                      p.value.rev = 1 - p.value) %>%
        lm(rank.rev ~ p.value.rev, data = .)
    slope = model$coefficients[['p.value.rev']]
    return(slope)
}

#' Convenience function to coerce single and mixed deltas to characters
flatten_to_chr = function(alist) {
    len = purrr::map_int(alist, length)
    return = purrr::map2_chr(alist, len,
                             ~ifelse(identical(.y, 1L),
                                     as.character(.x),
                                     'mixed'))
    return(return)
}
# flatten_to_chr(list(0, .2, list(0, .7)))

## Do many meta-analyses ----
#' Run the simulation: "Do many meta-analyses"
#'
#' Top-level function for running the simulation
#' @param NN Number of times to run the simulation for each condition
#' @param N Number of studies to conduct in each simulation run (vector)
#' @param delta Real effect size (vector)
#' @param n Sample size per group for each simulated study (vector)
#' @details The condition parameters `N`, `delta`, and `n` can be vectors; the function will automatically construct all combinations.  `NN` must have length 1.
#' @return A dataframe containing one row per simulation run and the following columns:
#'   \item{meta_idx}{Index of the simulation run within the condition}
#'   \item{N, delta, n}{Simulation conditions.  To pass mixed conditions, construct `delta` as a list of lists.}
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
                                            draw_studies),
                      delta_chr = flatten_to_chr(delta)) %>%
        dplyr::mutate(young_slope = purrr::map_dbl(studies,
                                                   young_slope),
                      schsp_slope = purrr::map_dbl(studies, schsp_slope),
                      qq_slope = purrr::map_dbl(studies, qq_slope),
                      qq_linear = purrr::map(studies, qq_linear)) %>%
        tidyr::unnest(qq_linear)
}


# Assess evidence ----
#' Calculate p-value
#'
#' Calculate the p-value, as the share of simulation runs in which `test_output` is true, given the condition `h`
#' @param h Bare (unquoted) expression for H0 (eg, `delta == 0.2`)
#' @param test_output Bare (unquoted) expression for the test output (eg, `aic_comp == 'quadratic'`)
#' @param dataf Data frame of simulation results, as returned by `many_metas()`
#' @return Dataframe with columns
#'   \item{h}{Expression for H0, as a string}
#'   \item{test_output}{Expression for the test output, as a string}
#'   \item{n_false}{Number of rows post-filtering where `test_output` is true}
#'   \item{n_true}{Number of rows post-filtering where `test_output` is false}
#'   \item{p}{Share of rows `n_true/(n_false + n_true)`}
#' @export
p_value = function(h, test_output, dataf) {
    h = rlang::enquo(h)
    test_output = rlang::enquo(test_output)

    # return(rlang::as_label(h))

    counted_df = dataf %>%
        dplyr::filter(!!h) %>%
        dplyr::mutate(test_output := !!test_output) %>%
        dplyr::count(test_output, wt = n) %>%
        dplyr::mutate(share = n / sum(n))

    if (identical(nrow(counted_df), 1L)) {
        ## If there's only 1 row, all the test_output values are the same
        if (any(counted_df$test_output)) {
            ## So if any test_output values are TRUE,
            ## then they're all TRUE
            p = counted_df %>%
            transmute(n_false = 0,
                      n_true = n,
                      share_true = share)
        } else {
            ## Otherwise they're all FALSE
            p = counted_df %>%
                transmute(n_false = n,
                          n_true = 0,
                          share_true = 1 - share)
        }
    } else {
        ## If there's more than one row, they have different values
        p = counted_df %>%
            tidyr::pivot_wider(names_from = test_output,
                               values_from = c(n, share),
                               names_repair = tolower)
    }

    return_df = p %>%
        dplyr::mutate(h = rlang::as_label(h),
                      test_output = rlang::as_label(test_output)) %>%
        dplyr::select(h, test_output, n_false, n_true, p = share_true)

    return(return_df)
}

#' Calculate the log likelihood ratio
#'
#' Calculate the log (base 10) likelihood ratio for the sets of hypotheses `h1` and `h2`.  Note that the interface is slightly different from `p_value()`: for the hypothesis and test output, that function takes bare expressions.  This function takes a list of calls, as created by `exprs()`.
#' @param h1 List of `call` (created by `exprs()`) for H1
#' @param h2 List of `call` for H2
#' @param test_output List of `call` for the test output
#' @param dataf Data frame of simulation results, as returned by `many_metas()`
#' @return Dataframe with columns
#'   \item{h1, h2, test_output}{Expression for H1, H2, and test output, as strings}
#'   \item{llr} Log (base 10) likelihood ratio L(H1; d) / L(H2; d)
#'   \item{n_false_h1, n_false_h2}{Number of rows post-filtering where `test_output` is true, for H1 and H2}
#'   \item{n_true_h1, n_true_h2}{Number of rows post-filtering where `test_output` is false, for H1 and H2}
likelihood_ratio = function(h1, h2, test_output, dataf) {
    h1_df = purrr::cross(tibble::lst(h1, test_output)) %>%
        purrr::map_dfr(~p_value(!!.[['h1']], !!.[['test_output']], dataf))
    h2_df = purrr::cross(tibble::lst(h2, test_output)) %>%
        purrr::map_dfr(~p_value(!!.[['h2']], !!.[['test_output']], combined_df))
    # return(h2_df)

    purrr::cross_df(lst(h1, h2)) %>%
        dplyr::mutate(across(c(h1, h2), ~purrr::map_chr(., rlang::as_label))) %>%
        dplyr::filter(h1 != h2) %>%
        dplyr::left_join(h1_df, by = c('h1' = 'h')) %>%
        dplyr::left_join(h2_df, by = c('h2' = 'h', 'test_output'),
                         suffix = c('_h1', '_h2')) %>%
        dplyr::rename(l_h1 = p_h1, l_h2 = p_h2) %>%
        dplyr::mutate(llr = log10(l_h1) - log10(l_h2)) %>%
        dplyr::select(h1, h2, test_output, llr, l_h1, l_h2,
                      dplyr::everything())
}


