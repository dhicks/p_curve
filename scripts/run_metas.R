#' ---
#' title: "Analysis for 'Young\\'s p-value plot does not provide evidence against air pollution hazards'"
#' author: "Daniel J. Hicks, University of California, Merced"
#' email: <hicks.daniel.j@gmail.com>
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#' ---
#' 
#' This document conducts the simulation described in the paper "Young's p-value plot does not provide evidence against air pollution hazards" by Daniel J. Hicks, and generates all plots and tables presented in that paper.  
#' 
#' The source code for this project is available at <https://github.com/dhicks/p_curve>.  The most recent Github version is automatically reproduced using Travis-CI, and the resulting `html` version of this file can be viewed at <https://dhicks.github.io/p_curve/>. 
#' 
#' After cloning the repository and installing the R dependencies listed below, the analysis, output files, and `html` version of this file can be recreated from the command line by calling `Rscript -e "rmarkdown::render('run_metas.R')"` within the `scripts` directory.  Alternavely, if `make` is installed, call `make script` from the command line within the top folder of the repository, or `make knit` from within the `scripts` folder.   
#' 
#' Note that reproducing the analysis creates a top-level `out` folder for plots and tables if this folder does not already exist, and quietly overwrites existing versions of the output files in this folder.  
#' 
#' The workhorse functions to conduct and analyze the simulation are all found in the `p.curve` folder, which is loaded as a local package using `devtools`.  Documentation for all of these functions can be constructed by calling `devtools::document(file.path('..', 'p.curve'))` from within the `scripts` folder, and then queried in the usual way, e.g., `?many_metas`. 
#' 

#' # Setup
#+ setup
library(devtools)
library(rlang)
library(tidyverse)
theme_set(theme_bw())
library(cowplot)
library(plotly)
library(knitr)

library(tictoc)

## Local package with simulation functions
load_all(file.path('..', 'p.curve'))
## Run the next line to recreate the documentation for this package
# document(file.path('..', 'p.curve'))

out_folder = file.path('..', 'out')
if (!dir.exists(out_folder)) {
    dir.create(out_folder)
}

write_plot = function(name, ...) {
    ggsave(file.path(out_folder, str_c(name, '.png')), ...)
}
plot_defaults = list(height = 4, width = 6, scale = 1.5)

#' # Run simulations
#+ run simulations
## Run simulations ----
## n = 60 gives us ~60% power to detect a moderate effect
power.t.test(delta = c(.2, .4, .6), sd = 1, sig.level = 0.05, n = 60)

## 5 effect sizes x 500 meta-studies x 20 studies x 25 samples
## 2-3 min
{
    set.seed(2020-10-17)
    tic()
    NN = 500  ## How many simulations to run? 
    combined_df = many_metas(NN = NN, 
                             delta = list(0, .2, .4, .6, list(0, .6)), 
                             N = 20, n = 60) %>% 
        mutate(delta_fct = as_factor(delta_chr))
    toc()
}

combined_df



#' # Model validation
#' Here we validate the model, showing that, at both the study level and meta-analysis level, the simulation gives an unbiased (ie, in expectation) estimate of the true effect, except for the mixed case
#' 
#+ model validation
## Model validation ----
## Distribution of effects estimates
combined_df %>%
    select(delta_fct, meta_idx, studies) %>%
    unnest(studies) %>%
    ggplot(aes(delta_fct, estimate2, color = delta_fct)) +
    geom_violin(draw_quantiles = .5) +
    labs(x = 'real effect',
         y = 'study-level effect estimate') +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

do.call(write_plot, c('estimates_study', plot_defaults))

combined_df %>%
    select(delta_fct, meta_idx, studies) %>%
    unnest(studies) %>%
    group_by(delta_fct) %>%
    summarize(across(estimate2, lst(mean, sd)), `NN*N` = n())

## Simulation-run/meta-analytic effects estimate
## This shows that, at the meta-analysis level (run of the simulation across N studies), the simulation gives an unbiased (in expectation) estimate of the effect, with less variance than the individual study estimates, except for the mixed case
combined_df %>%
    select(delta_fct, meta_idx, studies) %>%
    unnest(studies) %>%
    group_by(delta_fct, meta_idx) %>%
    summarize(mean = mean(estimate2)) %>%
    ungroup() %>%
    ggplot(aes(delta_fct, mean, color = delta_fct)) +
    geom_violin(draw_quantiles = .5) +
    labs(x = 'real effect',
         y = 'meta-analytic effect estimate') +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

do.call(write_plot, c('estimates_meta', plot_defaults))

## Because the mean of the mean is just the mean, this gives the same values as the table above
# combined_df %>%
#     select(delta_fct, meta_idx, studies) %>%
#     unnest(studies) %>%
#     group_by(delta_fct, meta_idx) %>%
#     summarize(across(estimate2, lst(mean, sd))) %>%
#     group_by(delta_fct) %>%
#     summarize(across(estimate2_mean, lst(mean)), NN = n())


#' # Sample plots
#+ sample plots
## Sample plots ----
# sample_slice = seq.int(from = 13, by = 3, length.out = 10)
set.seed(2020-07-23)
sample_slice = sample(1:NN, 5)
samples_par = list(height = 3.75, width = 6, scale = 2)

## Young
combined_df %>%
    group_by(delta) %>%
    slice(sample_slice) %>%
    ungroup() %>%
    select(delta_fct, meta_idx, studies) %>%
    unnest(studies) %>%
    young_curve(color = delta_fct) +
    facet_grid(rows = vars(delta_fct), cols = vars(meta_idx),
               as.table = TRUE,
               switch = 'y') +
    labs(x = 'rank (ascending)',
         y = 'p') +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

do.call(write_plot, c('fig_1_samples_young', samples_par))

young_composite(combined_df, alpha = .05, 
                color = delta_fct) + 
    facet_wrap(vars(delta_fct), scales = 'free_x') +
    labs(x = 'rank (ascending)',
         y = 'p') +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

do.call(write_plot, c('young_composite', samples_par))

## Sch-Sp
combined_df %>%
    group_by(delta) %>%
    slice(sample_slice) %>%
    ungroup() %>%
    select(delta_fct, meta_idx, studies) %>%
    unnest(studies) %>%
    schsp_curve(color = delta_fct) +
    facet_grid(rows = vars(delta_fct), cols = vars(meta_idx),
               as.table = TRUE,
               switch = 'y') +
    labs(x = '1 - p',
         y = 'rank (descending)') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

do.call(write_plot, c('samples_schsp', samples_par))

## Simonsohn et al.
combined_df %>%
    group_by(delta) %>%
    slice(sample_slice) %>%
    ungroup() %>%
    select(delta_fct, meta_idx, studies) %>%
    unnest(studies) %>%
    simonsohn_curve(color = delta_fct) +
    facet_grid(rows = vars(delta_fct), cols = vars(meta_idx),
               as.table = TRUE,
               switch = 'y') +
    labs(x = 'p-value',
         y = 'count') +
    scale_x_continuous(limits = c(1e-90, .05), 
                       breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

do.call(write_plot, c('samples_simonsohn', samples_par))


#' # Gaps
#+ gaps
## Gaps ----
combined_df %>% 
    select(delta_fct, meta_idx, gap, gappy) %>% 
    ggplot(aes(delta_fct, gap, color = delta_fct)) +
    geom_violin(draw_quantiles = .5, fill = NA) +
    geom_hline(yintercept = .125, alpha = .5) +
    labs(x = 'real effect', y = 'size of largest gap') +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

combined_df %>% 
    select(delta_fct, meta_idx, gap, gappy) %>% 
    ggplot(aes(delta_fct, fill = gappy)) +
    geom_bar(position = position_fill(reverse = TRUE)) +
    scale_fill_viridis_d(option = 'A') +
    scale_y_continuous(name = 'share', 
                       labels = scales::percent_format()) +
    labs(x = 'real effect', y = 'share')


#' # Slopes
#+ slopes
## Slopes ----
combined_df %>%
    rename(qq_slope = qq_estimate) %>% 
    select(delta_fct, meta_idx, matches('_slope$')) %>%
    rename('QQ' = qq_slope, 
           'Schweder and Spjøtvoll' = schsp_slope, 
           'Young' = young_slope) %>% 
    pivot_longer(cols = !c(delta_fct, meta_idx),
                 names_to = 'method', values_to = 'slope') %>%
    ggplot(aes(delta_fct, slope, color = delta_fct)) +
    # geom_beeswarm(alpha = .10) +
    geom_violin(draw_quantiles = .5, fill = NA) +
    labs(x = 'real effect') +
    scale_color_brewer(palette = 'Set1', guide = FALSE) +
    facet_wrap(vars(method), scales = 'free')

do.call(write_plot, c('slopes', plot_defaults))

combined_df %>%
    rename(qq_slope = qq_estimate) %>% 
    group_by(delta_fct) %>%
    summarize_at(vars(matches('_slope$')),
                 lst(median, sd)) %>% 
    pivot_longer(-delta_fct, 
                 names_to = c('plot', 'stat'), 
                 names_sep = '_slope_',
                 values_to = 'value') %>% 
    pivot_wider(names_from = stat, 
                values_from = value) %>% 
    mutate(plot = case_when(plot == 'young' ~ 'Young', 
                            plot == 'schsp' ~ 'Schweder and Spjøtvoll', 
                            plot == 'qq' ~ 'QQ',
                            TRUE ~ NA_character_)) %>% 
    select(plot, delta = delta_fct, everything()) %>% 
    arrange(plot, delta) %>% 
    kable(format = 'latex', booktabs = TRUE,
          digits = 2,
          caption = 'Slopes of QQ-plots, Schweder and Spjøtvoll\'s p-value plot, and Young\'s p-value plot, by condition (effect size $\\delta$).', 
          label = 'slopes') %>% 
    write_lines(file.path(out_folder, 'slopes.tex'))

## There's not an exact relationship between Young and Sch-Sp slopes,
## even when N (and so the max rank) is constant
ggplot(combined_df, aes(young_slope, schsp_slope)) +
    geom_point() +
    scale_y_log10()

## Whereas Young slope is a rescaling of QQ slope;
## the x-axis in the QQ plot is rank/max(rank), and rank is the x-axis in the Young plot
combined_df %>% 
    rename(qq_slope = qq_estimate) %>% 
    ggplot(aes(young_slope, qq_slope)) +
    geom_point()

## Distribution of calls for the slope
combined_df %>% 
    count(delta_fct, qq_t.comp, qq_tost.comp) %>% 
    pivot_longer(cols = c(qq_t.comp, qq_tost.comp), 
                 names_to = 'test', 
                 values_to = 'call') %>% 
    ggplot(aes(delta_fct, n, fill = fct_rev(call))) +
    geom_col(position = 'fill') +
    facet_wrap(vars(test)) +
    scale_fill_viridis_d(option = 'C', name = 'Inference') +
    xlab('real effect') +
    scale_y_continuous(labels = scales::percent_format(),
                       name = 'share of inferences')

## KS uniformity test
ggplot(combined_df, aes(delta_fct, fill = qq_ks.comp)) +
    geom_bar(position = 'fill') +
    scale_fill_viridis_d(option = 'C', name = 'KS inference') +
    xlab('real effect') +
    scale_y_continuous(labels = scales::percent_format(),
                       name = 'share of inferences')


#' # QQ linearity tests
#+ QQ linearity tests
## QQ linearity tests ----
# test_levels = c('AIC: linear' = 'linear',
#                 'AIC: quadratic' = 'quadratic',
#                 'F-test: non-sig.' = 'non-significant',
#                 'F-test: sig.' = 'significant')

combined_df %>%
    select(delta_fct, meta_idx, f_comp, aic_comp) %>%
    pivot_longer(c(f_comp, aic_comp),
                 names_to = 'test') %>%
    count(delta_fct, test, value) %>%
    mutate(test = fct_recode(test, 
                             'AIC' = 'aic_comp', 
                             'F-test' = 'f_comp')) %>%
    ggplot(aes(delta_fct, n, fill = fct_rev(value))) +
    geom_col(position = 'fill') +
    facet_wrap(vars(test), nrow = 3) +
    scale_fill_viridis_d(option = 'C', name = 'Inference') +
    xlab('real effect') +
    scale_y_continuous(labels = scales::percent_format(),
                       name = 'share of inferences')

do.call(write_plot, c('linearity', plot_defaults))

combined_df %>%
    select(delta_fct, meta_idx, f_comp, aic_comp) %>%
    pivot_longer(c(f_comp, aic_comp),
                 names_to = 'test') %>%
    count(delta_fct, test, value) %>%
    group_by(delta_fct, test) %>% 
    mutate(share = n/sum(n)) %>% 
    mutate(test = fct_recode(test, 
                             'AIC' = 'aic_comp', 
                             'F-test' = 'f_comp', 
                             'KS test' = 'ks_comp')) %>% 
    rename(`$\\delta$` = delta_fct) %>% 
    kable(format = 'latex', booktabs = TRUE, escape = FALSE,
          linesep = '',
          digits = 2, 
          caption = 'Distributions of outcomes for the linearity tests across conditions.  Columns ``linear" and ``non-linear" indicate counts of simulation runs in which tests indicate the plot was linear or non-linear, respectively.  ``Share" indicates the fraction of simulation runs in which the test indicated non-linearity.', 
          label = 'linearity') %>% 
    write_lines(file.path(out_folder, 'linearity.tex'))

#' # Combined accuracy table
#+ combined accuracy
## Combined accuracy table ----
combined_df %>%
    select(delta_fct, meta_idx, 
           gappy, 
           qq_t.comp, qq_tost.comp, qq_ks.comp,
           f_comp, aic_comp) %>%
    pivot_longer(c(gappy, 
                   qq_t.comp, qq_tost.comp, qq_ks.comp,
                   f_comp, aic_comp),
                 names_to = 'test') %>%
    mutate(exp_value = case_when(test %in% c('gappy') ~ 'not gappy',
                                 delta_fct == 0 & test %in% c('qq_t.comp', 'qq_tost.comp') ~ 'slope = 1',
                                 delta_fct == 0 & test %in% c('qq_ks.comp') ~ 'uniform',
                                 delta_fct == 0 & test %in% c('f_comp', 'aic_comp') ~ 'linear',
                                 delta_fct != 0 & test %in% c('qq_t.comp', 'qq_tost.comp') ~ 'slope ≠ 1',
                                 delta_fct != 0 & test %in% c('qq_ks.comp')~ 'non-uniform',
                                 delta_fct != 0 & test %in% c('f_comp', 'aic_comp') ~ 'non-linear',
                                 TRUE ~ NA_character_),
           correct_call = value == exp_value) %>%
    group_by(delta_fct, test) %>%
    summarize(accuracy = sum(correct_call) / n()) %>%
    ungroup() %>%
    pivot_wider(names_from = test, values_from = accuracy) %>%
    select(`real effect` = delta_fct, 
           Gap = gappy, 
           `T-test` = qq_t.comp, `TOST` = qq_tost.comp, `KS test` = qq_ks.comp, 
           AIC = aic_comp, `F-test` = f_comp)


# ## Do the QQ tests work better w/ null effects with larger samples? 
# ## No
# set.seed(2020-06-19*2)
# null_hi_power = many_metas(NN = 500, N = 20, delta = 0, 
#                            n = c(25, 250, 2500))
# 
# null_hi_power %>% 
#     select(delta, n, meta_idx, f_comp, aic_comp) %>% 
#     pivot_longer(c(f_comp, aic_comp), 
#                  names_to = 'test') %>% 
#     mutate(exp_value = case_when(delta == 0 & test == 'f_comp' ~ 'non-significant', 
#                                  delta != 0 & test == 'f_comp' ~ 'significant',
#                                  delta == 0 & test == 'aic_comp' ~ 'linear', 
#                                  delta != 0 & test == 'aic_comp' ~ 'quadratic',
#                                  TRUE ~ NA_character_), 
#            correct_call = value == exp_value) %>% 
#     group_by(delta, n, test) %>% 
#     summarize(accuracy = sum(correct_call) / n()) %>% 
#     ungroup() %>% 
#     arrange(delta, test, n)
# 
# ## Okay, so how about lots of studies? 
# ## Again, no. Accuracy apepars to go *down* with lots of studies
# set.seed(2020-06-19*3)
# null_many_reps = many_metas(NN = 200, N = c(20, 80, 140, 200), delta = 0, 
#                             n = 25)
# null_many_reps %>% 
#     select(delta, N, meta_idx, f_comp, aic_comp) %>% 
#     pivot_longer(c(f_comp, aic_comp), 
#                  names_to = 'test') %>% 
#     mutate(exp_value = case_when(delta == 0 & test == 'f_comp' ~ 'non-significant', 
#                                  delta != 0 & test == 'f_comp' ~ 'significant',
#                                  delta == 0 & test == 'aic_comp' ~ 'linear', 
#                                  delta != 0 & test == 'aic_comp' ~ 'quadratic',
#                                  TRUE ~ NA_character_), 
#            correct_call = value == exp_value) %>% 
#     group_by(delta, N, test) %>% 
#     summarize(accuracy = sum(correct_call) / n()) %>% 
#     ungroup() %>% 
#     arrange(delta, test, N)
# 
# ## Why would accuracy go down with more studies? 
# null_many_reps %>% 
#     group_by(N) %>% 
#     slice(1) %>% 
#     ungroup() %>% 
#     select(N, qq_slope, f_stat, f_p, f_comp, 
#            aic_linear, aic_quad)
# 
# null_many_reps %>% 
#     filter(meta_idx == 1) %>% 
#     select(-delta, -n) %>% 
#     unnest(studies) %>% 
#     qq_plot() +
#     facet_wrap(vars(N)) +
#     geom_smooth(method = 'lm', color = 'red') +
#     geom_smooth(method = 'lm', 
#                 formula = y ~ poly(x, 2, raw = TRUE), 
#                 color = 'blue')
# 
# null_many_reps %>% 
#     filter(meta_idx == 1) %>% 
#     select(-delta, -n) %>% 
#     unnest(studies) %>% 
#     ggplot(aes(p.value, color = N, group = N)) +
#     geom_density()


#' # Severity analysis
#+ severity analysis
## Severity analysis ----
## The cleanest way to specify and pass around the H0 and test output values is as `rlang` expressions.  
h_nought = exprs('0.2' = delta_fct == '0.2', 
                 'δ = 0.4' = delta_fct == '0.4', 
                 'δ = 0.6' = delta_fct == '0.6', 
                 'δ = 0.2, 0.4, 0.6' = delta_fct %in% c('0.2', '0.4', '0.6'), 
                 'δ is mixed' = delta_fct == 'mixed', 
                 'δ = 0' = delta_fct == '0', 
                 'δ is non-zero' = delta_fct %in% c('0.2', '0.4', '0.6', 'mixed'), 
                 'δ is not mixed' = delta_fct %in% c('0.2', '0.4', '0.6', '0'))
test_output = exprs(
    'ii-gap' = gappy == 'gappy',
    'iii-range' =  .9 < qq_estimate & qq_estimate < 1.1,
    'iii-T' = qq_t.comp == 'slope = 1',
    'iii-TOST' = qq_tost.comp == 'slope = 1',
    'iii-KS' = qq_ks.comp == 'uniform',
    'iv-AIC' = aic_comp == 'non-linear', 
    'iv-F' = f_comp == 'non-linear')

## These crosswalks let us attach short labels for hypotheses and outputs
h_xwalk = tibble(h_label = names(h_nought), 
                 h = map_chr(h_nought, as_label)) %>% 
    mutate(h_label = fct_inorder(h_label))
output_xwalk = tibble(output_label = names(test_output), 
                      test_output = map_chr(test_output, as_label)) %>% 
    mutate(output_label = fct_inorder(output_label))

## We first build a list containing each combination of hypothesis and output,
## then calculate the p-value for each combination. 
## The remaining steps just attach short labels and arrange the columns in a readable way
p_df = cross(lst(h_nought, test_output)) %>% 
    map_dfr(~p_value(!!.[['h_nought']], !!.[['test_output']], combined_df)) %>% 
    left_join(h_xwalk, by = 'h') %>% 
    rename(h_nought = h, h_nought_label = h_label) %>% 
    left_join(output_xwalk, by = 'test_output') %>% 
    select(h_nought_label, output_label, p, 
           h_nought, test_output, n_false, n_true)

ggplot(p_df, 
       aes(output_label, p)) +
    geom_point() +
    geom_area(data = tibble(output_label = 0:length(test_output) + 0.5),
              aes(y = .05), 
              alpha = .5) +
    labs(x = 'test output') +
    facet_wrap(vars(h_nought_label))

do.call(write_plot, c('fig_2_evidence_severity', plot_defaults))

p_df %>% 
    mutate(h_nought_label = str_replace(h_nought_label, 
                                        'δ', '$\\\\delta$')) %>% 
    select(`$H_0$` = h_nought_label, output = output_label, 
           p, false = n_false, true = n_true) %>% 
    arrange(`$H_0$`, output) %>% 
    kable(format = 'latex', booktabs = TRUE, longtable = TRUE, 
          digits = 2, 
          escape = FALSE,
          caption = 'Severity analysis results. $H_0$ indicates the null or rival hypothesis playing the role of $\\lnot H$. ``False" and ``true" indicate the number of simulation runs in which the test output is false and true, respectively, and p is calculated as the fraction of true runs.', 
          label = 'severity') %>% 
    write_lines(file.path(out_folder, 'severity.tex'))



#' # Likelihood analysis
#+ likelihood analysis
## Likelihood analysis ----
h1 = exprs('δ = 0' = delta_fct == '0', 
           'δ is mixed' = delta_fct == 'mixed')

h2 = exprs('δ = 0.2' = delta_fct == '0.2', 
           'δ = 0.4' = delta_fct == '0.4', 
           'δ = 0.6' = delta_fct == '0.6', 
           'δ = 0.2, 0.4, 0.6' = delta_fct %in% c('0.2', '0.4', '0.6'), 
           'δ is mixed' = delta_fct == 'mixed', 
           'δ = 0' = delta_fct == '0', 
           'δ is non-zero' = delta_fct %in% c('0.2', '0.4', '0.6', 'mixed'), 
           'δ is not mixed' = delta_fct %in% c('0.2', '0.4', '0.6', '0'))

## test_output is the same as for the severity analysis

## We can use the same crosswalks as above

## Calculate the likelihood ratios
llr_df = likelihood_ratio(h1, h2, test_output, combined_df) %>% 
    left_join(h_xwalk, by = c('h1' = 'h')) %>% 
    rename(h1_label = h_label) %>% 
    left_join(h_xwalk, by = c('h2' = 'h')) %>% 
    rename(h2_label = h_label) %>% 
    left_join(output_xwalk, by = 'test_output') %>% 
    filter((h1_label == 'δ is mixed' & str_detect(output_label, 'iv'))|
               (h1_label == 'δ = 0' & str_detect(output_label, 'iii'))) %>% 
    mutate(h1_label = fct_relevel(h1_label, 
                                  'δ = 0', 'δ is mixed')) %>% 
    select(h1_label, h2_label, output_label, llr, 
           h1, h2, test_output, everything())


ggplot(llr_df, 
       aes(output_label, llr, color = h2_label)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_rect(inherit.aes = FALSE,
              xmin = -Inf, xmax = Inf,
              ymin = 1, ymax = .5,
              alpha = .005) +
    geom_rect(inherit.aes = FALSE,
              xmin = -Inf, xmax = Inf,
              ymin = -1, ymax = -.5,
              alpha = .005) +
    geom_rect(inherit.aes = FALSE,
              xmin = -Inf, xmax = Inf,
              ymin = 1, ymax = 2,
              alpha = .01) +
    geom_rect(inherit.aes = FALSE,
              xmin = -Inf, xmax = Inf,
              ymin = -1, ymax = -2,
              alpha = .01) +
    geom_rect(inherit.aes = FALSE,
              xmin = -Inf, xmax = Inf,
              ymin = 2, ymax = Inf,
              alpha = .05) +
    geom_rect(inherit.aes = FALSE,
              xmin = -Inf, xmax = Inf,
              ymin = -Inf, ymax = -2,
              alpha = .05) +
    geom_point(size = 2, 
               position = position_jitter(width = .25, height = 0,
                                          seed = 2020-07-23)) +
    facet_wrap(vars(h1_label), scales = 'free_x') +
    ylim(-2, 2) +
    scale_color_viridis_d(option = 'C') +
    labs(color = 'H2', x = 'test output', y = 'log likelihood ratio')

ggplotly()

do.call(write_plot, c('fig_3_evidence_likelihood', plot_defaults))

## Test outcomes and H2 where the evidence is infinite
llr_df %>% 
    filter(is.infinite(llr)) %>% 
    select(h1_label, output_label, h2_label) %>% 
    arrange(h1_label, output_label, h2_label)

llr_df %>% 
    mutate(across(c(h1_label, h2_label), 
                  ~ str_replace(., 'δ', '$\\\\delta$'))) %>% 
    select(`$H_1$` = h1_label, `$H_2$` = h2_label, 
           output = output_label, llr, 
           `$L(H_1)$` = l_h1, `$L(H_2)$` = l_h2) %>% 
    arrange(`$H_1$`, `$H_2$`, output) %>% 
    kable(format = 'latex', booktabs = TRUE, longtable = TRUE, 
          digits = 2, 
          escape = FALSE,
          caption = 'Likelihood analysis results. ``llr" is the log likehood ratio.  Values $> 0.5$ indicate support for $H_1$.', 
          label = 'likelihood') %>% 
    write_lines(file.path(out_folder, 'likelihood.tex'))


#' # Power
#+ power
## Power ----
#' The supplemental analysis in this section shows how power influences the shape of Young's p-value plot and the outcomes of the KS tests for nonlinearity.  We examine three study designs, with 20%, 50%, and 80% power to detect a very small effect $\delta = .05$.  These are similar to the powers of the designs in the primary analysis (roughly 20%, 60%, and 90%).  

power.t.test(delta = .05, sd = 1, sig.level = .05, power = .2)
power.t.test(delta = .05, sd = 1, sig.level = .05, power = .5)
power.t.test(delta = .05, sd = 1, sig.level = .05, power = .8)

set.seed(2020-07-27)
tic()
power_analysis = many_metas(NN = 500, N = 25,
                            n = c(1002, 3075, 6281), 
                            delta = .05)
toc()


## Young's p-value plots
power_analysis %>% 
    group_by(n) %>%
    slice(1:5) %>%
    ungroup() %>%
    select(meta_idx, studies) %>% 
    unnest(studies) %>% 
    young_curve(color = as.factor(n)) +
    facet_grid(rows = vars(n), 
               cols = vars(meta_idx))

young_composite(power_analysis, alpha = .025) +
    facet_wrap(vars(n))

#' For the severely underpowered design, the KS test concludes that it is uniform nearly 75% of the time.  It almost always concludes that the other two designs are non-uniform.  
ggplot(power_analysis, aes(as.factor(n), 
                           fill = qq_ks.comp)) +
    geom_bar(position = 'fill')

#' Consequently, the KS-test can provide evidence for the zero hypothesis against the underpowered (50%) and adequately powered (80%) design.  But not against the severely underpowered design.  So the ability of the KS-test to provide evidence for the zero hypothesis depends on both real effect size and power.  
power_analysis %>% 
    split(.$n) %>% 
    map_dfr(~p_value(TRUE, qq_ks.comp == 'uniform', .), .id = 'n')

#' As power increases, both the t-test and TOST are increasingly likely to conclude that the slope is different from 1.  
power_analysis %>% 
    select(meta_idx, n, qq_t.comp, qq_tost.comp) %>% 
    pivot_longer(cols = c(qq_t.comp, qq_tost.comp), 
                 names_to = 'test', 
                 values_to = 'call') %>% 
    ggplot(aes(as.factor(n), fill = call)) +
    geom_bar(position = 'fill') +
    facet_wrap(vars(test))

#' As with the KS test, the t-test and TOST test can provide evidence for the zero hypothesis when power is greater, but not when power is low. 
power_analysis %>% 
    split(.$n) %>% 
    map_dfr(~p_value(TRUE, qq_t.comp == 'slope = 1', .), .id = 'n')
power_analysis %>% 
    split(.$n) %>% 
    map_dfr(~p_value(TRUE, qq_tost.comp == 'slope = 1', .), .id = 'n')


#' # Varying N #
#+ vary_N
## Varying N ----
#' The supplemental analysis in this section shows how varying the number of studies $N$ influences the shape of Young's p-value plot.  We use an underpowered study design with 60% power to detect a moderate effect $\delta = 0.4$. We use meta-analyses with 10, 50, 100, and 150 studies. 

power.t.test(delta = .4, sd = 1, sig.level = .05, power = .6)

set.seed(2020-10-16)
tic()
vary_N_analysis = many_metas(NN = 100, 
                             N = c(10, 50, 100, 150),
                             n = 62,
                             delta = 0.4)
toc()

#' Increasing $N$ simply effectively increases the meta-analysis' sampling from the underlying distribution of p-values for the study t-tests.  The result is that, as $N$ increases, the p-value curve gets filled in.  Among other things, this tends to decrease the size of the "gaps" in the blade of the "hockey stick." 
vary_N_analysis %>% 
    group_by(N) %>%
    slice(1:10) %>%
    ungroup() %>%
    select(meta_idx, studies, N) %>% 
    unnest(studies) %>% 
    young_curve(color = as.factor(N)) +
    # facet_grid(rows = vars(N), cols = vars(meta_idx), 
    #            scales = 'free_x')
    facet_wrap(vars(N, meta_idx), nrow = 3, 
               scales = 'free_x')

ggplot(vary_N_analysis, aes(as.factor(N), gap)) +
    geom_violin(draw_quantiles = .5) +
    geom_hline(yintercept = .125, alpha = .5)

ggplot(vary_N_analysis, aes(as.factor(N), fill = gappy)) +
    geom_bar(position = 'fill')


#' # Reproducibility
#+ reproducibility
## Reproducibility ----
sessionInfo()
