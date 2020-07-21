library(tidyverse)
theme_set(theme_bw())
library(broom)
library(cowplot)

# library(ggbeeswarm)
# library(tidytext)

library(tictoc)


## TODO: 
## - docs: ex of how to pass just homogeneous effects vs passing mixed effects
## - reorganize and rename
## - implement tests as tests
## - CI to knit

## Local package with simulation functions
devtools::load_all(file.path('..', 'p.curve'))
# devtools::document(file.path('..', 'p.curve'))
# source('test_scratch.R')


## Run simulations ----
## 5 effect sizes x 500 meta-studies x 20 studies x 25 samples
## 45 sec
set.seed(2020-06-25)
tic()
combined_df = many_metas(NN = 500, 
                         delta = list(0, .2, .4, .6, list(0, .6)), 
                         N = 20, n = 25) %>% 
    mutate(delta_fct = as_factor(delta_chr))
toc()

combined_df


## Model validation ----
## Distribution of effects estimates
## This shows that, at the study level, the simulation gives an unbiased (ie, in expectation) estimate of the effect, except for the mixed case
combined_df %>% 
    select(delta_fct, meta_idx, studies) %>% 
    unnest(studies) %>% 
    ggplot(aes(delta_fct, estimate2, color = delta_fct)) +
    geom_violin(draw_quantiles = .5) +
    labs(x = 'real effect', 
         y = 'study-level effect estimate') +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

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

## Because the mean of the mean is just the mean, this gives the same values as the table above
# combined_df %>% 
#     select(delta_fct, meta_idx, studies) %>% 
#     unnest(studies) %>% 
#     group_by(delta_fct, meta_idx) %>% 
#     summarize(across(estimate2, lst(mean, sd))) %>% 
#     group_by(delta_fct) %>% 
#     summarize(across(estimate2_mean, lst(mean)), NN = n())



## Sample plots ----
## Young
combined_df %>% 
    group_by(delta) %>% 
    slice(13, 17, 21, 25) %>% 
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

## Sch-Sp
combined_df %>% 
    group_by(delta) %>% 
    slice(13, 17, 21, 25) %>% 
    ungroup() %>% 
    select(delta_fct, meta_idx, studies) %>% 
    unnest(studies) %>% 
    schsp_curve(color = delta_fct) +
    facet_grid(rows = vars(delta_fct), cols = vars(meta_idx), 
               as.table = TRUE,
               switch = 'y') +
    labs(x = '1 - p', 
         y = 'rank (descending)') +
    scale_color_brewer(palette = 'Set1', guide = FALSE)

## Simonsohn et al.
combined_df %>% 
    group_by(delta) %>% 
    slice(13, 17, 21, 25) %>% 
    ungroup() %>% 
    select(delta_fct, meta_idx, studies) %>% 
    unnest(studies) %>% 
    simonsohn_curve(color = delta_fct) +
    facet_grid(rows = vars(delta_fct), cols = vars(meta_idx), 
               as.table = TRUE,
               switch = 'y') +
    labs(x = 'p-value', 
         y = 'count') +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_brewer(palette = 'Set1', guide = FALSE)



## Slopes ----
combined_df %>% 
    select(delta_fct, meta_idx, matches('slope')) %>% 
    pivot_longer(cols = matches('slope'), 
                 names_to = 'method', values_to = 'slope') %>% 
    ggplot(aes(delta_fct, slope, color = delta_fct)) +
    # geom_beeswarm(alpha = .10) +
    geom_violin(draw_quantiles = .5, fill = NA) +
    scale_color_brewer(palette = 'Set1', guide = FALSE) +
    facet_wrap(vars(method), scales = 'free')



combined_df %>% 
    group_by(delta_fct) %>% 
    summarize_at(vars(matches('slope')), 
                 lst(median, sd)) #%>% view()



## There's not an exact relationship between Young and Sch-Sp slopes, 
## even when N (and so the max rank) is constant
ggplot(combined_df, aes(young_slope, schsp_slope)) +
    geom_point() +
    scale_y_log10()

## Whereas Young slope is a rescaling of QQ slope; 
## the x-axis in the QQ plot is rank/max(rank), and rank is the x-axis in the Young plot
ggplot(combined_df, aes(young_slope, qq_slope)) +
    geom_point()


## QQ linearity tests ----
## TODO: this is awful
# combined_df %>% 
#     select(delta, meta_idx, f_comp, aic_comp) %>% 
#     pivot_longer(c(f_comp, aic_comp), 
#                  names_to = 'test') %>% 
#     group_by(test, delta) %>% 
#     count(value) %>% 
#     mutate(share = n / sum(n)) %>% 
#     ggplot(aes(test, share, fill = value)) +
#     geom_col() +
#     scale_fill_manual(values = c('linear' = 'red', 
#                                  'quadratic' = 'blue', 
#                                  'non-significant' = 'green', 
#                                  'significant' = 'yellow')) +
#     facet_wrap(vars(delta))

test_levels = c('AIC: linear' = 'linear', 
                'AIC: quadratic' = 'quadratic', 
                'F-test: non-sig.' = 'non-significant', 
                'F-test: sig.' = 'significant')

combined_df %>% 
    select(delta_fct, meta_idx, f_comp, aic_comp) %>% 
    pivot_longer(c(f_comp, aic_comp), 
                 names_to = 'test') %>% 
    count(delta_fct, test, value) %>% 
    mutate(value = fct_relevel(value, !!!test_levels), 
           value = fct_recode(value, !!!test_levels), 
           test = fct_recode(test, 'AIC' = 'aic_comp', 'F-test' = 'f_comp')) %>% 
    ggplot(aes(delta_fct, n, fill = value)) +
    geom_col(position = 'fill') +
    facet_wrap(vars(test), nrow = 2) +
    scale_fill_viridis_d() +
    xlab('real effect') +
    scale_y_continuous(labels = scales::percent_format(), 
                       name = 'share of inferences')


combined_df %>% 
    select(delta_fct, meta_idx, f_comp, aic_comp) %>% 
    pivot_longer(c(f_comp, aic_comp), 
                 names_to = 'test') %>% 
    mutate(exp_value = case_when(delta_fct == 0 & test == 'f_comp' ~ 'non-significant', 
                                 delta_fct != 0 & test == 'f_comp' ~ 'significant',
                                 delta_fct == 0 & test == 'aic_comp' ~ 'linear', 
                                 delta_fct != 0 & test == 'aic_comp' ~ 'quadratic',
                                 TRUE ~ NA_character_), 
           correct_call = value == exp_value) %>% 
    group_by(delta_fct, test) %>% 
    summarize(accuracy = sum(correct_call) / n()) %>% 
    ungroup() %>% 
    pivot_wider(names_from = test, values_from = accuracy) %>% 
    rename(AIC = aic_comp, `F-test` = f_comp)


## Do the QQ tests work better w/ null effects with larger samples? 
## No
set.seed(2020-06-19*2)
null_hi_power = many_metas(NN = 500, N = 20, delta = 0, 
                           n = c(25, 250, 2500))

null_hi_power %>% 
    select(delta, n, meta_idx, f_comp, aic_comp) %>% 
    pivot_longer(c(f_comp, aic_comp), 
                 names_to = 'test') %>% 
    mutate(exp_value = case_when(delta == 0 & test == 'f_comp' ~ 'non-significant', 
                                 delta != 0 & test == 'f_comp' ~ 'significant',
                                 delta == 0 & test == 'aic_comp' ~ 'linear', 
                                 delta != 0 & test == 'aic_comp' ~ 'quadratic',
                                 TRUE ~ NA_character_), 
           correct_call = value == exp_value) %>% 
    group_by(delta, n, test) %>% 
    summarize(accuracy = sum(correct_call) / n()) %>% 
    ungroup() %>% 
    arrange(delta, test, n)

## Okay, so how about lots of studies? 
## Again, no. Accuracy apepars to go *down* with lots of studies
set.seed(2020-06-19*3)
null_many_reps = many_metas(NN = 200, N = c(20, 80, 140, 200), delta = 0, 
                            n = 25)
null_many_reps %>% 
    select(delta, N, meta_idx, f_comp, aic_comp) %>% 
    pivot_longer(c(f_comp, aic_comp), 
                 names_to = 'test') %>% 
    mutate(exp_value = case_when(delta == 0 & test == 'f_comp' ~ 'non-significant', 
                                 delta != 0 & test == 'f_comp' ~ 'significant',
                                 delta == 0 & test == 'aic_comp' ~ 'linear', 
                                 delta != 0 & test == 'aic_comp' ~ 'quadratic',
                                 TRUE ~ NA_character_), 
           correct_call = value == exp_value) %>% 
    group_by(delta, N, test) %>% 
    summarize(accuracy = sum(correct_call) / n()) %>% 
    ungroup() %>% 
    arrange(delta, test, N)

## Why would accuracy go down with more studies? 
null_many_reps %>% 
    group_by(N) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(N, qq_slope, f_stat, f_p, f_comp, 
           aic_linear, aic_quad)

null_many_reps %>% 
    filter(meta_idx == 1) %>% 
    select(-delta, -n) %>% 
    unnest(studies) %>% 
    qq_plot() +
    facet_wrap(vars(N)) +
    geom_smooth(method = 'lm', color = 'red') +
    geom_smooth(method = 'lm', 
                formula = y ~ poly(x, 2, raw = TRUE), 
                color = 'blue')

null_many_reps %>% 
    filter(meta_idx == 1) %>% 
    select(-delta, -n) %>% 
    unnest(studies) %>% 
    ggplot(aes(p.value, color = N, group = N)) +
    geom_density()
