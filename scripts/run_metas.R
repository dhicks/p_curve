library(tidyverse)
theme_set(theme_bw())
library(ggbeeswarm)
library(broom)
library(cowplot)

library(tictoc)


## TODO: 
## - "mixture" simulations
## - reorganize and rename
## - finishing documenting p-curve.R
## - implement tests as tests

## Local package with simulation functions
devtools::load_all(reset = TRUE, recompile = TRUE, file.path('..', 'p.curve'))
?qq_linear
# source('test_scratch.R')


## Run simulations ----
## 4 effect sizes x 500 meta-studies x 20 studies x 25 samples
## 45 sec
set.seed(2020-06-25)
tic()
combined_df = many_metas(delta = c(0, .2, .4, .6), 
                         NN = 5, N = 20, n = 25) %>% 
    mutate(delta = as_factor(delta))
toc()


## Slopes ----
combined_df %>% 
    select(delta, meta_idx, matches('slope')) %>% 
    pivot_longer(cols = matches('slope'), 
                 names_to = 'method', values_to = 'slope') %>% 
    ggplot(aes(delta, slope, color = delta)) +
    # geom_beeswarm(alpha = .10) +
    geom_violin(draw_quantiles = .5, fill = NA) +
    facet_wrap(vars(method), scales = 'free')

## TODO: do this with facets rather than nested cowplots
combined_df %>% 
    group_by(delta) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(qq_plot = map(studies, qq_plot), 
           schsp_curve = map(studies, schsp_curve), 
           young_curve = map(studies, young_curve)) %>% 
    select(delta, meta_idx, qq_plot, schsp_curve, young_curve) %>% 
    pivot_longer(cols = c(qq_plot, schsp_curve, young_curve), 
                 names_to = 'method', values_to = 'plot') %>% 
    group_by(method) %>% 
    summarize(plot = list(plot_grid(plotlist = plot, 
                                    labels = delta))) %>% 
    ungroup() %>% 
    summarize(plot = list(plot_grid(plotlist = plot, 
                                    labels = method))) %>% 
    pull(plot)



combined_df %>% 
    group_by(delta) %>% 
    summarize_at(vars(matches('slope')), 
                 lst(median, sd)) #%>% view()


# ggplot(combined_df, aes(schsp_slope, color = true_effect)) +
#     geom_density()
# 
# combined_df %>% 
#     group_by(true_effect) %>% 
#     summarize_at(vars(schsp_slope), 
#                  funs(median, sd))

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
combined_df %>% 
    select(delta, meta_idx, f_comp, aic_comp) %>% 
    pivot_longer(c(f_comp, aic_comp), 
                 names_to = 'test') %>% 
    group_by(test, delta) %>% 
    count(value) %>% 
    mutate(share = n / sum(n)) %>% 
    ggplot(aes(test, share, fill = value)) +
    geom_col() +
    scale_fill_manual(values = c('linear' = 'red', 
                                 'quadratic' = 'blue', 
                                 'non-significant' = 'green', 
                                 'significant' = 'yellow')) +
    facet_wrap(vars(delta))

combined_df %>% 
    select(delta, meta_idx, f_comp, aic_comp) %>% 
    pivot_longer(c(f_comp, aic_comp), 
                 names_to = 'test') %>% 
    mutate(exp_value = case_when(delta == 0 & test == 'f_comp' ~ 'non-significant', 
                                 delta != 0 & test == 'f_comp' ~ 'significant',
                                 delta == 0 & test == 'aic_comp' ~ 'linear', 
                                 delta != 0 & test == 'aic_comp' ~ 'quadratic',
                                 TRUE ~ NA_character_), 
           correct_call = value == exp_value) %>% 
    group_by(delta, test) %>% 
    summarize(accuracy = sum(correct_call) / n())


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
