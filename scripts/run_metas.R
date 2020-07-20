library(tidyverse)
theme_set(theme_bw())
library(ggbeeswarm)
library(broom)
library(cowplot)

library(tictoc)


## TODO: 
## - "mixture" simulations
## - reorganize and rename
## - implement tests as tests

## Local package with simulation functions
devtools::load_all(file.path('..', 'p.curve'))
# devtools::document(file.path('..', 'p.curve'))
# source('test_scratch.R')


## Run simulations ----
## 4 effect sizes x 500 meta-studies x 20 studies x 25 samples
## 45 sec
set.seed(2020-06-25)
tic()
combined_df = many_metas(NN = 500, 
                         delta = c(0, .2, .4, .6), 
                         N = 20, n = 25) %>% 
    mutate(delta = as_factor(delta))
toc()

combined_df

## Sample plots ----
## Young
combined_df %>% 
    group_by(delta) %>% 
    slice(13, 17, 21, 25) %>% 
    ungroup() %>% 
    select(meta_idx, studies) %>% 
    unnest(studies) %>% 
    mutate(delta = as.factor(delta)) %>% 
    young_curve(color = delta) +
    facet_grid(rows = vars(delta), cols = vars(meta_idx), 
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
    select(meta_idx, studies) %>% 
    unnest(studies) %>% 
    mutate(delta = as.factor(delta)) %>% 
    schsp_curve(color = delta) +
    facet_grid(rows = vars(delta), cols = vars(meta_idx), 
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
    select(meta_idx, studies) %>% 
    unnest(studies) %>% 
    mutate(delta = as.factor(delta)) %>% 
    simonsohn_curve(color = delta) +
    facet_grid(rows = vars(delta), cols = vars(meta_idx), 
               as.table = TRUE,
               switch = 'y') +
    labs(x = 'p-value', 
         y = 'count') +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_brewer(palette = 'Set1', guide = FALSE)



## Slopes ----
combined_df %>% 
    select(delta, meta_idx, matches('slope')) %>% 
    pivot_longer(cols = matches('slope'), 
                 names_to = 'method', values_to = 'slope') %>% 
    ggplot(aes(delta, slope, color = delta)) +
    # geom_beeswarm(alpha = .10) +
    geom_violin(draw_quantiles = .5, fill = NA) +
    scale_color_brewer(palette = 'Set1', guide = FALSE) +
    facet_wrap(vars(method), scales = 'free')



combined_df %>% 
    group_by(delta) %>% 
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

test_levels = c('AIC: linear' = 'linear', 
                'AIC: quadratic' = 'quadratic', 
                'F-test: non-sig.' = 'non-significant', 
                'F-test: sig.' = 'significant')

combined_df %>% 
    select(delta, meta_idx, f_comp, aic_comp) %>% 
    pivot_longer(c(f_comp, aic_comp), 
                 names_to = 'test') %>% 
    count(delta, test, value) %>% 
    mutate(value = fct_relevel(value, !!!test_levels), 
           value = fct_recode(value, !!!test_levels), 
           test = fct_recode(test, 'AIC' = 'aic_comp', 'F-test' = 'f_comp')) %>% 
    ggplot(aes(delta, n, fill = value)) +
    geom_col(position = 'fill') +
    facet_wrap(vars(test), nrow = 2) +
    scale_fill_viridis_d() +
    scale_y_continuous(labels = scales::percent_format(), 
                       name = 'share of inferences')


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
