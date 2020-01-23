library(tidyverse)
library(broom)

delta = 0  ## true difference of means
n = 25 ## sample size of each study
N = 20 ## number of studies conducted


## Functions to draw a sample ("run one study"), do a t-test, and draw multiple samples ("do a meta-analysis") ----
draw_samples = function(delta, n) {
    group1 = rnorm(n, mean = 0)
    group2 = rnorm(n, mean = delta)
    return(list(group1, group2))
}

# foo = draw_samples(delta, n)

t_test = function(samples) {
    test = t.test(samples[[1]], samples[[2]])
    return(tidy(test))
}

# t_test(foo)

draw_studies = function(N, delta, n, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    tibble(study_id = 1:N) %>% 
        mutate(delta = delta, 
               n = n) %>% 
        mutate(samples = map2(delta, n, draw_samples)) %>% 
        mutate(test = map(samples, t_test)) %>% 
        unnest(test) %>% 
        mutate(rank = rank(p.value), 
               p.uniform = rank / (max(rank) + 1))
}

# studies = draw_studies(N, delta, n)
# studies_small = draw_studies(N, delta = .1, n)


## Young's version of the p-curve ----
young_curve = function(studies, draw_alpha = TRUE) {
    plot = ggplot(studies, aes(rank, p.value)) +
        geom_point()
    if (draw_alpha) {
        plot = plot + geom_hline(yintercept = .05)
    }
    return(plot)
}

# young_curve(studies)
# young_curve(studies_small)

young_slope = function(studies) {
    model = lm(p.value ~ rank, data = studies)
    slope = model$coefficients[['rank']]
    return(slope)
}


## QQ plot of p-values against a uniform distribution ----
qq_plot = function(studies) {
    ggplot(studies, aes(p.uniform, p.value)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0)
}

# qq_plot(studies)
# qq_plot(studies_small)

qq_slope = function(studies) {
    model = lm(p.value ~ p.uniform, data = studies)
    slope = model$coefficients[['p.uniform']]
    return(slope)
}
# qq_slope(studies)

qq_linear = function(studies, alpha = .05) {
    ## Checks linearity of the QQ plot by comparing linear and quadratic fits, using (a) F-test and (b) AICs
    model_linear = lm(p.value ~ p.uniform, data = studies)
    model_quad = lm(p.value ~ p.uniform + I(p.uniform^2), 
                    data = studies)
    
    f_comp = (anova(model_linear, model_quad)[['Pr(>F)']][[2]]) < alpha
    aic_comp = AIC(model_quad) < AIC(model_linear)
    return(tibble(f_comp, aic_comp))
}


## Simonsohn et al.'s p-curve ----
# studies_med$studies[[2]] %>%
#     # filter(p.value < .05) %>%
#     ggplot(aes(p.value)) +
#     stat_bin(binwidth = .01, geom = 'point') +
#     stat_bin(binwidth = .01, geom = 'line')

## Schweder and Spjøtvoll's p-curve ----
schsp_curve = function(studies) {
    ggplot(studies, aes(1 - p.value, 
                        max(rank) - rank)) +
        geom_point()
}
# schsp_curve(studies)

## Schweder and Spjøtvoll slope
schsp_slope = function(studies) {
    model = studies %>% 
        mutate(rank.rev = max(rank) - rank, 
               p.value.rev = 1 - p.value) %>% 
        lm(rank.rev ~ p.value.rev, data = .)
    slope = model$coefficients[['p.value.rev']]
    return(slope)
}
# schsp_slope(studies)


## Do many meta-analyses ----
many_metas = function(NN, ## how many meta-studies
                      N, delta, n) {
    tibble(meta_id = 1:NN) %>% 
        mutate(N, delta, n) %>% 
        mutate(studies = pmap(list(N, delta, n), draw_studies)) %>% 
        mutate(young_slope = map_dbl(studies, young_slope), 
               schsp_slope = map_dbl(studies, schsp_slope),
               qq_slope = map_dbl(studies, qq_slope), 
               qq_linear = map(studies, qq_linear)) %>% 
        unnest(qq_linear)
}

studies_null = many_metas(250, N, delta = 0, n)
studies_small = many_metas(250, N, delta = 0.2, n)
studies_med = many_metas(250, N, delta = 0.5, n)

combined_df = bind_rows(list(null = studies_null, 
                             small = studies_small, 
                             medium = studies_med),
                        .id = 'true_effect')

combined_df %>% 
    select(true_effect, meta_id, matches('slope')) %>% 
    pivot_longer(cols = matches('slope'), 
                 names_to = 'method', values_to = 'slope') %>% 
    ggplot(aes(slope, color = true_effect)) +
    geom_density() +
    facet_wrap(vars(method), scales = 'free')

# qq_plot(studies_med$studies[[1]])
# qq_plot(studies_small$studies[[1]])
# qq_plot(studies_null$studies[[1]])

# schsp_curve(studies_small$studies[[1]])
# young_curve(studies_med$studies[[1]], draw_alpha = FALSE)

combined_df %>% 
    group_by(true_effect) %>% 
    summarize_at(vars(matches('slope')), 
                 list(~median(.), ~sd(.))) #%>% view()


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
    geom_point()

## Whereas Young slope is a rescaling of QQ slope; 
## the x-axis in the QQ plot is rank/max(rank), and rank is the x-axis in the Young plot
ggplot(combined_df, aes(young_slope, qq_slope)) +
    geom_point()

## How good are the QQ linearity checks? 
combined_df %>% 
    select(true_effect, meta_id, f_comp, aic_comp) %>% 
    pivot_longer(c(f_comp, aic_comp), 
                 names_to = 'test') %>% 
    group_by(test, true_effect) %>% 
    count(value) %>% 
    mutate(share = n / sum(n)) %>% 
    ggplot(aes(test, share, fill = value)) +
    geom_col() +
    facet_wrap(vars(true_effect))
