many_metas(delta = c(0, .2, .4, .6), 
           5, N = 14, n = 25) %>% 
    mutate(delta = as_factor(delta)) %>% 
    unnest(studies, names_repair = 'unique') %>% 
    # young_curve() +
    # schsp_curve() +
    facet_grid(rows = vars(delta...3), cols = vars(meta_idx), 
               as.table = FALSE,
               switch = 'y')


many_metas(delta = c(0, .2, .4, .6), 
           5, N = 14, n = 25) %>% 
    mutate(delta = as_factor(delta)) %>% 
    unnest(studies, names_repair = 'unique') %>% 
    ggplot(aes(p.value)) +
    geom_freqpoly(aes(y = after_stat(count)), 
                  binwidth = .01, 
                  boundary = 0.05,
                  pad = FALSE) +
    geom_rug() +
    geom_vline(xintercept = .05, linetype = 'dashed') +
    coord_cartesian(xlim = c(0, .06), 
                    expand = FALSE) +
    scale_y_continuous(name = 'Number of p-values') +
    facet_grid(rows = vars(delta...3), cols = vars(meta_idx), 
               as.table = FALSE,
               scales = 'free',
               switch = 'y')

