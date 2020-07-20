library(testthat)


delta = 0  ## true difference of means
n = 25 ## sample size of each study
N = 20 ## number of studies conducted

set.seed(2020-06-19)
samples_test = draw_samples(delta, n)
expect_equivalent(t.test(samples_test[[1]], 
                         samples_test[[2]])$statistic,
                  -0.09119745,
                  tolerance = .001, scale = 1)

t_test_exp = structure(list(estimate = -0.0243059191254392, 
                            estimate1 = -0.243587663558127, 
                            estimate2 = -0.219281744432687, 
                            statistic = c(t = -0.0911974540978669), 
                            p.value = 0.927772654726848, 
                            parameter = c(df = 41.6739457624658), 
                            conf.low = -0.562289172102317, 
                            conf.high = 0.513677333851439, 
                            method = "Welch Two Sample t-test",
                            alternative = "two.sided"), 
                       row.names = c(NA, -1L), 
                       class = c("tbl_df", "tbl", "data.frame"))

expect_equivalent(t_test(samples_test)$p.value[[1]], 
                  t_test_exp$p.value[[1]], 
                  tolerance = .0001, scale = 1)

set.seed(2020-06-19)
studies_null = draw_studies(N, delta = 0, n)

set.seed(2020-06-19)
studies_small = draw_studies(N, delta = .1, n)


## Young's p-curve ----
# young_curve(studies_null)
# young_curve(studies_small)

## These confirm that the seed is reproducing the results
expect_equivalent(young_slope(studies_null), 
                  0.04874478, tolerance = .0001, scale = 1)
expect_equivalent(young_slope(studies_small), 
                  0.04965476, tolerance = .0001, scale = 1)

## Confirm that young_slope() can recover slope correctly
slope_target = .629
data.frame(rank = 1:50) %>% 
    mutate(p.value = slope_target*rank) %>% 
    young_slope() %>% 
    expect_equal(slope_target, tolerance = .0001)


## QQ plot ----
# qq_plot(studies_null)
# qq_plot(studies_small)

## Confirm that qq_slope() can recover slope correctly
slope_target = .629
data.frame(rank = 1:50) %>% 
    mutate(p.uniform = rank / (max(rank) + 1),
           p.value = slope_target*p.uniform) %>% 
    qq_slope() %>% 
    expect_equal(slope_target, tolerance = .0001)


## Check qq_linear output
qq_linear_expect = structure(list(f_stat = 1.02990276511243, 
                                  alpha = 0.05, 
                                  f_p = 0.324408357117809, 
                                  f_comp = "non-significant", 
                                  aic_linear = -69.1177453155886, 
                                  aic_quad = -68.294111319479, 
                                  aic_comp = "linear"), 
                             row.names = c(NA, -1L), 
                             class = c("tbl_df", "tbl", "data.frame"))

expect_equal(qq_linear(studies_null)$f_stat[[1]], 
             qq_linear_expect$f_stat[[1]])
expect_equal(qq_linear(studies_null)$aic_quad[[1]], 
             qq_linear_expect$aic_quad[[1]])


## Schweder and Spjøtvoll's p-curve ----
# schsp_curve(studies_null)
# schsp_slope(studies_null)

## Confirm that schsp_slope() can recover slope correctly
## NB is 1/slope of p.value ~ rank
slope_target = 1/.629
data.frame(rank = 1:50) %>% 
    mutate(p.value = 1/slope_target*rank) %>% 
    schsp_slope() %>% 
    expect_equal(slope_target, tolerance = .0001)


## many_metas() ----
mm_test = many_metas(NN = 20, 
                 N = seq(10, 20, by = 10), 
                 n = seq(5, 10, by = 5), 
                 delta = c(0, 1))
expect_equal(nrow(mm_test), 
             20*2^3)


