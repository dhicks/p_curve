---
title: "Young's p-value plot does not provide evidence against air pollution epidemiology"
author: "Daniel J. Hicks"
institution: "University of California, Merced"
email: hicks.daniel.j@gmail.com
abstract: foo
bibliography: Young.bib
---

# Introduction #

In a recent series of papers, statistician S. Stanley Young and collaborators *[cites]* have criticized air pollution epidemiology, arguing that a graphical method they call a "p-value plot" undermines these studies and, by extension, the scientific basis for regulating air pollution *[cite page]*.  The current paper uses simulation methods to show that the p-value plot method does not have the statistical properties required to support the criticisms leveled by Young and collaborators.  Therefore Young's p-value plot does not provide evidence against air pollution epidemiology.  

*[background: 
- air pollution epidemiology
- Young's two arguments]*


# Methods #

I first discuss Young's p-value plot, the graphical method used as the basis of the second argument by Young and collaborators against air pollution epidemiology.  I then turn to the simulation methods used in the current study to analyze the statistical properties of Young's p-value plot.  

Every peer-reviewed publication describing the method has been written by a collection of authors.  In each of these works the method is referred to as a "p-value plot," without attributing it to any individual.  However, Young is the only author who appears consistently across all of these papers, and the method needs to be distinguished from other, nominally similar methods.  I therefore refer to the method as "Young's p-value plot."  

## Young's p-value plot ##

Young and collaborators do not provide a formal definition of their p-value plot method, and at least one key term in the methods appears to be used in a non-standard way.  They do attempt to justify their method by associating it with two other graphical methods, but these methods are not equivalent to each other and do not appear to be equivalent to Young's p-value plot.  

\renewcommand{\P}{\ensuremath{\mathbb{P}}}

All of these methods take as input a set of $N$ p-values $\P = \{p_1, p_2, \ldots, p_N\}$, nominally produced by applying a given hypothesis statistical test to $N$ replications of a given study design, each replication drawing samples of size $n$ from a given population.  This corresponds to the simplest case of meta-analysis.  Thus the p-values in $\P$ are nominally samples from a single underlying distribution $p_i \sim P$.  Note that, if the null hypothesis $H_0$ for the study design is true, then $P$ is the uniform distribution on $[0,1]$.  

### Schweder and Spjøtvoll's p-value plot ###

Young and collaborators frequently cite the "p-value plot" presented in @SchwederPlotsPvaluesEvaluate1982.  For this p-value plot, let $rank_{desc}(p_i)$ be the (1-indexed) *descending rank* of $p_i \in \P$, i.e., $rank_{desc}(p_i)$ is the number of p-values $p_j \in \P$ greater than or equal to $p_i$.  The largest p-value has descending rank 1, and the smallest p-value has descending rank $N$.  Then Schweder and Spjøtvoll's p-value plot plots the graph $(1-p_i, rank_{desc}(p_i))$.  See figure *[example + cite]*

@SchwederPlotsPvaluesEvaluate1982 briefly argue that the relationship between $1-p_i$ and $rank_{desc}(p_i)$ should be approximately linear "when [$p_i$] is not too small" and that, from left to right, "often, the plot will not show a clearcut break but rather a gradual bend" away from linearity.  Note that this means Schweder and Spjøtvoll's p-value plot generally ignores small p-values, which roughly correspond to the statistically significant p-values.  They illustrate their method with example data and a straight line "drawn by visual fit" rather than regression analysis.    

### Simonsohn, Nelson, and Simmons' p-curve ###

Young and collaborators also echo concerns about the replication crisis unfolding in social psychology and certain areas of biomedical research *[cite]*.  In particular, they appeal to concerns about "p-hacking," in which researchers conduct variant analyses (e.g., including or excluding various controls in a regression analysis) and report only the analyses with statistically significant p-values.  Note that, other than Young's p-value plot, Young and collaborators provide no empirical evidence of p-hacking in environmental epidemiology, and to my knowledge no such evidence has been published.  

Young and collaborators associate Young's p-value plot with a method developed to detect p-hacking, called a "p-curve" [@SimonsohnPcurveKeyFiledrawer2014].  The intuition behind the p-curve is that p-hacking will tend to produce an excess number of p-values "just below" the conventional $0.05$ threshold for statistical significance.  Formally, they divide the interval $[0, 0.05]$ into 5 bins, divided by the thresholds $0.01, 0.02, 0.03, 0.04, 0.05$, then calculate $N_b$, the number of p-values in bin $b$.  The p-curve is the graph $(b_t, N_b)$, where $b_t$ is the threshold for bin $b$.  @SimonsohnPcurveKeyFiledrawer2014 then formally test for p-hacking by applying statistical tests of the null hypothesis that the restricted distribution $P|_{p < 0.05}$ is uniform.  See figure *[example + cite]*

Note that (1) the p-curve is equivalent to a histogram on the interval $[0, 0.5]$ with binwidth $0.01$ and (2) the p-curve only includes statistically significant p-values.  Thus the p-curve and Schweder and Spjøtvoll's p-value plot not only produce different kinds of plots from the input data $\P$, but actually direct their attention to different — typically disjoint — subsets of $\P$.  The two kinds of plots cannot be equivalent to each other.  At no point do Young and collaborators acknowledge this fundamental difference between the two methods that they cite as support for their own method.  

### Young's p-value plot ###

For Young's p-value plot, let $rank_{asc}(p_i)$ be the (1-indexed) *ascending rank* of $p_i \in \P$, i.e., $rank_{asc}(p_i)$ is the number of p-values $p_j \in \P$ less than or equal to $p_i$.  The smallest p-value has ascending rank 1, and the largest p-value has ascending rank $N$.  Without loss of generality, if \P is already in ascending order $p_1 < p_2 < \cdots < p_N$, then $rank_{asc}(p_i) = i$.  And Young's p-value plot is the graph $(i, p_i)$.  See figure *[example + cite]*

Note that Young's p-value plot is a rescaled QQ-plot of $\P$ against the uniform distribution, with the theoretical quantiles $q_i = \frac{i}{N} = \frac{rank_{asc}(p_i)}{N}$.  In addition, $rank_{desc}(p_i) = N - rank_{asc}(p_i) + 1$, and so Young's p-value does have a 1-1 mathematical relationship to Schweder and Spjøtvoll's p-value plot.  However, the nonlinear transformation (switching axes) means that a regression line fit to Schweder and Spjøtvoll's p-value plot will not necessarily translate to a regression line fit to Young's p-value plot.  

Young and collaborators give the following method to analyze Young's p-value plot:  "Evaluation of a p-value plot follows a logical path. Is the p-value plot homogeneous? If the points roughly [sic] on a 45-degree, they are homogeneous and consistent with randomness; a lessor slope with all points roughly on a line indicates a consistent effect even if some of the individual p-values are not considered statistically significant. If the effects differ, one from another, beyond chance, then the effects are heterogeneous" [@YoungReliabilityEnvironmentalEpidemiology2019 48].  None of the key terms here are given formal definitions.  

A "45-degree line" typically refers to the graph of an identity function or the line $y = x$, which forms a 45-degree angle with respect to both the x- and y-axes.  This interpretation makes sense for a QQ-plot of p-values against the uniform distribution on $[0, 1]$, where both axes are on the scale $[0, 1]$.  This interpretation does not make sense for Young's p-value plot, where the x-axis (rank) is on the scale $[1, N]$ and the y-axis (p-value) is on the scale $[0, 1]$.  Instead, Young and collaborators seem to use the term "45-degree line" when the p-value plot is approximately linear [for example, @YoungCombinedBackgroundInformation2019 23].  This interpretation is consistent with one figure by Young and collaborators [@YoungCombinedBackgroundInformation2019 18], which seems to indicate the use of a F-test to compare a quadratic regression against a linear regression, which is interpreted as indicating that there is a "bilinear effect."  Context suggests that "bilinear effect" is synonymous with "heterogeneous" in the quotation above; this language suggests a case in which the p-values \P are produced by two separate processes, e.g., the effect under study is null in one population and non-null in another.  (Note that this is conceptually quite different from a single non-linear effect.)  This "heterogeneity" also appears to be connected to "a hockey stick pattern, a number of small p-values followed by p-values falling roughly on a 45-degree line" [@YoungCombinedBackgroundInformation2019 9].  On the other hand, this F-test method is not used in other analyses.  

Alternatively, the "45-degree line" might mean that the slope of the regression line is $1/N$, corresponding to the fact that Young's p-value plot is a rescaled QQ-plot.  (In the corresponding QQ-plot, a slope of 1 with good fit would indicate that the p-values are approximately uniformly distributed, which in turn would indicate that the null hypothesis of the underlying study design is true.)  However, other than the single use of a F-test cited in the previous paragraph, Young and collaborators do not apply regression or indeed any other formal analysis to the p-value plots.  Instead, analyses appear to be based on (1) a visual judgment of whether the curve is linear or not, and (2) the number or fraction of statistically significant p-values.  

Young and collaborators do not appear to have validated the p-value plot method in any peer-reviewed studies.  @YoungCombinedBackgroundInformation2019 (8-9) briefly report a simulation study, specified as:  "We simulate 15 p-values ten times and produce the p-value plots, Figure 1. Each of the 10 simulated meta-analysis studies has 15 base papers .... The goal is to get an impression of the variability of a p-value plot under ideal conditions where there is nothing going on and there is no multiple testing or multiple modeling in the base studies."  In my notation, each run of the simulation has $N = 15$ p-values drawn from a uniform distribution, i.e., the null hypothesis is true.  They conduct $NN = 10$ simulation runs, but do not report any aggregate or distributional statistics; instead they report the smallest (rank = 1) p-value and the p-value plot for each simulation run.  Note that they do not report any simulations of cases where the null hypothesis is false or the set of p-values \P is produced by two distinct processes.  They do not report making the simulation code available anywhere for checking reproducibility or extending their analysis.  


## The current simulation study ##

The current study is an automatically-reproduced simulation study of Young's p-value plot, systematically testing it across three kinds of cases, using formal analyses of the graphs, and comparing it to Schweder and Spjøtvoll's p-value plot and Simonsohn, Nelson, and Simmons' p-curve.  

Each run of the simulation is composed of $N$ studies.  To make the study parameters easy to understand, each study is based on a two-sample t-test.  Two samples, each of size $n$, are drawn from Gaussian distributions with mean 0 and $\delta$, respectively, and common standard deviation $\sigma$.  In the basic "single process" case, $\delta$ is the same for all studies in the simulation run, and so $\delta$ corresponds to the true effect size of the (single) process being studied.  

A *condition* of the simulation fixes the values of these parameters.  Conditions can be systematically varied to compare, e.g., different effect sizes, and multiple runs $NN$ in each condition allow us to analyze the frequentist properties of Young's p-value plot within and across conditions.  In particular, 4 different effect sizes are used — corresponding to null, small, moderate, and moderate-large effects — and each condition is simulated with 500 runs.  500 runs per condition was chosen because this should produce good estimates of central tendency and variation without requiring significant computation time.  Note that 500 runs is probably too few to produce good estimates of extreme values, e.g., 5th percentiles.  The primary conditions examined in the current study are summarized in [@tbl:params].  

| parameter | meaning                             | value(s) |
|:----------|:------------------------------------|---------:|
| $\delta$  | study effect size                   | 0, .2, .4, .6 |
| $\sigma$  | s.d. of study samples               | 1        |
| $n$       | study sample size                   | 25       |
| $N$       | num. studies in each run            | 20       |
| $NN$      | num. runs for each condition        | 500      |

Table: Base parameter values used in the current simulation study {#tbl:params}





