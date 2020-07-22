---
title: "Young's p-value plot does not provide evidence against air pollution hazards"
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
- Young's two arguments
- two hypotheses:  1. null effect (the pollutant has no effect on the population); 2. heterogeneous effects (the population has two subpopulations; the pollutant has an effect on only one of the two subpopulations)]*

For the purpose of evaluating Young's p-value plot, simulations offer a number of advantages over exact methods.  First, for the target audience of this paper, well-written R code is likely to be more accessible than an exact analysis of the statistical properties of meta-analyses.  This makes it easier for readers to scrutinize the details of my argument.  Second, the simulation is modular, making it more extensible or modifiable by readers than an exact analysis.  The simulated data generating process and its analysis can be modified by simply rewriting two short functions and re-running the simulation.  An exact analysis of a different data-generating process might require a completely different formal approach.  Third, it may be algebraically intractable to give exact characterizations of constructions such as plots except in the limit (as sample size and/or number of studies goes to infinity).  


# Methods #

I first discuss Young's p-value plot, the graphical method used as the basis of the second argument by Young and collaborators against air pollution epidemiology.  I then turn to the simulation methods used in the current study to analyze the statistical properties of Young's p-value plot.  

Every peer-reviewed publication describing the method has been written by a collection of authors.  In each of these works the method is referred to as a "p-value plot," without attributing it to any individual.  However, Young is the only author who appears consistently across all of these papers, and the method needs to be distinguished from other, nominally similar methods.  I therefore refer to the method as "Young's p-value plot."  

## Young's p-value plot ##

Young and collaborators do not provide a formal definition of their p-value plot method, and at least one key term in the methods appears to be used in a non-standard way.  They do attempt to justify their method by associating it with two other graphical methods, but these methods are not equivalent to each other and do not appear to be equivalent to Young's p-value plot.  

\renewcommand{\P}{\ensuremath{\mathbb{P}}}

All of these methods take as input a set of $N$ p-values $\P = \{p_1, p_2, \ldots, p_N\}$, nominally produced by applying a given hypothesis statistical test to $N$ replications of a given study design, each replication drawing samples of size $n$ from a given population.  This corresponds to the simplest case of meta-analysis.  Thus the p-values in $\P$ are nominally samples from a single underlying distribution $p_i \sim P$.  Note that, if the null hypothesis $H_0$ for the study design is true, then $P$ is the uniform distribution on $[0,1]$.  

### Schweder and Spjøtvoll's p-value plot ###

\newcommand{\SchSp}{Schweder and Spjøtvoll}

Young and collaborators frequently cite the "p-value plot" presented in @SchwederPlotsPvaluesEvaluate1982.  For this p-value plot, let $rank_{desc}(p_i)$ be the (1-indexed) *descending rank* of $p_i \in \P$, i.e., $rank_{desc}(p_i)$ is the number of p-values $p_j \in \P$ greater than or equal to $p_i$.  The largest p-value has descending rank 1, and the smallest p-value has descending rank $N$.  Then \SchSp's p-value plot plots the graph $(1-p_i, rank_{desc}(p_i))$.  See figure *[example + cite]*

@SchwederPlotsPvaluesEvaluate1982 briefly argue that the relationship between $1-p_i$ and $rank_{desc}(p_i)$ should be approximately linear "when [$p_i$] is not too small" and that, from left to right, "often, the plot will not show a clearcut break but rather a gradual bend" away from linearity.  Note that this means \SchSp's p-value plot generally ignores small p-values, which roughly correspond to the statistically significant p-values.  They illustrate their method with example data and a straight line "drawn by visual fit" rather than regression analysis.    

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

The current study is an automatically-reproduced simulation study of Young's p-value plot, systematically testing it across five kinds of cases, using formal analyses of the graphs, and comparing it to Schweder and Spjøtvoll's p-value plot and Simonsohn, Nelson, and Simmons' p-curve.  

Each run of the simulation is composed of $N$ studies.  To make the study parameters easy to understand, each study is based on a two-sample t-test.  Two samples, each of size $n$, are drawn from Gaussian distributions with mean 0 and $\delta$, respectively, and common standard deviation $\sigma$.  In the basic "single process" case, $\delta$ is the same for all studies in the simulation run, and so $\delta$ corresponds to the true effect size of the (single) process being studied.  

A *condition* of the simulation fixes the values of these parameters.  Conditions can be systematically varied to compare, e.g., different effect sizes, and multiple runs $NN$ in each condition allow us to analyze the frequentist properties of Young's p-value plot within and across conditions.  In particular, 5 different effect sizes are used — corresponding to null, small, moderate, moderate-large effects, and a "mixed" or heterogeneous condition in which half the population has a null response and half the population has a moderate-large response — and each condition is simulated with 500 runs.  500 runs per condition was chosen because this should produce good estimates of central tendency and variation without requiring significant computation time.  Note that 500 runs is probably too few to produce good estimates of extreme values, e.g., 5th percentiles.  The primary conditions examined in the current study are summarized in [@tbl:params].  

| parameter | meaning                             | value(s) |
|:----------|:------------------------------------|---------:|
| $\delta$  | real effect size                    | 0, .2, .4, .6, (0, .6) |
| $\sigma$  | s.d. of samples                     | 1        |
| $n$       | study sample size                   | 25       |
| $N$       | num. studies in each run            | 20       |
| $NN$      | num. runs for each condition        | 500      |

Table: Base parameter values used in the current simulation study.  The real effect size of (0, .6) indicates the "mixed" or heterogeneous condition, in which half the population have no response and half the population have a strong response. {#tbl:params}

The simulation can also simulate a "mixed," multi-process, or heterogeneous condition, in which subsets of the population have different responses to the intervention or exposure.  The current study uses a (0, .6) mixture, meaning that half of the population has a zero effect or no response and half the population has a medium-strong response of .6.  These values were chosen to create a relatively dramatic or striking mixture that is less likely to look like a single-process or homogeneous case.  For this mixed case, the simulation draws $n$ control/non-exposed group samples (as with the single process case) and then, for the intervention/exposed group, draws $n/2$ samples from the no response subpopulation and $n/2$ samples from medium-strong response subpopulation.  This guarantees that the sample is balanced across the two (latent) subpopulations.  I expect that unbalanced samples will tend to produce less striking p-curves, as unbalanced samples will tend to look more like the majority subpopulation.  (This approach to modeling a "mixed" condition can be considered a constrained and simplified version of a Gaussian mixture *[cite]*.  So a generalization of the simulation used here might start instead by drawing from a standard Gaussian mixture, and treating single-process conditions as limiting cases with a single mixture component.)   

After generating the data, it is straightforward to construct the three graphs (defined here formally as a set of ordered pairs).  I then analyze the graphs three ways.  First, graphs can be plotted and assessed visually.  All three graphs can be visualized; but visual assessment cannot be scaled up to thousands of replications or reproduced computationally.  

Second, slopes can be calculated using a simple univariate linear regression.  Both Young's p-value plot and \SchSp's p-value plot are intended to be analyzed using the slope.  Simonsohn et al.'s p-curve is *not* intended to be analyzed this way, so I do not analyze its slope.  Because Young's p-value plot is a rescaled QQ-plot, I do construct the QQ-plot and analyze its slope.  

Third, linearity can be tested by fitting two regression models — one linear, one quadratic — and selecting the model with the better fit.  I use two standard methods of model selection to compare regression models, AIC and an F-test *[cites]*.  AIC is a statistic of a regression model *[related to divergence]* that penalizes models with "too many" terms; in standard practice the model with the lower AIC is considered a better fit.  (Note that, while AIC is a random variable and hence has a standard error, this error is not taken into account in standard practice.)  The F-test tests the null hypothesis that the linear and quadratic regression models explain the same variance in the data.  I use the conventional 0.05 threshold for statistical significance in this test.  

\SchSp's p-value plot is not intended to be analyzed with respect to linearity; they assume that for "small" p-values the plot will not be linear *[cite]*.  Because Young's p-value plot is a rescaled QQ-plot, AIC and F-test assessments of linearity will be the same.  Because the QQ-plot's domain is stable (not dependent on $N$), I assess linearity using the QQ-plot rather than Young's p-value plot.  

The results of the slope and linearity analyses are random variables, determined entirely by the input data (and free parameters such as the threshold of statistical significance used in the F-test of linearity).  So the simulation can be used to evaluate the frequentist properties of Young's p-value plot — and the way these properties vary with the real effect — by examining their distributions across hundreds of runs of the simulation and the five different conditions.  


## Evidence ##

Again, Young and collaborators use Young's p-value plot to provide evidence supporting two skeptical hypotheses about air pollution hazards.  Using the simulation method, I formally evaluate the evidentiary value of this method, using two quantifiable conceptions of evidence that have been widely discussed in the philosophy of science literature.  These two conceptions of evidence are typically treated as philosophical rivals and associated with different philosophies of probability and statistics.  Using both allows me to remain agnostic to these disagreements.  

The severity conception of evidence is associated with philosopher of science Deborah Mayo's reconceptualization of frequentist hypothesis testing *[cites]*. [I note that Young offered a positive blurb of @MayoStatisticalInferenceSevere2018, stating that "Her severity requirement [sic] demands that the scientist provide a sharp question and related data. Absent that, the observer should withhold judgment or outright reject."]  In its most recent form [@MayoStatisticalInferenceSevere2018], the severity conception of evidence involves two "severity requirements":  

- *Severity (weak)*: One does not have evidence for a claim if nothing has been done to rule out ways the claim may be false. If data $x$ agree with a claim $C$ but the method used is practically guaranteed to find such agreement, and had little or no capability of finding flaws with $C$ even if they exist, then we have bad evidence, no test (BENT). [@MayoStatisticalInferenceSevere2018 5]
- *Severity (strong)*: We have evidence for a claim $C$ just to the extent it survives a stringent scrutiny. If $C$ passes a test that was highly capable of finding flaws or discrepancies from $C$, and yet none or few are found, then the passing result, $x$, is evidence for $C$. [@MayoStatisticalInferenceSevere2018 14]

On this conception of evidence, a test or analytical method $T$ with observed output $t$ can provide evidence supporting a hypothesis $H$ only if, counterfactually, $T$ would have given a different output if $H$ were false.  Hypothesis testing assesses this counterfactual in the form of a p-value using a mathematical model of the data-generating process: $p = pr(T = t | \lnot H) = L(\lnot H; T = t)$, where the role of $\lnot H$, "$H$ is not true," is usually played by a null hypothesis $H_0$.  A small p-value indicates that the counterfactual is probably true, that is, if $H$ were false then $T$ would probably have given a different output.  (This isn't sufficient for $T=t$ to provide evidence; other background conditions must apply, such as an absence of p-hacking and that the instruments used to produce the data were working correctly *[cite]*.)  On the other hand, a large p-value indicates that the test "is practically guaranteed" to produce this output, and so in this case by the weak severity principle "we have bad evidence, no test." 

The severity conception of evidence can be applied to Young and collaborator's skeptical hypotheses about air pollution hazards as follows.  The claims $H$ are the zero hypothesis $\delta = 0$ and mixture hypothesis $\delta = \{0, .6\}$.  The method or test $T$ is Young's p-value plot, analyzed visually or in terms of slopes and linearity; the outputs $t$ that I will consider are (i) the "hockey stick" pattern for the visual analysis; (ii) "gaps" in the visual analysis; (iii) slope of $1\pm0.1$ (on the QQ-plot), corresponding to the standard meaning of a "45-degree line"; and (iv) non-linearity inferences from either AIC or the F-test.  To generate p-values (or informally assess severity for the visual analysis), we need to specify a hypothesis for $\lnot H$.  I will use each of the following:  

a. $\delta = 0.2$
b. $\delta = 0.4$
c. $\delta = 0.6$
d. $a \vee b \vee c$ (that is, any one of the non-zero effects)
e. $\delta = \{0, .6\}$ for $H: \delta = 0$ and vice-versa (that is, the other skeptical hypothesis)
f. $d \vee e$ (the disjunction of any one of the non-zero effects and the other skeptical hypothesis)

In each case, insofar as the p-value is large $p > .05$, this implies that Young's p-value plot does not provide evidence to support the zero or mixture hypotheses.  

The likelihood conception of evidence is not strongly associated with any one statistician or philosopher of science, though it is associated with an approach to Bayesian statistics where the likelihood ratio is equivalent to the "Bayes factor" (if the prior probabilities are the same) *[cites]*.  

Formally, the likelihood conception of evidence compares two rival hypotheses $H_1$ and $H_2$ using some data $d$.  The likelihood ratio is defined as
\[ K(H_1, H_2; d) = \frac{L(H_1; d)}{L(H_2; d)} = \frac{pr(d | H_1)}{pr(d | H_2)}.]
If $K > 1$, then the evidence favors $H_1$; and $K < 1$ then the evidence favors $H_2$.  Sometimes $\log K$ is used to create symmetry between $H_1$ and $H_2$.  On one common interpretive scale, $\abs(\log_{10} K) < 0.5$ is "Not worth more than a bare mention" (minimal evidence); $0.5 < \abs(\log_{10} K) < 1$ is "substantial"; $1 < \abs(\log_{10} K) < 2$ is "strong"; and $2 < \abs(\log_{10} K)$ is "decisive".  *[https://www.jstor.org/stable/2291091?seq=1]*

To apply the likelihood conception of evidence to Young and collaborator's skeptical hypotheses about air pollution, $H_1$ will be the zero or mixture hypothesis, the rival hypothesis $H_2$ will be the hypotheses (a-f), and the data $d$ will be the analysis outputs (i-iv).  In each case, insofar as $K < 1$, this implies that Young's p-value plot does not provide evidence to support the zero or mixture hypotheses.  



# Results #

## Model validation ##

## Visual analyses ##

## Computationally reproducible analyses ##



# Discussion #

# Conclusion #

