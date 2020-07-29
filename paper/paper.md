---
title: "Young's p-value plot does not provide evidence against air pollution hazards"
author: "Daniel J. Hicks"
institute: "University of California, Merced"
email: hicks.daniel.j@gmail.com
abstract: "*[goes here]*"
bibliography: Young.bib

header-includes:
  - \usepackage[nolists]{endfloat}
  - \DeclareDelayedFloatFlavour*{longtable}{table}
---

\newcommand{\SchSp}{Schweder and Spjøtvoll}

<https://ehp.niehs.nih.gov/authors/research-articles>

# Introduction #

In a recent series of papers [@YoungCerealinducedGenderSelection2009; @YoungAssessingGeographicHeterogeneity2013; @YouPM2OzoneIndicators2018; @YoungReliabilityEnvironmentalEpidemiology2019; @YoungAmbientAirPollution2019; @YoungCombinedBackgroundInformation2019; @YoungEvaluationMetaanalysisAir2019; @YoungReExtendedMortality2020], statistician S. Stanley Young and collaborators have criticized epidemiological studies and meta-analyses, especially those that find hazards associated with particulate matter (PM) exposure.  Young and collaborators have used a graphical method they call a "p-value plot," claiming that this method reveals zero effects, heterogeneity, and p-hacking in the primary air pollution literature and so implies that meta-analyses are unreliable.  The current paper uses simulation methods to show that the p-value plot method does not have the statistical properties required to support the criticisms leveled by Young and collaborators.  Therefore Young's p-value plot does not provide evidence against hazards of air pollution.  

In their papers criticizing air pollution epidemiology, Young and collaborators have made three kinds of arguments.  The most common argument observes that analysts of complex data sets must have a number of decisions about which hypotheses to investigate and how to construct models, then relates the large number of possible combinations of decisions to the traditional frequentist problem of correcting for multiple comparisons [@YoungCerealinducedGenderSelection2009; @YoungDemingDataObservational2011; @YoungReModelingAssociation2013; @PeaceReliabilityNutritionalMetaanalysis2017; @YoungAirQualityEnvironmental2017; @YoungAirQualityAcute2017; @YouPM2OzoneIndicators2018; @YoungAirPollutionMortality2018; @YoungReliabilityEnvironmentalEpidemiology2019; @YoungEvaluationMetaanalysisAir2019; @YoungReExtendedMortality2020].  For a critical assessment of this argument, see *[me]*.  A second argument deploys a statistical method that Young and collaborators call "local control," which first clusters geographic units and then evaluates hazards separately within each cluster [@ObenchainAdvancingStatisticalThinking2013; @YoungBiasResponseHeterogeneity2015; @ObenchainLocalControlStrategy2017; @ObenchainLowlevelRadonExposure2019].  This method is not examined here.  

The third argument made by Young and collaborators examines a collection of p-values, either extracted from the primary studies in a meta-analysis or similar aggregate of studies [@YoungReliabilityEnvironmentalEpidemiology2019; @YoungAmbientAirPollution2019; @YoungCombinedBackgroundInformation2019; @YoungEvaluationMetaanalysisAir2019; @YoungReExtendedMortality2020] or from testing multiple distinct hypotheses on a single complex dataset [@YoungCerealinducedGenderSelection2009; @YoungAssessingGeographicHeterogeneity2013; @YouPM2OzoneIndicators2018]. The p-values are visualized using what Young and collaborators call a "p-value plot," and features of the plot are interpreted as indicating that there is no underlying effect [@YoungCerealinducedGenderSelection2009; @YouPM2OzoneIndicators2018; @YoungCombinedBackgroundInformation2019; @YoungReExtendedMortality2020], a heterogeneous mixture of zero and nonzero effects [@YoungReliabilityEnvironmentalEpidemiology2019; @YoungAmbientAirPollution2019; @YoungCombinedBackgroundInformation2019; @YoungEvaluationMetaanalysisAir2019], or that the authors of the primary studies have engaged in p-hacking [@YoungEvaluationMetaanalysisAir2019]. 

However, Young and collaborators have provided only a minimal analysis of the statistical properties of the p-value plot method, in a set of notes that has been made available on the arXiv, an open paper repository, but has not undergone peer review [@YoungCombinedBackgroundInformation2019].  They have attempted to justify their method by frequently citing plots of p-values developed by other authors [@SchwederPlotsPvaluesEvaluate1982; @SimonsohnPcurveKeyFiledrawer2014].  But these other plots are designed to answer different questions and have different visual and statistical properties.  

Despite this lack of justification, Young and collaborators have drawn strong conclusions using their p-value plot method, claiming that "causality of PM10/PM2.5 on heart attacks is not supported" [@YoungReliabilityEnvironmentalEpidemiology2019], and Young has advocated that "regulation of PM2.5 should be abandoned altogether" [@YoungSuggestionsEPA2017].  Perhaps because most of the papers using the p-value plot method were published in 2019, they have not been highly cited.  But a slightly older paper by Young and colleagues using the first, potential comparisons argument [@YoungAirQualityAcute2017] was cited in the scientific review of US EPA's Ozone Integrated Science Assessment [@USEPACleanAirScientificAdvisoryCommitteeReviewEPAIntegrated2020, p.D-50].  

Therefore, the motivation for the current paper is to formally evaluate the statistical properties and evidentiary value of the p-value plot method used by Young and collaborators, before it can influence regulatory decisionmaking.  

For this purpose, simulations offer a number of advantages over exact methods.  First, for the target audience of this paper, well-written R code is likely to be more accessible than an exact analysis of the statistical properties of meta-analyses.  This makes it easier for readers to scrutinize the details of my analysis.  Second, the simulation is modular, making it more extensible or modifiable by readers than an exact analysis.  The simulated data generating process and its analysis can be modified by simply rewriting two short functions and re-running the simulation.  An exact analysis of a different data-generating process might require a completely different formal approach.  Third, it may be algebraically intractable to give exact characterizations of constructions such as plots except in the limit (as sample size and/or number of studies goes to infinity).  


# Methods #

I first discuss Young's p-value plot, the graphical method used as the basis of the third argument by Young and collaborators against air pollution epidemiology.  I then turn to the simulation methods used in the current study to analyze the statistical properties of Young's p-value plot.  

With one exception, every peer-reviewed publication describing the method has been written by a collection of authors.  In each of these works the method has been referred to as a "p-value plot," without attributing it to any individual.  However, Young is the only author who appears across all of these papers, and the method needs to be distinguished from other, nominally similar methods.  I therefore refer to the method as "Young's p-value plot."  

## Young's p-value plot ##

Young and collaborators have not provided a formal definition of their p-value plot method, and at least one key term in the methods appears to be used in a non-standard way.  They have attempted to justify their method by associating it with two other graphical methods, but these methods are not equivalent to each other and do not appear to be equivalent to Young's p-value plot.  

\renewcommand{\P}{\ensuremath{\mathbb{P}}}

All three graphical methods take as input a set of $N$ p-values $\P = \{p_1, p_2, \ldots, p_N\}$, nominally produced by applying a given hypothesis statistical test to $N$ replications of a given study design, each replication drawing samples of size $n$ from a given population.  This corresponds to the simplest case of meta-analysis.  Thus the p-values in $\P$ are nominally samples from a single underlying distribution $p_i \sim P$.  Note that, if the real effect is zero $\delta = 0$, then $P$ is the uniform distribution on $[0,1]$.  

Young and collaborators have also applied Young's p-value plot to collections of p-values produced by, for example, conducting hundreds of possible analyses on a single dataset [@YoungCerealinducedGenderSelection2009] or all partial correlations between pairs of covariates [@YoungAssessingGeographicHeterogeneity2013].  This method resembles "multiverse analysis," which systematically examines how methodological decisions shape the outputs of statistical analysis [@SteegenIncreasingTransparencyMultiverse2016].  However, the point of multiverse analysis is to identify patterns in the associations between decisions and outputs.  Young and collaborators have drawn much less nuanced conclusions, such as that "all of the p-values are simply the result of chance" [@YoungCerealinducedGenderSelection2009].  Because simulating these applications would require modeling complex multivariate distribution, I do not consider these applications here.  

### Schweder and Spjøtvoll's p-value plot ###

Young and collaborators have frequently cited the "p-value plot" presented in @SchwederPlotsPvaluesEvaluate1982.  For this p-value plot, let $rank_{desc}(p_i)$ be the (1-indexed) *descending rank* of $p_i \in \P$, i.e., $rank_{desc}(p_i)$ is the number of p-values $p_j \in \P$ greater than or equal to $p_i$.  The largest p-value has descending rank 1, and the smallest p-value has descending rank $N$.  Then \SchSp's p-value plot plots the graph $(1-p_i, rank_{desc}(p_i))$.  (Examples of \SchSp's p-value plot are included in the supplemental materials.)

@SchwederPlotsPvaluesEvaluate1982 briefly argue that the relationship between $1-p_i$ and $rank_{desc}(p_i)$ should be approximately linear "when [$p_i$] is not too small" and that, from left to right, "often, the plot will not show a clearcut break but rather a gradual bend" away from linearity.  Note that this means, as a method, \SchSp's p-value plot generally ignores small p-values, which roughly correspond to the statistically significant p-values.  They illustrate their method with example data and a straight line "drawn by visual fit" rather than regression analysis.    

### Simonsohn, Nelson, and Simmons' p-curve ###

Young and collaborators have also echoed concerns about the replication crisis unfolding in social psychology and certain areas of biomedical research *[cites]*.  In particular, they appeal to concerns about "p-hacking," in which researchers conduct variant analyses (e.g., including or excluding various controls in a regression analysis) and report only the analyses with statistically significant p-values [@SimonsohnPcurveKeyFiledrawer2014].  Note that, other than Young's p-value plot, Young and collaborators provide no empirical evidence of p-hacking in environmental epidemiology, and to my knowledge no such evidence has been published.  *[general study]*

Young and collaborators have attempted to associate Young's p-value plot with a method developed to detect p-hacking, called a "p-curve" [@SimonsohnPcurveKeyFiledrawer2014].  The intuition behind the p-curve is that p-hacking will tend to produce an excess number of p-values "just below" the conventional $0.05$ threshold for statistical significance.  Formally, Simonsohn et al.'s p-curve first divides the interval $[0, 0.05]$ into 5 bins at the thresholds $0.01, 0.02, 0.03, 0.04, 0.05$, then calculates $N_b$, the number of p-values in bin $b$.  The p-curve is the graph $(b_t, N_b)$, where $b_t$ is the threshold for bin $b$.  The method then formally tests for p-hacking by applying statistical tests of the null hypothesis that the restricted distribution $P|_{p < 0.05}$ is uniform.  (Examples of Simonsohn et al.'s p-curve are included in the supplemental materials.)

Note that the p-curve is equivalent to a histogram on the interval $[0, 0.5]$ with binwidth $0.01$ and the p-curve only includes statistically significant p-values.  Thus the p-curve and \SchSp's p-value plot not only produce different kinds of plots using the same input data $\P$, but actually direct their attention to different — typically disjoint — subsets of $\P$.  The two kinds of plots cannot be equivalent to each other.  At no point have Young and collaborators acknowledged this fundamental difference between the two methods that they cite as support for their own method.  

### Young's p-value plot ###

For Young's p-value plot, let $rank_{asc}(p_i)$ be the (1-indexed) *ascending rank* of $p_i \in \P$, i.e., $rank_{asc}(p_i)$ is the number of p-values $p_j \in \P$ less than or equal to $p_i$.  The smallest p-value has ascending rank 1, and the largest p-value has ascending rank $N$.  Without loss of generality, if \P is already in ascending order $p_1 < p_2 < \cdots < p_N$, then $rank_{asc}(p_i) = i$.  And Young's p-value plot is the graph $(i, p_i)$.  

Note that Young's p-value plot is a rescaled QQ-plot of $\P$ against the uniform distribution, with the theoretical quantiles $q_i = \frac{i}{N} = \frac{rank_{asc}(p_i)}{N}$.  In addition, $rank_{desc}(p_i) = N - rank_{asc}(p_i) + 1$, and so given the number of studies $N$ Young's p-value does have a 1-1 mathematical relationship to \SchSp's p-value plot.  However, because the axes are switched, a regression line fit to \SchSp's p-value plot will not necessarily correspond to a regression line fit to Young's p-value plot.  Thus there is not a 1-1 relationship between the slopes of the two plots.  

Young and collaborators have given the following method to analyze Young's p-value plot:  "Evaluation of a p-value plot follows a logical path. Is the p-value plot homogeneous? If the points roughly [sic] on a 45-degree, they are homogeneous and consistent with randomness; a lessor slope with all points roughly on a line indicates a consistent effect even if some of the individual p-values are not considered statistically significant. If the effects differ, one from another, beyond chance, then the effects are heterogeneous" [@YoungReliabilityEnvironmentalEpidemiology2019 48].  

A "45-degree line" typically refers to the graph of an identity function or the line $y = x$, which forms a 45-degree angle with respect to both the x- and y-axes.  This interpretation makes sense for a QQ-plot of p-values against the uniform distribution on $[0, 1]$, where both axes are on the scale $[0, 1]$; here a slope of 1 indicates that the underlying distribution $P$ is uniformly distributed, which in turn indicates that the real effect is zero.  Strictly speaking, this interpretation does not make sense for Young's p-value plot, where the x-axis (rank) is on the scale $[1, N]$ and the y-axis (p-value) is on the scale $[0, 1]$, giving a maximum possible slope of $1/N$.  I will assume that a "45-degree line" generally means that the slope of the equivalent QQ-plot is 1.  Note that Young and collaborators only evaluate slopes visually; they do not fit regression models or use any other quantitative methods to measure slopes.  

Young and collaborators have frequently identified a "hockey stick" pattern [@YoungReliabilityEnvironmentalEpidemiology2019; @YoungAmbientAirPollution2019; @YoungCombinedBackgroundInformation2019; @YoungEvaluationMetaanalysisAir2019], "small p-values to the lower left ... and then points ascending in a roughly 45-degree line," which here seems to mean that the pattern is linear rather than having a particular slope [@YoungEvaluationMetaanalysisAir2019 5].  Young and collaborators have taken the "hockey stick"  to indicate some combination of heterogeneous effects, p-hacking, and researcher misconduct.  On two occasions Young and collaborators have formally tested non-linearity by comparing a linear and quadratic regression using an F-test [@YoungReliabilityEnvironmentalEpidemiology2019 50; @YoungCombinedBackgroundInformation2019 18].  More often they appear to have identified the "hockey stick" pattern by visual inspection alone.  

Young and collaborators do not appear to have validated the p-value plot method in any peer-reviewed studies.  @YoungCombinedBackgroundInformation2019 (8-9) briefly reported a simulation study:  "We simulate 15 p-values ten times and produce the p-value plots, Figure 1. Each of the 10 simulated meta-analysis studies has 15 base papers .... The goal is to get an impression of the variability of a p-value plot under ideal conditions where there is nothing going on and there is no multiple testing or multiple modeling in the base studies."  In my notation, each run of the simulation had $N = 15$ p-values drawn from a uniform distribution, i.e., the zero hypothesis is true.  They conducted $NN = 10$ simulation runs, but do not report any aggregate or distributional statistics; instead they report the smallest (rank = 1) p-value and the p-value plot for each simulation run.  Note that they did not report any simulations of cases where the real effect was non-zero or heterogeneous.  They also do not report making the simulation code available anywhere for checking reproducibility or extending their analysis.  


## The current simulation study ##

The current study is an automatically-reproduced simulation study of Young's p-value plot, systematically testing it across five kinds of cases, using formal analyses of the graphs, and comparing it to Schweder and Spjøtvoll's p-value plot and Simonsohn, Nelson, and Simmons' p-curve.  

Each run of the simulation is composed of $N$ studies.  To make the study parameters easy to understand, each study is based on a two-sample t-test.  Two samples, each of size $n$, are drawn from Gaussian distributions with mean 0 and $\delta$, respectively, and common standard deviation $\sigma$.  In the basic "single process" case, $\delta$ is the same for all studies in the simulation run, and so $\delta$ corresponds to the true effect size of the (single) process being studied.  

A *condition* of the simulation fixes the values of these parameters.  Conditions can be systematically varied to compare, e.g., different effect sizes, and multiple runs $NN$ in each condition allow us to analyze the frequentist properties of Young's p-value plot within and across conditions.  For the primary analysis of the current paper, 5 different effect sizes are used — corresponding to zero, small, moderate, moderate-large effects, and a "mixed" or heterogeneous condition in which half the population has a zero response and half the population has a moderate-large response — and each condition is simulated with 500 runs.  500 runs per condition was chosen because this should produce good estimates of central tendency and variation with just a few minutes of computation time.  Note that 500 runs is probably too few to produce good estimates of extreme values, e.g., 5th percentiles.  The primary conditions examined in the current study are summarized in [@tbl:params].  

| parameter | meaning                             | value(s) |
|:----------|:------------------------------------|---------:|
| $\delta$  | real effect size                    | 0, .2, .4, .6, (0, .6) |
| $\sigma$  | s.d. of samples                     | 1        |
| $n$       | study sample size                   | 60       |
| $N$       | num. studies in each run            | 20       |
| $NN$      | num. runs for each condition        | 500      |

Table: Base parameter values used in the current simulation study.  The real effect size of (0, .6) indicates the "mixed" or heterogeneous condition, in which half the population have no response and half the population have a strong response. {#tbl:params}

The simulation can also simulate a "mixed," multi-process, or heterogeneous condition, in which subsets of the population have different responses to the intervention or exposure.  The current study uses a (0, .6) mixture, meaning that half of the population has a zero effect or no response and half the population has a medium-strong response of .6.  These values were chosen to create a relatively dramatic or striking mixture that is less likely to look like a single-process or homogeneous case.  For this mixed case, the simulation draws $n$ control/non-exposed group samples (as with the single process case) and then, for the intervention/exposed group, draws $n/2$ samples from the no response subpopulation and $n/2$ samples from medium-strong response subpopulation.  This guarantees that the sample is balanced across the two (latent) subpopulations.  I expect that unbalanced samples will tend to produce less striking p-curves, as unbalanced samples will tend to look more like the majority subpopulation.  (This approach to modeling a "mixed" condition can be considered a constrained and simplified version of a Gaussian mixture *[cite]*.  So a generalization of the simulation used here might start instead by drawing from a standard Gaussian mixture, and treating single-process conditions as limiting cases with a single mixture component.)   

After generating the data, it is straightforward to construct Young's p-value plot.  I then analyze the plot three ways.  First, graphs can be plotted and assessed visually.  All three graphs can be visualized; but visual assessment cannot be scaled up to thousands of replications or reproduced computationally.  

Second, slopes can be calculated using a simple univariate linear regression.  If the slope of the QQ-plot is approximately 1, this indicates that the underlying distribution of p-values $P$ is uniform, which in turn indicates that the zero hypothesis is true.  Because the slope of Young's p-value plot is $N$ times the slope of the QQ-plot, analysis of the QQ-plot is simpler than analyzing Young's p-value plot directly.  

After fitting a linear regression to the QQ-plot, I formally assess its slope in two ways, first by considering whether it is in the range $1\pm0.1$ and second by using to Z-test to evaluate whether it is statistically significantly different from 0.  In addition, I use a Kolmogorov-Smirnov test (KS-test) to compare the observed set of p-values \P to the uniform distribution.  I use the conventional 0.05 threshold for statistical significance for both the Z-test and KS-test.  

The third analytical approach is to evaluate the linearity of the plot.  Young and collaborators take non-linearity as evidence of some combination of heterogeneity, p-hacking, and researcher misconduct.  Linearity of a plot can be tested by fitting two regression models to its graphs — one linear, one quadratic — and selecting the model with the better fit.  I use two standard methods of model selection to compare regression models, AIC and an F-test.  I use the conventional 0.05 threshold for statistical significance in the F-test.  As with slope, for simplicity I analyze the QQ-plot rather than Young's p-value plot directly.  

The results of the slope and linearity analyses are random variables, determined entirely by the input data (and free parameters such as the threshold of statistical significance used in the F-test of linearity).  So the simulation can be used to evaluate the frequentist properties of Young's p-value plot — and the way these properties vary with the real effect — by examining their distributions across hundreds of runs of the simulation and the five different conditions.  


## Evidence ##

Again, Young and collaborators use Young's p-value plot to provide evidence supporting two skeptical hypotheses about air pollution hazards.  Using the simulation method, I formally evaluate the evidentiary value of this method, using two quantifiable conceptions of evidence that have been widely discussed in the philosophy of science literature.  These two conceptions of evidence are typically treated as philosophical rivals and associated with different philosophies of probability and statistics.  Using both allows me to remain agnostic to these disagreements.  

### Severity ###

The severity conception of evidence is associated with philosopher of science Deborah Mayo's reconceptualization of frequentist hypothesis testing [@MayoErrorGrowthExperimental1996; @MayoStatisticalInferenceSevere2018]. [Note that Young offered a positive blurb of @MayoStatisticalInferenceSevere2018, stating that "Her severity requirement [sic] demands that the scientist provide a sharp question and related data. Absent that, the observer should withhold judgment or outright reject."]  In its most recent form, the severity conception of evidence involves two "severity requirements":  

- *Severity (weak)*: One does not have evidence for a claim if nothing has been done to rule out ways the claim may be false. If data $x$ agree with a claim $C$ but the method used is practically guaranteed to find such agreement, and had little or no capability of finding flaws with $C$ even if they exist, then we have bad evidence, no test (BENT). [@MayoStatisticalInferenceSevere2018 5]
- *Severity (strong)*: We have evidence for a claim $C$ just to the extent it survives a stringent scrutiny. If $C$ passes a test that was highly capable of finding flaws or discrepancies from $C$, and yet none or few are found, then the passing result, $x$, is evidence for $C$. [@MayoStatisticalInferenceSevere2018 14]

On this conception of evidence, a test or analytical method $T$ with observed output $t$ can provide evidence supporting a hypothesis $H$ only if, counterfactually, $T$ would have given a different output if $H$ were false.  Hypothesis testing assesses this counterfactual in the form of a p-value using a mathematical model of the data-generating process: $p = pr(T = t | \lnot H) = L(\lnot H; T = t)$, where the role of $\lnot H$ (read "$H$ is false") is played by a null hypothesis $H_0$, which often states that some quantity of interest is zero.  A small p-value indicates that the counterfactual is probably true, that is, if $H$ were false then $T$ would probably have given a different output.  (This isn't sufficient for $T=t$ to provide evidence; other background conditions must apply, such as an absence of p-hacking and that the instruments used to produce the data were working correctly *[cite]*.)  On the other hand, a large p-value indicates that the test "is practically guaranteed" to produce this output, and so in this case by the weak severity principle "we have bad evidence, no test." 

(Note that, throughout this paper, I am distinguishing a zero hypothesis — that some real effect is zero — from a null hypothesis — the alternate or rival hypothesis used to calculate a p-value.  Null hypotheses can claim that some quantity of interest is equal to a non-zero value.  For example, if we are using a QQ-plot to test whether a set of values is uniformly distributed, the null hypothesis is that the slope is 1.)

Severity can also be evaluated informally when a p-value cannot be calculated.  Consider visual features of plots, such as the "hockey stick" pattern that Young and collaborators interpret as evidence of a mixed or heterogeneous process.  If this pattern is common in cases where there is only a single processes, then again the weak severity principle implies that this visual pattern does not provide evidence of a mixed process.  

The severity conception of evidence can be applied to Young and collaborator's skeptical hypotheses about air pollution hazards as follows.  The claims $H$ are the zero hypothesis $\delta = 0$ and mixture hypothesis $\delta = \{0, .6\}$.  The method or test $T$ is Young's p-value plot, analyzed visually or in terms of slopes and linearity; the outputs $t$ that I will consider are (i) the "hockey stick" pattern for the visual analysis; (ii) "gaps" in the visual analysis; (iii) slope of approximately 1 (on the QQ-plot), corresponding to the standard meaning of a "45-degree line," assessed based on whether the slope (as estimated using linear regression) falls into the range $1\pm0.1$, non-statistically significant results for a Z-test against the null hypothesis that the slope is exactly 1, and non-statistically significant results for the KS-test of uniformity; and (iv) non-linearity inferences using AIC and the F-test.  

[@Tbl:outputs] summarize the outputs examined in the current study.  

|     | output            | determined using      | taken as evidence for |
|----:|:------------------|:----------------------|:----------------------|
|   i | "hockey stick"    | visual inspection     | mixed effect          |
|  ii | "gaps"            | visual inspection     | p-hacking or other problems |
| iii | slope $\approx$ 1 | range $1\pm0.1$       | zero effect           |
|     |                   | Z-test not stat. sig. ||
|     |                   | KS-test not stat. sig.||
|  iv | non-linearity     | AIC: quadratic        | mixed effect
|     |                   | F-test: stat. sig.    ||

Table: Outputs of Young's p-value plot examined using the simulation.  "Outputs" are features of plots that Young and collaborators take as evidence for critical assessments of air pollution epidemiological studies.  The "determined using" column indicates how these outputs are identified as present/absent in the current study.  {#tbl:outputs}

To generate p-values (or informally assess severity for the visual analysis), we need to specify a null hypothesis for $\lnot H$.  I will use each of the following:  

a. weak effect: $\delta = 0.2$
b. moderate effect: $\delta = 0.4$
c. strong-moderate effect: $\delta = 0.6$
d. non-zero effect: $a \vee b \vee c$ (read $\vee$ as "or")
e. mixed effect: $\delta = \{0, .6\}$ for $H: \delta = 0$ and vice-versa (that is, the other skeptical hypothesis)
f. any other effect: $d \vee e$ (any of the non-zero effects or the other skeptical hypothesis)

In each case, insofar as the p-value is large $p > .05$, this implies that Young's p-value plot does not provide evidence to support the zero or mixture hypotheses.  

### Likelihood ###

The likelihood conception of evidence is not strongly associated with any one statistician or philosopher of science, though it is associated with an approach to Bayesian statistics where the likelihood ratio is equivalent to the "Bayes factor" (if the prior probabilities are the same) *[cites]*.  

Formally, the likelihood conception of evidence compares two rival hypotheses $H_1$ and $H_2$ using some data $d$.  The likelihood ratio is defined as
$$ K(H_1, H_2; d) = \frac{L(H_1; d)}{L(H_2; d)} = \frac{pr(d | H_1)}{pr(d | H_2)}. $$
If $K > 1$, then the evidence favors $H_1$; and $K < 1$ then the evidence favors $H_2$.  Sometimes $\log K$ is used to create symmetry between $H_1$ and $H_2$.  On one common interpretive scale, $\left|\log_{10} K\right| < 0.5$ is "not worth more than a bare mention" (not supporting either hypothesis); $0.5 < \left|\log_{10} K\right| < 1$ is "substantial"; $1 < \left|\log_{10} K\right| < 2$ is "strong"; and $2 < \left|\log_{10} K\right|$ is "decisive".  *[https://www.jstor.org/stable/2291091?seq=1]*

To apply the likelihood conception of evidence to Young and collaborator's skeptical hypotheses about air pollution, $H_1$ will be the zero or mixture hypothesis, the rival hypothesis $H_2$ will be the hypotheses (a-f), and the data $d$ will be the analysis outputs (i-iv).  (For simplicity, the same dichotomous frequentist test outputs are used, e.g., statistically significant or not, rather than continuous-valued likelihoodist or Bayesian alternatives.)  In each case, insofar as $K < 0.5$, this implies that Young's p-value plot does not provide evidence to support the zero or mixture hypotheses.  

## Reproducibility ##

The simulation, analysis, and outputs (figures and tables included here) are publicly available and automatically reproduced.  Code is available at <https://github.com/dhicks/p_curve> *[version]* and the automatically reproduced analysis can be viewed at <https://dhicks.github.io/p_curve/>.  

The simulation and analysis are both written in R *[cite]* and make extensive use of the `tidyverse` suite of packages for R *[cite]* *[versions]*



# Results #

Figures and tables included in the automatically reproduced analysis validate the simulation, showing that both the individual studies and meta-analyses (average within each simulation run) produce unbiased estimates of the real effect size $\delta$.  The only exception is the mixed case, where half of the population experiences zero effect ($\delta = 0$) and half the population experiences a medium-strong effect ($\delta = 0.6$).  In this case the average effect is 0.3, and on average the simulation gives this estimate.  


## Visual analyses ##

[@Fig:samples_young] shows 10 examples of Young's p-value plot across the five conditions.  See the supplemental material for examples of \SchSp's and Simonsohn et al.'s plots.  

![Examples of Young's p-value plot, drawn at random from the simulation results.  Rows and colors correspond to conditions or real effects ($\delta$), from zero (0) to moderate-strong (0.6) and a mixed condition $\delta = \{0.0, 0.6\}$.  Columns correspond to indices for the simulation runs that produced these results, and are not meaningful.  (In particular, there is no relationship between simulation run $j$ in condition $a$ and simulation run $j$ in condition $b$.)  In Young's p-value plot, each point corresponds to a single p-value in the meta-analysis (simulation run); the x-axis is the ascending rank of the p-value in the set \P, and the y-axis is the p-value itself.](../out/samples_young.png){#fig:samples_young}

As preliminary observations, note that there are both substantial qualitative differences within effect sizes (within rows) as well as substantial qualitative similarity across consecutive effect sizes (comparing adjacent rows).  Larger effect sizes have more statistically significant results, resulting in a series of small p-values that gradually bend up.  Even small and zero effects can look nonlinear (\# 97 for $\delta = 0$; \# 90 for $\delta = 0.2$).  The mixed condition plots (bottom row) are not obviously distinct from the small $\delta = 0.2$ and moderate effects $\delta = 0.4$ plots (second and third rows).  Across all effect sizes, plots tend to have both statistically significant and insignificant p-values.  

I focus on two visual patterns that Young and collaborators frequently discuss in their critiques of air pollution epidemiology:  (i) the "hockey stick" pattern; (ii) "gaps" in the plot.  

Young and collaborators take the "hockey stick" pattern to be evidence of a mixed-process or heterogeneous case.  The "hockey stick" comprises a more-or-less flat series of small p-values on the left (supposedly corresponding to the mixture component with a real effect) and a second series of steeply increasing p-values on the right (supposedly corresponding to the mixture component with zero effect).  This pattern is clearly visible in all of the example plots for the moderate, strong-moderate, and mixed effects (bottom three rows) and arguably is also present in most of the plots for the small effect ($\delta = 0.2$, second row).  Further, the "hockey stick" pattern in the mixed condition (bottom row) seems to be indistinguishable from the pattern in the moderate condition.  Because the "hockey stick" pattern appears in plots where there is only a single process, the weak severity criterion implies that *the hockey stick pattern does not provide evidence of a mixed-process or heterogeneous case*.  

Young and collaborators take visual "gaps" in the plot to be evidence of p-hacking, publication bias, or other questionable researcher practices *[cites]*.  These gaps are quite common in the example plots produced by the simulation, across all conditions.  Note that the simulation does not include p-hacking, publication bias, or other questionable researcher practices.  Thus the weak severity criterion implies that *gaps in the plot do not provide evidence of p-hacking, publication bias, or other questionable researcher practices*.  



## Computationally reproducible analyses ##

Two outputs of the analysis of Young's p-value plot can be quantified and reproduced computationally:  (iii) a slope of approximately 1 (corresponding to the "45-degree line" as it is typically understood) and (iv) non-linearity.  Young and collaborators present compatibility with a "45-degree line" as evidence for the zero effect hypothesis *[cites]* and non-linearity as evidence for the mixed-effect hypothesis *[cites]*.  
I use the QQ-plot corresponding to Young's p-value plot, assessing a slope of 1 in two ways and non-linearity three ways.  

### Severity analysis

@Fig:evidence_severity show the results of the severity analysis as p-values; see the supplemental materials for a table version of these results.  By the weak severity criterion, when the results of the severity analysis are $> .05$, the test output does not provide evidence against the null hypothesis and supporting the target hypothesis.  

![Results of the severity analysis for outputs (iii) slope of 1 and (iv) non-linearity.  Severity analysis results are reported as p-values:  small values (conventionally $< .05$) indicate a severe test with respect to the null hypothesis.  Panels correspond to null hypotheses, and y-axis values correspond to the severity assessment (as a p-value) for the output with respect to the given null hypothesis.  The gray region indicates conventionally small values, so outputs provide evidence only if their points are in the gray region.  Only a few of these outputs provide evidence to support the claims made by Young and collaborators.](../out/evidence_severity.png){#fig:evidence_severity}

The zero effect hypothesis is supposedly supported by a slope of approximately 1 on the QQ-plot.  The range method (whether the slope is within the range $1\pm0.1$, output iii-range) and the KS-test (output iii-KS) are severe against the null hypotheses of moderate and moderate-strong effect ($\delta = 0.4$ and $\delta = 0.6$); the Z-test method (output iii-Z) is only severe against moderate-strong effects ($\delta = 0.6$).  That is, the "45 degree line" test is capable of distinguishing zero effects from moderate and stronger effects.  But it cannot distinguish zero effects from weak effects ($\delta = 0.2$), mixed effects, or any condition that includes these two. 

The mixed-effect hypothesis, or heterogeneity, is supposedly supported by non-linearity.  Two quantified versions of this output are examined here, comparing linear and quadratic regressions using AIC (iv-AIC) and an F-test (iv-F).  The AIC and F-test evaluations of non-linearity do not provide severe tests against any of the alternative conditions.  That is, neither of the tests of linearity considered here provide evidence for heterogeneity.  



### Likelihood analysis

@Fig:evidence_likelihood show the results of the likelihood analysis; see the supplemental materials for a table version of these results.  Log likelihood ratios are reported, so results above $0.5$ support $H_1$ and results below $-0.5$ support $H_2$.  So, by the weak severity criterion, when the results of the likelihood analysis are $< 0.5$, the test output does not provide evidence supporting the target hypothesis of zero or mixed effect.  

![Results of the likelihood analysis for outputs (iii) slope of approximately 1 and (iv) non-linearity.  Each point gives the log likelihood ratio for a target vs. rival hypothesis, given an output.  Panels correspond to target hypothesis $H_1$, either a mixed effect $\delta = \{0, 0.6\}$ or a zero effect $\delta = 0$; point color corresponds to rival hypothesis $H_2$.  Position on the x-axis indicates the output; points have been jittered horizontally to reduce overlap. Position on the y-axis indicates the strength of the evidence that the output provides to the hypotheses:  values greater than 0.5 support $H_1$ over $H_2$ and values less than -0.5 support $H_2$ over $H_1$.  (Points at the plot margins have infinite value due to division by zero.)  Shaded regions indicate support for either hypothesis, in order from lightest to darkest:  none, "substantial," "strong," "decisive."  An interactive version of this plot is included in the automatic reproduction of the analysis for this paper.](../out/evidence_likelihood.png){#fig:evidence_likelihood}
 
The zero effect hypothesis (left panel) is supposedly supported by a slope of approximately 1.  All three methods provide decisive support for zero effect against moderate-strong effects, and substantive or better support against moderate effects.  The range test just barely provides substantive support against the disjunction of non-zero effects ($\delta = 0.2, 0.4, 0.6$), and the KS-test also provides substantive support against mixed effects and non-zero effects (either mixed or one of the non-zero effects). 

For the mixed effect hypothesis, all of the points for AIC and the F-test are in the "not worth mentioning" range, and so neither method provides evidence to support heterogeneity.  



# Discussion #

This simulation analysis finds that Young's p-value plot does not provide evidence for heterogeneity or p-hacking based on the "hockey stick" shape, "gaps" in the plot, or AIC or F-tests of non-linearity.  The method can provide evidence for zero effects based on a slope of 1 when compared to moderate or stronger effect; and depending on the conception of evidence and exact method used, can provide evidence for a zero effect against weak and mixed effects.  However, the evidence against weak and mixed effects is at best "substantial," not "strong."  

In the meta-analyses criticized by Young and collaborators, the estimated short-term effects of air pollution are small on a relative risk scale.  @NawrotPublicHealthImportance2011 estimated the effect for air pollution (increase of 10 $\mu g/m^3$ PM_{10}) on non-fatal myocardial infarction at 1.02 (95% CI 1.01–1.02); the point estimates for six pollutants reported by @MustaficMainAirPollutants2012 (increase of 10 $\mu g/m^3$ for all except carbon monoxide) were all in the range 1.003-1.048; and @LiuAmbientParticulateAir2019 estimated effects for PM_{10} (increase of 10 $\mu g/m^3$) on all-cause mortality of 1.044 (95% CI 1.039-1.050).  A precise conversion from risk ratios to Cohen's $d$ (equivalent to $\delta$ used in the simulations) is beyond the scope of this paper; however, using the rule of thumb that the risk ratio is approximately equal to the log odds when the outcome is rare (here, non-fatal myocardial infarction in the general public) *[https://sma.org/southern-medical-journal/article/odds-ratios-and-risk-ratios-whats-the-difference-and-why-does-it-matter/]* and the conversion factor $\sqrt{3}/\pi$ between the log odds and Cohen's $d$, a risk ratio of 1.05 is roughly equivalent to $d=0.03$, which in turn is roughly an order of magnitude smaller than the weak effect condition examined in the current study.  So, even if Young and collaborators had applied the KS-test to the p-values in these studies and the test had not been statistically significant, this would not necessarily provide even "substantial" evidence of a zero effect of PM on human health.  

Further, Young and collaborators did not claim zero effects in their criticisms of these air pollution epidemiological studies *[cite]*.  Instead they claimed to find heterogeneity, p-hacking, and publication bias.  The simulation results indicate that Young's p-value plot is incapable of providing evidence for any of these hypotheses.  The features that Young and collaborators point to — the "hockey stick" shape, "gaps," non-linearity — are readily produced by moderate and stronger effects, and can even appear in zero and weak-effect conditions.  Young and collaborators have drawn incorrect conclusions, which could have been avoided easily had they validated their methods.  

In this light, it is notable that Young and collaborators have been able to publish multiple papers in peer-reviewed journals using an unvalidated and ultimately flawed method.  Young and collaborators regularly raise concerns about the inability of peer review to catch flawed scientific research *[cites]*, and the evident purpose of Young's work in this area has been to demonstrate these failures of peer review.  Ironically, Young's p-value plot provides exactly this kind of demonstration.  



## Limitations ##
 
The most important limitation of the current study is that it deploys a simplified model of scientific data collection and analysis.  Actual epidemiological studies are often much more complex than a two-sample t-test with no confounders, temporal or spatial patterns, or sampling issues.  I noted above that a Gaussian process might be a more adequate way of representing heterogeneous effects.  However, the modular design of the simulation means that future research can extend the simulation to deploy more complex models.  In addition, the fact that Young's p-value plot cannot support the kinds of claims made by Young and collaborators for this oversimplified model strongly suggests that it also cannot support these kinds of claims for more complex real-world cases.  

Another limitation is that I explore only a limited set of meta-analysis parameters.  In particular, the primary analysis does not vary the number of studies $N$, and does not thoroughly explore the apparent interaction between sample size $n$, effect $\delta$, and power of the study design.  It also does not explore variation in study design, e.g., different sample sizes across studies.  

Similarly, I explore only a few ways of translating the visual patterns identified by Young and collaborators into reproducible, quantified outputs.  For example, other tests against the uniform distribution could be used to assess non-linearity.  However, I have attempted to relate the tests used here to standard statistical practice.  

Finally, the simulation does not include any representation of publication bias or p-hacking.  The simulation outputs can be used to support a severity analysis of these hypotheses — by calculating p-values for null hypotheses of no publication bias and no p-hacking — they cannot be used in the same way to support a likelihood analysis.  


# Acknowledgments #

No particular funding was allocated for this project.  

# Data Sharing #

All code to conduct the simulation and analyze the results is available at <https://github.com/dhicks/p_curve>.  The simulation and analysis are automatically reproduced, and the resulting report is available at <https://dhicks.github.io/p_curve/>. 
