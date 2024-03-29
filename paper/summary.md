---
title: "The p-value plot does not provide evidence against air pollution hazards"
author: "Daniel J. Hicks"
<!--abstract: >
  **Background**:  A number of papers by Young and collaborators have criticized epidemiological studies and meta-analyses of air pollution hazards using a graphical method that the authors call a p-value plot, claiming to find zero effects, heterogeneity, and p-hacking.  However, the p-value plot method has not been validated in a peer-reviewed publication. The aim of this study was to investigate the statistical and evidentiary properties of this method. 
  **Methods**: A simulation was developed to create studies and meta-analyses with known real effects $\delta$, integrating two quantifiable conceptions of evidence from the philosophy of science literature. The simulation and analysis is publicly available and automatically reproduced.  
  **Results**: In this simulation, the plot did not provide evidence for heterogeneity or p-hacking with respect to any condition. Under the right conditions, the plot can provide evidence of zero effects; but these conditions are not satisfied in any actual use by Young and collaborators. 
  **Conclusions**: The p-value plot does not provide evidence to support the skeptical claims about air pollution hazards made by Young and collaborators.-->
bibliography: Young.bib
csl: epidemiology.csl
notes-after-punctuation: true
header-includes:
  - \usepackage{booktabs}
  - \usepackage{longtable}
  - \usepackage{setspace}
  - \doublespacing
tblPrefix:
  - "table"
  - "tables"
---

\raggedright
\setlength{\parskip}{1em}
\renewcommand{\P}{\ensuremath{\mathbb{P}}}

| Type: Review
| Title: The p-value plot does not provide evidence against air pollution hazards
| Authors:  Daniel J. Hicks (University of California, Merced)
| Contact: Daniel J. Hicks: 5200 N Lake Road, University of California, Merced, Merced, CA, 95343, USA, <dhicks4@ucmerced.edu>
| Running head: P-value plot
| COI: The author has no conflicts of interest to disclose
| Sources of financial support:  No specific funding has been used in the development of this project
| Data and code:  The complete code and manuscript for this paper is available at <https://github.com/dhicks/p_curve>. The automatically-reproduced analysis can be viewed at <https://dhicks.github.io/p_curve/>. 

\newpage

# Abstract #

  **Background**:  A number of papers by Young and collaborators have criticized epidemiological studies and meta-analyses of air pollution hazards using a graphical method that the authors call a p-value plot, claiming to find zero effects, heterogeneity, and p-hacking.  However, the p-value plot method has not been validated in a peer-reviewed publication. The aim of this study was to investigate the statistical and evidentiary properties of this method. 
  **Methods**: A simulation was developed to create studies and meta-analyses with known real effects $\delta$, integrating two quantifiable conceptions of evidence from the philosophy of science literature. The simulation and analysis is publicly available and automatically reproduced.  
  **Results**: In this simulation, the plot did not provide evidence for heterogeneity or p-hacking with respect to any condition. Under the right conditions, the plot can provide evidence of zero effects; but these conditions are not satisfied in any actual use by Young and collaborators. 
  **Conclusions**: The p-value plot does not provide evidence to support the skeptical claims about air pollution hazards made by Young and collaborators.  


| 
| **What this paper adds**:  This paper uses a simulation approach to examine the statistical and evidentiary properties of the p-value plot, a graphical method that has been used to criticize air pollution epidemiology.  These properties have not been examined in previous peer-reviewed publications.  The results show that the method is incapable of providing evidence to support claims of p-hacking and statistical heterogeneity.  While the method can produce evidence of a zero effect, the method only has this ability under certain conditions.  These conditions are identified, and it is observed that the published criticisms do not satisfy these conditions. 

# Introduction # {#sec:intro}

In numerous recent works [@YoungCerealinducedGenderSelection2009; @YoungAssessingGeographicHeterogeneity2013; @YouPM2OzoneIndicators2018; @YoungReliabilityEnvironmentalEpidemiology2019; @YoungAmbientAirPollution2019; @YoungCombinedBackgroundInformation2019; @YoungEvaluationMetaanalysisAir2019; @YoungReExtendedMortality2020; @YoungPM2AllcauseMortality2020; @YoungParticulateMatterExposure2020; @KindzierskiEvaluationMetaAnalysisAmbient2021; @YoungShiftingSandsUnsound2021], statistician S. Stanley Young and collaborators have criticized epidemiological studies and meta-analyses of the harmful effects of air pollution.  These authors have used a graphical method they call a p-value plot, claiming that this method reveals that zero effects, heterogeneity, and p-hacking are widespread in environmental epidemiology.  

Young and collaborators have drawn highly skeptical conclusions about the hazards of air pollution using their p-value plot method, claiming that "causality of PM10/PM2.5 on heart attacks is not supported" [@YoungReliabilityEnvironmentalEpidemiology2019] and "There is no convincing evidence of an effect of PM2.5 on all-cause mortality" [@YoungPM2AllcauseMortality2020]; Young has advocated that "regulation of PM2.5 should be abandoned altogether" [@YoungSuggestionsEPA2017].  While recent papers in this body of work have not yet been highly cited in the academic literature, one paper by these authors [@YoungAirQualityAcute2017] was cited in the scientific review of US EPA's Ozone Integrated Science Assessment [@USEPACleanAirScientificAdvisoryCommitteeReviewEPAIntegrated2020]<!-- , D-50 -->; this review was conducted while Young was serving on US EPA's Science Advisory Board [@USEPAMembersScienceAdvisory2017].  <!-- In addition, at least some of Young's work on air pollution has been funded by the American Petroleum Institute [@YoungAirQualityEnvironmental2017; @YoungAirQualityAcute2017; @YouPM2OzoneIndicators2018; @YoungReliabilityEnvironmentalEpidemiology2019; @YoungEvaluationMetaanalysisAir2019]. -->  So it is highly plausible that this body of work could be cited in the future.  

However, Young and collaborators have provided only a minimal analysis of the statistical properties of the p-value plot method, in a set of publicly available but non-peer-reviewed notes [@YoungCombinedBackgroundInformation2019].  They have sometimes attempted to justify their method by citing plots of p-values developed by other authors [@SchwederPlotsPvaluesEvaluate1982; @SimonsohnPcurveKeyFiledrawer2014].  But these other plots are designed to answer different questions, are constructed in different ways, and have different statistical properties. 

The aim of this study was to formally evaluate the evidentiary value of the p-value plot method as used by Young and collaborators.  Numerical simulations were chosen for accessibility, extensibility, and speed.  


# Methods #

## The p-value plot ##

The p-value plot is constructed from collections of p-values, often extracted from the primary studies in a meta-analysis.  Features of the plots are interpreted as indicating that there is no underlying effect [@YoungCerealinducedGenderSelection2009; @YouPM2OzoneIndicators2018; @YoungCombinedBackgroundInformation2019; @YoungReExtendedMortality2020; @YoungPM2AllcauseMortality2020], a heterogeneous mixture of zero and nonzero effects [@YoungReliabilityEnvironmentalEpidemiology2019; @YoungAmbientAirPollution2019; @YoungCombinedBackgroundInformation2019; @YoungEvaluationMetaanalysisAir2019], or that the authors of the primary studies have engaged in p-hacking [@YoungEvaluationMetaanalysisAir2019; @YoungShiftingSandsUnsound2021].  (Young and collaborators often refer to the replication crisis literature, but this literature is not immediately relevant to environmental epidemiology.[@HicksOpenScienceReplication2021]) 

Lack of exposition is a major initial challenge in evaluating this method.  The method appears to rely almost entirely on visual inspection, with only vague characterizations of how features of the plots should be interpreted.  These interpretations are rarely given a clear justification, and at least one key term in the exposition is used in a non-standard way.  Therefore, in this section, I develop a more precisely-defined and automatically reproducible set of methods.  Because we wouldn't expect less rigorous methods to provide better evidence, if my rigorous methods do not provide evidence to support the skeptical claims about air pollution made by Young and collaborators, then we should conclude that purely visual methods do not provide evidence either.  

The p-value plot does not appear to have been validated in any peer-reviewed studies.  Young and collaborators do give references to two other graphical methods [@SchwederPlotsPvaluesEvaluate1982; @SimonsohnPcurveKeyFiledrawer2014].  However, these other methods are substantially different from each other and the p-value plot used by Young and collaborators.  The p-curve [@SimonsohnPcurveKeyFiledrawer2014] is a histogram of p-values below the conventional 0.05 threshold, analyzed in terms of its skew, whereas the p-value plot contains points for each individual p-value in a collection (including those above 0.05) and is analyzed in terms of a "hockey stick" shape, "gaps," slope, and linearity.  Simulation methods have been used to validate the p-curve @SimonsohnPcurveKeyFiledrawer2014.  Schweder and Spjøtvoll's p-value plot[@SchwederPlotsPvaluesEvaluate1982] corresponds more closely to the p-value plot used by Young and collaborators, and is also analyzed in terms of its slope; but the two slopes are not in 1-1 correspondence and Schweder and Spjøtvoll assume the p-values come from tests of different hypotheses rather than replications of a test of a single hypothesis.  (The supplement provides formal definitions for all three plots, shows examples using simulated data for a range of real effects that include heterogeneous cases, and explains the differences in detail.)  So citations to these other graphical methods are insufficient to validate the method used by Young and collaborators.  Two works[@YoungCombinedBackgroundInformation2019; @YoungShiftingSandsUnsound2021] include a handful of examples generated using simulated data with zero real effect $\delta = 0$.  However these works do not report any simulations of cases where the real effect was non-zero or heterogeneous, do not report making the simulation code available anywhere for checking reproducibility or extending the analysis, and have not undergone peer review.  

To formally define the p-value plot, we begin with a set of $N$ p-values $\P = \{p_1, p_2, \ldots, p_N\}$, (nominally) produced by applying a given statistical hypothesis test to $N$ replications of a given study design, each replication drawing samples of size $n$ from the given population.  This corresponds to the simplest case of meta-analysis.  Thus the p-values in $\P$ are nominally samples from a single underlying distribution $p_i \sim P$.  Note that, if the real effect is zero $\delta = 0$, then $P$ is the uniform distribution on $[0,1]$.  

Next, let $rank_{asc}(p_i)$ be the (1-indexed) *ascending rank* of $p_i \in \P$, i.e., $rank_{asc}(p_i)$ is the number of p-values $p_j \in \P$ less than or equal to $p_i$.  The smallest p-value has ascending rank 1, and the largest p-value has ascending rank $N$.  Without loss of generality, if \P is already in ascending order $p_1 < p_2 < \cdots < p_N$, then $rank_{asc}(p_i) = i$.  The p-value plot is the graph $(i, p_i)$.  Note that this is equivalent to a rescaled QQ-plot of $\P$ against the uniform distribution, with the theoretical quantiles $q_i = \frac{i}{N} = \frac{rank_{asc}(p_i)}{N}$.  

Young and collaborators explain their interpretation of the plot as follows:  "Evaluation of a p-value plot follows a logical path. Is the p-value plot [sic] homogeneous? If the points roughly [sic] on a 45-degree, they are homogeneous and consistent with randomness; a lessor slope with all points roughly on a line indicates a consistent effect even if some of the individual p-values are not considered statistically significant. If the effects differ, one from another, beyond chance, then the effects are heterogeneous" [@YoungReliabilityEnvironmentalEpidemiology2019 <!-- 48 -->].  

A "45-degree line" typically refers to the graph of an identity function or the line $y = x$, which forms a 45-degree angle with respect to both the x- and y-axes.  This interpretation makes sense for the equivalent QQ-plot, where both axes are on the scale $[0, 1]$; here a slope of 1 indicates that the underlying distribution $P$ is uniformly distributed, which in turn indicates that the real effect is zero.  Strictly speaking, this interpretation does not make sense for the p-value plot as defined, where the x-axis (rank) is on the scale $[1, N]$ and the y-axis (p-value) is on the scale $[0, 1]$.  I will assume that a "45-degree line" typically means that the slope of the equivalent QQ-plot is 1.  It appears that slopes are only evaluated visually; there are no reports of fitting regression models or using any other quantitative methods to measure slopes.  

There are frequent claims that the p-value plot contains a "hockey stick" or "bilinear" pattern [@YoungReliabilityEnvironmentalEpidemiology2019; @YoungAmbientAirPollution2019; @YoungCombinedBackgroundInformation2019; @YoungEvaluationMetaanalysisAir2019; @YoungShiftingSandsUnsound2021 <!-- 106ff -->], which has "small p-values to the lower left ... and then points ascending in a roughly 45-degree line" [@YoungEvaluationMetaanalysisAir2019 <!-- 5 -->].  In this context, "45-degree line" seems to mean that the right-hand side of the plot is linear, even if it does not have a slope of 1.  This non-linear "hockey stick" pattern is taken to indicate some combination of heterogeneous effects, p-hacking, and researcher misconduct.  On two occasions a formal test for non-linearity was conducted by comparing a linear and quadratic regression using an F-test [@YoungReliabilityEnvironmentalEpidemiology2019 <!-- 50 -->; @YoungCombinedBackgroundInformation2019 <!-- 18 -->].  More often there is no explanation of how the "hockey stick" pattern was determined to be present or absent.  

## Simulation design ##

Each run of the simulation comprises $N$ studies, collected together as though for a single meta-analysis.  To facilitate interpretation, each study is based on a two-sample t-test.  Two samples, each of size $n$, are drawn from Gaussian distributions with means $\mu_1 = 0$ and $\mu_2 = \delta$, respectively, and common standard deviation $\sigma = 1$.  (This data generating process is just slightly more complicated for the heterogeneous case; see below.)  $\delta$ corresponds to the true effect size as measured by Cohen's $d = (\mu_2 - \mu_1) / \sigma$ [@SawilowskyNewEffectSize2009]. 

Parameter settings can be systematically varied to compare, e.g., different effect sizes, and multiple runs in each condition allow us to analyze the statistical properties of the p-value plot within and across conditions.  For the primary analysis, 7 different effect sizes are used — corresponding to conventional thresholds for zero, very small, small, moderate, large, and very large effects [@SawilowskyNewEffectSize2009], and a "mixed" or heterogeneous condition.  All studies have the same sample size, $n = 26$.  Compared to the convention of 80% power, this sample size makes the studies severely underpowered to detect the very small (power 3%) and small effects (11%), somewhat underpowered for the moderate effect (42%), adequately powered for the large effect (81%), and overpowered for the very large effect (99%). Each condition is simulated with 500 runs.  The primary conditions examined in the current study are summarized in [@tbl:params].  

| parameter | meaning                             | value(s) |
|:----------|:------------------------------------|---------:|
| $\delta$  | real effect size                    | $0, .01, .2, .5, .8, 1.2, \{0, .8\}$ |
| $\sigma$  | s.d. of samples                     | 1        |
| $n$       | study sample size                   | 26       |
| $N$       | num. studies in each run            | 20       |

Table: **Parameter values used in the current simulation study.**  The real effect size of $\{0, .8\}$ indicates the "mixed" or heterogeneous condition, in which half the population have no response and half the population have a strong response. 500 runs in each combination of parameters are used for the results reported below. {#tbl:params}

In the "mixed" or heterogeneous condition, subsets of the population have different responses to the intervention or exposure.  I use a $\{0, .8\}$ mixture, meaning that one subpopulation has a zero effect or no response to the intervention and the other has a strong response of .8.  These values were chosen to create a mixture that is less likely to look like a homogeneous case; a condition with a $\{.3, .5\}$ mixture would be difficult to distinguish from a homogeneous .4 condition.  

When given this kind of mixture, in the simulation, an individual study draws its sample from one of the subpopulations, selected uniformly at random.  So, in expectation, half of the studies will find a zero effect and half a strong effect, though in any given set of studies there will be variation in the ratios of the two subpopulations.  For simplicity, the simulation currently does not support a continuous mixture.  

The simulation does not have a way to represent p-hacking, publication bias, or researcher misconduct; however, this means that all conditions represent cases in which these factors are absent.  

After generating the data, plots are constructed and analyzed (i) visually for the "hockey stick" pattern, and (ii-iv) using computationally reproducible analyses of gaps, slopes, and linearity.  

To assess "gaps" in the plot (ii), I calculate the largest difference in consecutive p-values $\max(p_{i+1} - p_i)$.  I classify a p-value plot as "gappy" if this largest difference is greater than a threshold value, .125.  That is, if there is at least one visual gap of at least 12.5% in the p-value plot, the plot is considered "gappy."  This threshold was chosen to capture the sense of a visual "discontinuity" in the sequence of p-values.  A continuous measure, size of the largest gap, is reported in the suppement.  

Slopes can be calculated using a simple univariate linear regression.  If the slope of the QQ-plot is approximately 1, this indicates that the underlying distribution of p-values $P$ is uniform, which in turn indicates that the zero hypothesis is true.  Because the slope of the p-value plot is $N$ times the slope of the QQ-plot, analysis of the QQ-plot is simpler than analyzing the p-value plot directly.  

I assess the fitted slope in several ways, to judge whether it is "approximately 1" (iii).   First, I consider simply whether the slope is in the range $1\pm0.1$.  This range was chosen to capture the sense of visual difference from 1.  Next, I apply three hypothesis tests:  a t-test against the null hypothesis that the slope is exactly 1; an equivalence test, the TOST procedure [@LakensEquivalenceTestingPsychological2018], against the null that the slope is outside the range $1\pm0.1$; and a Kolmogorov-Smirnov test (KS-test) against the null of the uniform distribution.  I use the conventional 0.05 threshold for statistical significance for all three tests.  When the estimate is not statistically significant, the tests as implemented in the simulation accept the null hypothesis.  While strictly problematic, this approach simplifies the presentation and analysis of results and aligns with the way hypothesis tests are often used in practice.  Because we are interested in these tests independently — as though the other analyses were not conducted — there is no need to correct for multiple comparisons.  Distributions of slopes are also reported in the supplement.  

Finally I evaluate the linearity of the plot (iv).  Linearity of a plot can be tested by fitting two regression models — one linear, one quadratic — and selecting the model with the better fit.  I use AIC and an F-test for model selection.  I use the conventional 0.05 threshold for statistical significance in the F-test, and accept the null (linear model) when the F-test is not significant.  As with slope, for simplicity I analyze the QQ-plot rather than the p-value plot directly.  In the supplement, area under the curve (AUC) is reported as a continuous measure of linearity, where non-linear QQ-plots plots have lower AUC.  


## Measuring evidence

I use two quantifiable conceptions of evidence that have been widely discussed in the philosophy of science literature, severity and likehood ratios.  For brevity, only severity analysis results are reported here; likelihood ratio results are included in the supplement.  

The severity conception of evidence is associated with philosopher of science Deborah Mayo's reconceptualization of frequentist hypothesis testing [@MayoErrorGrowthExperimental1996; @MayoStatisticalInferenceSevere2018]. Young provided a positive blurb for Mayo's 2018 book @MayoStatisticalInferenceSevere2018, stating that "Her severity requirement [sic] demands that the scientist provide a sharp question and related data. Absent that, the observer should withhold judgment or outright reject." In one work,@YoungShiftingSandsUnsound2021 Young and collaborators cite Mayo @MayoStatisticalInferenceSevere2018 multiple times, and repeatedly characterize the p-value plot as a "severe test."  

I use Mayo's *weak severity criterion*: 

> One does not have evidence for a claim if nothing has been done to rule out ways the claim may be false. If data $x$ agree with a claim [$H$] but the method used is practically guaranteed to find such agreement, and had little or no capability of finding flaws with [$H$] even if they exist, then we have bad evidence, no test (BENT). [@MayoStatisticalInferenceSevere2018 <!-- 5 -->]

On this conception of evidence, a test or analytical method $T$ with observed output $t$ can provide evidence supporting a target hypothesis $H$ only if $T$ would have given a different output if $H$ were false.  (Note that the test $T$ can give any kind of output:  $T = t$ might mean that an test statistic is greater than or equal to some value, or that a plot has some purely qualitative visual feature.)  Hypothesis testing assesses this counterfactual condition using the p-value.  In its most general form, the p-value can be defined as $p = pr(T = t | \lnot H)$, where the role of $\lnot H$ ("not H," the logical negation of the hypothesis of interest) is played by a null hypothesis $H_0$. (Throughout this paper, I am distinguishing a *zero hypothesis* — that some effect is exactly zero — from a *null hypothesis* — the alternate or rival hypothesis used to calculate a p-value.  For instance, the TOST procedure as used in this study has the null hypothesis that the true slope is in the set $\{0.9, 1.1\}$.)  A "small" p-value indicates that the counterfactual is probably true, that is, if $H$ were false then $T$ would probably have given a different output.  On the other hand, a "large" p-value indicates that the test "is practically guaranteed" to produce this output, and so in this case by the weak severity principle "we have bad evidence, no test." 

Severity can also be evaluated qualitatively when a p-value cannot be calculated.  Consider visual features of plots, such as the "hockey stick" pattern that is interpreted as evidence of heterogeneity.  It's not clear how to quantitatively determine whether this visual pattern is present in a given plot.  However, if this pattern is qualitatively common in homogeneous cases, then the weak severity principle implies that this visual pattern does not provide evidence of heterogeneity.  

The severity conception of evidence can be applied to the skeptical claims about air pollution hazards as follows.  The claims $H$ are the zero hypothesis $\delta = 0$, mixture hypothesis $\delta = \{0, .8\}$, and the hypothesis that researchers have engaged in p-hacking or other questionable research practices.  The method or test $T$ is the p-value plot; the outputs $t$ that I will consider are detailed above: (i) the "hockey stick" pattern; (ii) "gaps" in the plot; (iii) slope of approximately 1 (on the QQ-plot); and (iv) non-linearity inferences.  [@Tbl:outputs] summarizes the outputs examined in the current study.  

|     | output            | determined using      | taken as evidence for |
|----:|:------------------|:----------------------|:----------------------|
|   i | "hockey stick"    | visual inspection     | mixed effect          |
|  ii | "gaps"            | visual inspection     | p-hacking or other problems |
|     |                   | largest gap $> .125$  ||
| iii | slope $\approx$ 1 | range $1\pm0.1$       | zero effect           |
|     |                   | T-test not stat. sig. ||
|     |                   | TOST test stat. sig.  || 
|     |                   | KS-test not stat. sig.|| 
|  iv | non-linearity     | AIC: quadratic        | mixed effect
|     |                   | F-test: stat. sig.    ||

Table: **Outputs of the p-value plot examined using the simulation.**  "Outputs" are features of plots that are taken as evidence for critical assessments of air pollution epidemiological studies.  The "determined using" column indicates how these outputs are identified as present/absent in the current study.  {#tbl:outputs}

In the current simulation study, there are p-values at multiple levels.  First there are the p-values produced by the t-tests in the primary studies, the $N$ replications bundled together and plotted in a single p-value plot.  I will call these *primary p-values*.  Second, some of the tests conducted using the p-value plot are themselves statistical hypothesis tests, such as the F-test of linearity.  These tests involve *meta-level p-values*.  Meta-level p-values are tracked by the simulation, but not reported directly here.  Instead, meta-level p-values are compared to the conventional 0.05 threshold, and this comparison is used to produce a meta-level test output such as "non-linear."  Finally, to assess the validity of the p-value plot method, we use the simulation results to estimate the probability of observing a certain meta-level test output value — such as "non-linear" — when the true effect size satisfies a given null hypothesis — such as $\delta = 0.5$.  Because the resulting p-value is being used to assess the validity of the p-value plot method, I will call it a *validation p-value*.  When the validation p-value is "large" (greater than the conventional 0.05), the weak severity criterion implies that this particular use of the p-value plot does not provide evidence for the target hypothesis (heterogeneity or a mixed effect) with respect to the given null hypothesis.  

To generate validation p-values (or qualitatively assess severity for the visual analysis), we need to specify the null hypothesis that plays the role of $\lnot H$.  I will use each of the following:  

a. very small: $\delta = .01$
b. small effect: $\delta = 0.2$
c. moderate effect: $\delta = 0.5$
d. strong effect: $\delta = 0.8$
e. very strong effect: $\delta = 1.2$
f. greater than zero: $a \vee \cdots \vee e$
g. mixed effect: $\delta = \{0, .8\}$ for $H: \delta = 0$ and vice-versa (that is, the other skeptical hypothesis)
h. any other effect: $f \vee g$ (any of the non-zero effects or the other skeptical hypothesis)

In each case, when the validation p-value is "large" $p > .05$, the simulation results indicate that this test output is common in the null hypothesis case, and so the weak severity criterion implies that the p-value plot would not provide evidence to support the skeptical claim made by Young and collaborators.  

## Reproducibility ## {#sec:reproducibility}

The simulation, analysis, and outputs (figures and tables) are publicly available and automatically reproduced.  Code is available at <https://github.com/dhicks/p_curve> and the automatically reproduced analysis can be viewed at <https://dhicks.github.io/p_curve/>. 

The simulation and analysis were both written in R version 4.1.0 [@RCoreTeamLanguageEnvironmentStatistical2018] and make extensive use of the `tidyverse` suite of packages, version 1.3.1 [@WickhamTidyverseEasilyInstall2019].  The TOSTER package version 0.3.4 [@LakensTOSTERTwoOneSided2018] was used to conduct the TOST analysis.  Because the software on the virtual machine used to automatically reproduce the analysis is updated each time the analysis is re-run, software versions reported online may be different from those reported here. 


# Results 

[@Fig:samples_young] shows 35 examples of the p-value plot across the seven conditions, and [@fig:composite_young] shows the p-value plot across all runs of the simulation.  

![**Examples of the p-value plot.**  Drawn at random from the simulation results.  Rows and colors correspond to conditions or real effects ($\delta$), from zero (0) to very strong (1.2) and a mixed condition $\delta = \{0.0, 0.8\}$.  Columns correspond to indices for the simulation runs that produced these results, and are not meaningful.  Each point corresponds to a single p-value in the meta-analysis (simulation run); the x-axis is the ascending rank of the p-value in the set \P, and the y-axis is the p-value itself.](fig_1_samples_young.png){ #fig:samples_young width=6in height=4in }

![**Composite of the p-value plot.**  Each individual line is the p-value plot for a single run of the simulation; all simulation runs are shown here.  Panels correspond to conditions or real effects ($\delta$), from zero (0) to strong (1.2) and a mixed condition $\delta = \{0.0, 0.8\}$.](fig_2_young_composite.png){ #fig:composite_young width=6in height=4in }

There are substantial qualitative differences within effect sizes as well as similarity across consecutive effect sizes.  Except for the very large (overpowered) effect size, there tends to be both statistically significant and insignificant primary p-values.  Larger effect sizes have more statistically significant primary results, resulting in a series of small primary p-values that gradually bend up.  However, the top row of [@fig:samples_young] and the first panels of [@fig:composite_young] indicates that even zero effects can look non-linear.  Comparing the composite plots ([@fig:composite_young]) for the mixed condition and the moderate effects condition, on average the mixed condition tends to produce a sharper bend upwards.  However, [@fig:samples_young] indicates that an individual moderate effects plot can have a sharp bend (index 290) and a mixed effects plot can have a gradual bend (index 66).  

## Visual analysis

I focus on two visual patterns that are frequently discussed in the critiques of air pollution epidemiology:  (i) the "hockey stick" pattern; (ii) "gaps" in the plot.  

The "hockey stick" pattern is taken to be evidence of heterogeneity.  The "hockey stick" comprises a more-or-less flat series of small p-values on the left (supposedly corresponding to the mixture component with a real effect) and a second series of steeply increasing p-values on the right (supposedly corresponding to the mixture component with zero effect).  This pattern is visible in all of the example plots for the moderate, strong, and mixed effects (rows 4, 5, and 7), and arguably also with small effects (row 3).  Because the "hockey stick" pattern appears in plots where there is only a single homogeneous effect, the weak severity criterion implies that the hockey stick pattern does not provide evidence of a mixed or heterogeneous case.  

Visual "gaps" in the plot are taken to be evidence of p-hacking, publication bias, or other questionable researcher practices [@YoungAssessingGeographicHeterogeneity2013; @YoungCombinedBackgroundInformation2019; @YoungEvaluationMetaanalysisAir2019].  These gaps are common in [@fig:samples_young], across all conditions except the very strong effect (where almost all values are below 0.05).  Note that the simulation does not include p-hacking, publication bias, or other questionable researcher practices.  Thus the weak severity criterion implies that gaps in the plot do not provide evidence of p-hacking, publication bias, or other questionable researcher practices.  


## Computationally reproducible analyses ##

@Fig:evidence_severity shows the results of the severity analysis as validation p-values; see the supplement for a table version of these results.  By the weak severity criterion, when the results of the severity analysis are greater than .05 (above the dashed line), the test output does not provide evidence in support of the target hypothesis.  

![**Results of the severity analysis for outputs (ii) gaps in the plot, (iii) slope of 1, and (iv) non-linearity.**  Severity analysis results are reported as validation p-values:  small values (conventionally $< .05$) indicate a severe test with respect to the null hypothesis.  Panels correspond to null hypotheses, and y-axis values correspond to the severity assessment (as a validation p-value) for the output with respect to the given null hypothesis.  The dashed line indicates $p = 0.05$; *points below this line indicate severe tests*.](fig_3_evidence_severity.png){#fig:evidence_severity width=6in height=4in}

Unspecified problems of p-hacking, publication bias, or other questionable research practices are supposedly supported by gaps in the plot.  The validation p-value for the presence of these gaps (output ii-gap) is greater than .25 for every null hypothesis except the very strong effect, indicating that gaps are quite common.  In all conditions except the large and very large effects, the majority of p-value plots are "gappy."  Thus gaps in the plot do not provide evidence of p-hacking, publication bias, or other questionable researcher practices.  The distribution of the size of the largest gap across each real effect size is reported in the supplement.  Most plots have gaps larger than .125 across all conditions except the very strong effect.  

The zero effect hypothesis is supposedly supported by a slope of approximately 1 on the QQ-plot.  All methods are severe against large and very large effects; only the TOST and KS tests are severe against moderate effects; and no methods are severe against very small or small effects.  The TOST test also may be severe against the greater-than-zero and non-zero hypotheses.  So the "45-degree line" may or may not provide evidence for zero effects, depending on the particular null hypothesis being tested and particular test used. 

This means that, in order to determine whether the "45-degree line" provides evidence for zero effects, the analyst must be explicit about both the test method used and the null hypothesis being tested.  A slope of .92 or statistically significant TOST test (meta-level p-value less than 0.05), for example, will not provide evidence if the null hypothesis is a very small effect.  In addition, while these results recommend using the TOST test, this test requires an explicit equivalence interval, the range of values that are "approximately 1."  But Young and collaborators are not explicit in any of these requisite ways. They do not report calculating a slope, using a regression line, or other method, much less conducting some further analysis of that slope. I take it that the lack of detail is a good reason to think that they have simply assessed the slope visually.  A visual assessment will be less sensitive than calculating a slope and determining whether it is in the range $1\pm.1$, which can only provide evidence against a relatively large null hypothesis.  So in cases where a small effect is a live possibility, insofar as the p-value plot is interpreted using visual judgment alone, the "45-degree line" does not provide evidence of a zero effect.  The distribution of slopes across each real effect size is reported in the supplement.  Slopes in the range $1\pm.1$ are common for very small, small, and mixed effects.  

The mixed-effect hypothesis, or heterogeneity, is supposedly supported by non-linearity.  Two quantified versions of this output are examined here, comparing linear and quadratic regressions using AIC (iv-AIC) and an F-test (iv-F).  The AIC and F-test evaluations of non-linearity do not provide severe tests against any of the alternative conditions.  That is, neither of the tests of linearity considered here provide evidence for heterogeneity.  In the supplement, the area under the curve (AUC) of the QQ-plot is used as a continuous measure of non-linearity, and distributions of AUC values are reported across each real effect size.  Lower AUC values indicate further deviation from the line $x=y$.  The AUC distribution for the mixed effect overlaps substantially with the distributions for all other effect sizes except very large.  


# Discussion 

This simulation analysis finds that the p-value plot does not provide evidence for heterogeneity or p-hacking based on the "hockey stick" shape, "gaps" in the plot, or AIC or F-tests of non-linearity.  The method can provide evidence for zero effects based on a slope of 1, depending on what rival or null hypothesis is considered and how the plot is analyzed.  In general, producing evidence for zero effects against small effects requires using the TOST approach; visual inspection alone will not be sufficient.  This approach requires setting an explicit range of values within which the slope is considered "approximately 1."  

In the meta-analyses criticized by Young and collaborators, the estimated short-term effects of air pollution are small or very small on a relative risk scale. For example,   Nawrot et al. estimated the effect for air pollution (increase of 10 $\mu g/m^3$ PM$_{10}$) on non-fatal myocardial infarction to be 1.02 (95% CI 1.01–1.02)@NawrotPublicHealthImportance2011; the point estimates for six pollutants reported by Mustafić et al. (increase of 10 $\mu g/m^3$ for all except carbon monoxide) were all in the range 1.003-1.048 @MustaficMainAirPollutants2012; Liu et al. estimated effects for PM$_{10}$ (increase of 10 $\mu g/m^3$) on all-cause mortality of 1.044 (95% CI 1.039-1.050) @LiuAmbientParticulateAir2019; and Orellano et al. estimated effects for PM$_{2.5}$ on all-cause morality of 1.0065 (95% CI 1.0044–1.0086) @OrellanoShorttermExposureParticulate2020.  A precise conversion from risk ratios to Cohen's $d$ is beyond the scope of this paper.  However, using a rule of thumb that the risk ratio is approximately equal to the log odds when the outcome is rare [@VieraOddsRatiosRisk2008] and the conversion factor $\sqrt{3}/\pi$ between the log odds and Cohen's $d$, a risk ratio of 1.05 is roughly equivalent to $d=0.03$, which is only slightly larger than the very small effect condition examined here. (An additional analysis in the automatically reproduced analysis document examines conditions with a very small real effect sizes of $\delta = 0.05$ and varying power to detect these very small effects.)  Visual inspection of the p-value plot alone is incapable of producing evidence against very small epidemiological effects of air pollution.  

More often, Young and collaborators have claimed to find evidence of heterogeneity, p-hacking, and publication bias [@YoungReliabilityEnvironmentalEpidemiology2019; @YoungAmbientAirPollution2019; @YoungEvaluationMetaanalysisAir2019].  The simulation results indicate that the p-value plot is incapable of providing evidence for any of these claims, using either visual inspection or either of the quantitative approaches examined here.  The features that Young and collaborators point to — the "hockey stick" shape, "gaps," non-linearity — are readily produced by moderate and stronger effects, and can even appear in zero and very small effect conditions.  

All together, the p-value plot method cannot support the skeptical claims about air pollution epidemiology made by Young and collaborators.  
