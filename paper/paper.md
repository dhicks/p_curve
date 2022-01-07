---
title: "Supplement to \"The p-value plot does not provide evidence against air pollution hazards\""
author: "Daniel J. Hicks"
bibliography: Young.bib
header-includes:
  - \usepackage{booktabs}
  - \usepackage{longtable}
tblPrefix:
  - "table"
  - "tables"
---

\renewcommand{\P}{\ensuremath{\mathbb{P}}}
\newcommand{\SchSp}{Schweder and Spjøtvoll}




# Supplemental methods #

## The three plots ##

Young and collaborators frequently cite two other graphical methods [@SchwederPlotsPvaluesEvaluate1982; @SimonsohnPcurveKeyFiledrawer2014].

All three graphical methods take as input a set of $N$ p-values $\P = \{p_1, p_2, \ldots, p_N\}$.  In @SchwederPlotsPvaluesEvaluate1982 these are taken from separate tests of $N$ different hypotheses.  In @SimonsohnPcurveKeyFiledrawer2014 and the works by Young and collaborators, the p-values are (nominally) produced by applying a given statistical hypothesis test to $N$ replications of a given study design, each replication drawing samples of size $n$ from a given population.  This corresponds to the simplest case of meta-analysis.  Thus the p-values in $\P$ are nominally samples from a single underlying distribution $p_i \sim P$.  Note that, if the real effect is zero $\delta = 0$, then $P$ is the uniform distribution on $[0,1]$.  

## Schweder and Spjøtvoll's p-value plot ##

Young and collaborators have frequently cited the "p-value plot" presented in @SchwederPlotsPvaluesEvaluate1982.  For this p-value plot, let $rank_{desc}(p_i)$ be the (1-indexed) *descending rank* of $p_i \in \P$, i.e., $rank_{desc}(p_i)$ is the number of p-values $p_j \in \P$ greater than or equal to $p_i$.  The largest p-value has descending rank 1, and the smallest p-value has descending rank $N$.  Then Schweder and Spjøtvoll's p-value plot plots the graph $(1-p_i, rank_{desc}(p_i))$.  See @fig:samples_schsp. 

![Examples of \SchSp's p-value plot, drawn at random from the simulation results [@SchwederPlotsPvaluesEvaluate1982].  Rows and colors correspond to conditions or real effects ($\delta$), from zero (0) to moderate-strong (0.6) and a mixed condition $\delta = \{0.0, 0.6\}$.  Columns correspond to indices for the simulation runs that produced these results, and are not meaningful.  (In particular, there is no relationship between simulation run $j$ in condition $a$ and simulation run $j$ in condition $b$.) In \SchSp's p-value plot, each point corresponds to a single p-value in the meta-analysis (simulation run); the y-axis is the p-value itself and the x-axis is the descending rank of the p-value.](../out/samples_schsp.png){#fig:samples_schsp}

@SchwederPlotsPvaluesEvaluate1982 give a brief (and rather informal) argument that the relationship between $1-p_i$ and $rank_{desc}(p_i)$ should be approximately linear "when [$p_i$] is not too small" and that, from left to right, "often, the plot will not show a clearcut break but rather a gradual bend" away from linearity.  Their argument concludes that "the slope of that straight line is an estimate of ... the number of true null hypotheses" in \P, and so the former can be used to estimate the latter [@SchwederPlotsPvaluesEvaluate1982 494].  Schweder and Spjøtvoll's p-value plot is designed around different assumptions and to answer a different question than either of the other two plots.  Also, Schweder and Spjøtvoll's p-value plot generally ignores "small" p-values, which roughly correspond to the statistically significant p-values.  They illustrate their method with example data and a straight line "drawn by visual fit" rather than regression analysis.  


## Simonsohn, Nelson, and Simmons' p-curve ##

Young and collaborators have regularly echoed concerns about the replication crisis unfolding in social psychology and certain areas of biomedical research [for example, @YoungReliabilityEnvironmentalEpidemiology2019 50].  In particular, they appeal to concerns about p-hacking [@SimonsohnPcurveKeyFiledrawer2014].  Note that, other than Young's p-value plot, Young and collaborators have provided no specific[^phacking] empirical[^counting] evidence of p-hacking in environmental epidemiology, and to my knowledge no such evidence has been published [@HicksOpenScienceReplication2021]. 

[^phacking]: Young and collaborators do cite @HeadExtentConsequencesPHacking2015, which used text mining methods to examine p-values reported in "all Open Access papers available in the PubMed database," classified at the journal level into 22 disciplines.   @HeadExtentConsequencesPHacking2015 did not report the full size of their sample, but did include tens of thousands of p-values from "Medical and health sciences," which likely includes epidemiology but also other fields with very different methods, e.g., small-n animal model experiments and industry-funded clinical trials.  While their statistical tests did find evidence of p-hacking in "Medical and health sciences," they did not consider subfields.  Based on the methodological and funding diversity of the fields covered, and the lack of information about the distribution of subfields within their sample, we should be hesitant about drawing inferences to subfields.  

Young and collaborators have frequently associated Young's p-value plot with a method developed to detect p-hacking, called a "p-curve" [@SimonsohnPcurveKeyFiledrawer2014; for a comparison of several methods to detect p-hacking, see @McShaneAdjustingPublicationBias2016].  The intuition behind the p-curve is that p-hacking will tend to produce an excess number of p-values "just below" the conventional $0.05$ threshold for statistical significance.  Formally, Simonsohn et al.'s p-curve first divides the interval $[0, 0.05]$ into 5 bins at the thresholds $0.01, 0.02, 0.03, 0.04, 0.05$, then calculates $N_b$, the number of p-values in bin $b$.  The p-curve is the graph $(b_t, N_b)$, where $b_t$ is the threshold for bin $b$.  The method then formally tests for p-hacking by applying statistical tests of the null hypothesis that the restricted distribution $P|_{p < 0.05}$ is uniform.  See @fig:samples_simonsohn. 

![Examples of Simonsohn et al.'s p-curve, drawn at random from the simulation results [@SimonsohnPcurveKeyFiledrawer2014].  Rows and colors correspond to conditions or real effects ($\delta$), from zero (0) to moderate-strong (0.6) and a mixed condition $\delta = \{0.0, 0.6\}$.  Columns correspond to indices for the simulation runs that produced these results, and are not meaningful.  (In particular, there is no relationship between simulation run $j$ in condition $a$ and simulation run $j$ in condition $b$.) Simonsohn et al.'s p-curve is restricted to p-values below the conventional 0.05 threshold; empty plots correspond to cases in which no p-values were below the threshold.  These p-values are binned at the thresholds $0.01, 0.02, 0.03, 0.04, 0.05$, and each point corresponds to the number of p-values in the given bin.  This plot is equivalent to a histogram.](../out/samples_simonsohn.png){#fig:samples_simonsohn}

Note that the p-curve is a histogram on the interval $[0, 0.5]$ with binwidth $0.01$, and that it only includes statistically significant p-values.  Thus the p-curve and Schweder and Spjøtvoll's p-value plot not only produce different kinds of plots, but actually direct their attention to different — typically disjoint — subsets of p-values.  The two kinds of plots cannot be equivalent to each other.  At no point have Young and collaborators acknowledged this fundamental difference between the two methods that they cite as support for their own method.  

## Young's p-value plot ##

For Young's p-value plot, let $rank_{asc}(p_i)$ be the (1-indexed) *ascending rank* of $p_i \in \P$, i.e., $rank_{asc}(p_i)$ is the number of p-values $p_j \in \P$ less than or equal to $p_i$.  The smallest p-value has ascending rank 1, and the largest p-value has ascending rank $N$.  Without loss of generality, if \P is already in ascending order $p_1 < p_2 < \cdots < p_N$, then $rank_{asc}(p_i) = i$.  And Young's p-value plot is the graph $(i, p_i)$.  

Statistically-minded readers might have already noted that Young's p-value plot is a rescaled QQ-plot of $\P$ against the uniform distribution, with the theoretical quantiles $q_i = \frac{i}{N} = \frac{rank_{asc}(p_i)}{N}$.  It is not equivalent to the other two plots, and so cannot be validated by references to them.  First, $rank_{desc}(p_i) = N - rank_{asc}(p_i) + 1$, and so for a fixed number of studies $N$ Young's p-value does have a 1-1 mathematical relationship to Schweder and Spjøtvoll's p-value plot.  In geometric terms, Young's p-value plot swaps the axes of Schweder and Spjøtvoll's p-value plot and reverses the direction of the ranking.  However, a regression line fit to Schweder and Spjøtvoll's p-value plot will not deterministically correspond to a regression line fit to Young's p-value plot; see discussion in the supplemental results.  In addition, the properties that Young and collaborators use in their analysis of their p-value plots do not correspond to the properties used by Schweder and Spjøtvoll (which, again, are justified with informal arguments rather than a formal analysis).  So, even if Schweder and Spjøtvoll's p-value plot can be considered validated for certain purposes (it should be clear that I'm skeptical on this point), this does not validate Young's p-value plot as used by Young and collaborators.  Second, Simonsohn, Nelson, and Simmons' p-curve constructs a histogram on a subset of p-values; this is a completely different construction from Young's p-value plot, and so citations to the former also do not validate the latter.  

## Likelihood ##

The likelihood conception of evidence is not strongly associated with any one statistician or philosopher of science, though it can be associated with one approach to Bayesian statistics [@KassBayesFactors1995; @RomeijnPhilosophyStatistics2017].  

Formally, the likelihood conception of evidence compares two rival hypotheses $H_1$ and $H_2$ using some data $d$.  The likelihood ratio is defined as
$$ K(H_1, H_2; d) = \frac{L(H_1; d)}{L(H_2; d)} = \frac{pr(d | H_1)}{pr(d | H_2)}. $$
If $K > 1$, then the evidence favors $H_1$; and $K < 1$ then the evidence favors $H_2$.  Sometimes $\log K$ is used to create symmetry between $H_1$ and $H_2$.  On one common interpretive scale, $\left|\log_{10} K\right| < 0.5$ is "not worth more than a bare mention," not supporting either hypothesis; $0.5 < \left|\log_{10} K\right| < 1$ is "substantial"; $1 < \left|\log_{10} K\right| < 2$ is "strong"; and $2 < \left|\log_{10} K\right|$ is "decisive" [@KassBayesFactors1995]. 

To apply the likelihood conception of evidence to Young and collaborators' skeptical claims about air pollution, $H_1$ will be the zero or mixture hypothesis, the rival hypothesis $H_2$ will be the hypotheses (a-h), and the data $d$ will be the analysis outputs (i-iv).  (For simplicity, the same dichotomous frequentist test outputs are used, e.g., statistically significant or not, rather than continuous-valued likelihoodist or Bayesian alternatives.)  In each case, insofar as $K < 0.5$, this implies that the p-value plot does not provide evidence to support the zero or mixture hypotheses.  



# Supplemental results #

## Gaps 

@Fig:gaps shows the distribution of sizes of the largest gap across effect sizes.  

![**Distribution of gaps**:  Each violin plot shows the distribution of sizes of the largest gap for each real effect size.  The horizontal line indicates the threshold .125 for a plot to be considered "gappy."  Except for the very large effect size, most plots across almost all effect sizes have gaps.](../out/gaps.png){#fig:gaps}

## Slopes

@Fig:slopes and @fig:slopes_qq show the distribution of slopes of regression lines fit to the QQ-plot, \SchSp's p-value plot, and Young's p-value plot across effect sizes.  @fig:slopes_qq uses the same data as the left panel of @fig:slopes, with a wider aspect ratio for readability and dashed lines indicating the thresholds for a slope of "approximately 1."  

![**Distribution of slopes**:  Each violin plot shows the distribution of slopes for the linear regression fit to each plot across each real effect size.  Note that the slopes for Young's p-value plot are rescaled from the QQ-plot.](../out/slopes.png){#fig:slopes}

![**Slopes for the QQ-plot**: Data are the same as in the left panel of @fig:slopes.  These slopes are used for the slope analysis that supposedly gives evidence of zero effects.  The dashed lines indicate the $1\pm.1$ threshold used for "approximately 1."](../out/slopes_qq.png){#fig:slopes_qq}

@Fig:slopes_scatter shows the relationship between the slopes of Young's and \SchSp's p-value plots.  For an ordinary least-squares regression of $y$ against $x$, the slope of the fitted regression line is $r s_x / s_y$, where $r$ is the estimated correlation coefficient between $x$ and $y$ and the $s$ are the estimated standard deviations.  So the slope of $x$ and $y$ (swapping the axes) is $r s_y / s_x$.  The ratio of the first slope to the reciprocal of the second slope is $r^2$.  This might suggest a deterministic relationship between the slopes of Schweder and Spjøtvoll's and Young's p-value plots.  However, the estimated correlation coefficient $r$ is calculated from observations drawn from the random variables $X$ and $Y$, and so $r$ itself is an observation drawn from a random variable.  That is, the observed value of $r$ will vary between different iterations of the study.  So the relationship between the slopes of the two p-value plots is noisy.  See @fig:slopes. 

![**Relationship between slopes of Young's and \SchSp's p-value plots**:  Scatterplot of the slopes for the two p-value plots, log-log scale.  The two plots are in a 1-1 relationship with each other, by reversing an axis and swapping the x and y axes.  But the slopes of the regression lines fit to each plot are not in a 1-1 relationship.](../out/slopes_scatter.png){#fig:slopes_scatter}

## Linearity

@Fig:auc shows the distribution of area under the curve (AUC) of the QQ-plot across each real effect size.  AUC has its maximum (expected) value of 0.5 when the real effect size $\delta = 0$, and decreases as the curve of the QQ-plot bends away from linearity.  So AUC can be used as a continuous measure of linearity.  

![**Area under the curve (AUC) of the QQ-plot**:  Each violin plot shows the distribution of area under the curve (AUC) across each real effect size.  When the real effect size is $\delta = 0$, distribution of p-values is uniform and the expected QQ-plot is a straight line $y=x$.  The area under this perfectly linear QQ-plot is .5, as indicated by the median in the corresponding violin plot.  The AUC distribution for the very small effect is almost identical to the zero effect, and expected AUC decreases as real effect size increases.  The distribution for the mixed or heterogenous effect overlaps with the distributions for almost every other effect size, except for very large effects.](../out/auc.png){#fig:auc}


## Severity analysis

\input{../out/severity.tex}


## Likelihood analysis

@Fig:evidence_likelihood_zero and @fig:evidence_likelihood_mix show the results of the likelihood analysis; see the supplemental materials for a tables and interactive versions of these results.  Log likelihood ratios are reported, so results above $0.5$ support $H_1$ and results below $-0.5$ support $H_2$.  So, by the weak severity criterion, when the results of the likelihood analysis are $< 0.5$, the test output does not provide evidence supporting the target hypothesis of zero or mixed effect.  

![**Results of the likelihood analysis for $H_1: \delta = 0$.**  Each point gives the log likelihood ratio for $H_1$ vs. rival hypothesis, given an output.  Each panel represents one comparison of $H_1$ against a rival hypothesis $H_2$.  Position on the y-axis indicates the strength of the evidence that the output provides to the hypotheses:  *greater values indicate more support for $H_1$ over $H_2$*.  (Points at the plot margins have infinite value due to division by zero.)  Shaded regions indicate the degree of support for one hypothesis against the other, in order from lightest to darkest:  none, "substantial," "strong," "decisive."  An interactive version of this plot is included in the automatic reproduction of the analysis for this paper.](fig_4_evidence_likelihood_zero.png){#fig:evidence_likelihood_zero  width=6in height=4in }

![**Results of the likelihood analysis for $H_1: \delta$ is mixed.**  Interpretation is the same as @fig:evidence_likelihood_zero.](fig_5_evidence_likelihood_mix.png){#fig:evidence_likelihood_mix  width=6in height=4in }

For "gaps" in the plot, calculating the likelihood ratio would require simulation conditions that included p-hacking and other questionable research practices.  Because the simulation does not currently support these kinds of conditions, likelihood analysis cannot be used for this output.  
 
The zero effect hypothesis is supposedly supported by a slope of approximately 1.  All four methods provide "decisive" support for a zero effect against strong and very strong effects, and all except the T-test provide "substantial" or better support against moderate effects.  The TOST approach provides stronger or equally strong evidence, compared to the other approaches, across all of the rival hypotheses.  When the rival hypothesis includes small or very small effects, other approaches either do not provide evidence to support zero effects.  So, as with the severity analysis, *whether and to what degree the "45-degree line" might provide evidence for zero effects depends on the choice of rival hypothesis and analytical approach used*. In addition, and again as in the severity analysis, the visual assessment apparently used by Young and collaborators will provide weaker evidence than the range test, which does not provide evidence against rivals that include small effects.  *So, against rival hypotheses that include small effects, insofar as Young and collaborators are relying on visual judgment, the "45-degree line" does not provide evidence of a zero effect.* 

For the mixed effect hypothesis, all of the points for both AIC and the F-test are in the "not worth mentioning" or no evidence range, and so *neither method provides evidence to support heterogeneity*.  This is the same conclusion reached by the severity analysis. 

\input{../out/likelihood.tex}


# References #

::: {#refs}
:::
