---
title: "Supplemental materials for 'Young's p-value plot does not provide evidence against air pollution hazards'"
author: "Daniel J. Hicks"
institute: "University of California, Merced"
email: hicks.daniel.j@gmail.com
bibliography: Young.bib

header-includes:
  - \usepackage{longtable}
  - \usepackage{booktabs}
---

\newcommand{\SchSp}{Schweder and Spjøtvoll}

# Sample plots #

[@Fig:samples_schsp] shows examples of \SchSp's p-value plot; [@fig:samples_simonsohn] shows examples of Simonsohn et al.'s p-curve.  Note that these plots are quite different from each other and Young's p-value plot.  In particular, while \SchSp's p-value plot is a reflection of Young's p-value plot, the exchange of axes means that regression lines fit to one do not necessarily correspond to regression lines fit to the other, and so there is not a 1-1 relationship between their slopes.  

![Examples of \SchSp's p-value plot, drawn at random from the simulation results [@SchwederPlotsPvaluesEvaluate1982].  Rows and colors correspond to conditions or real effects ($\delta$), from zero (0) to moderate-strong (0.6) and a mixed condition $\delta = \{0.0, 0.6\}$.  Columns correspond to indices for the simulation runs that produced these results, and are not meaningful.  (In particular, there is no relationship between simulation run $j$ in condition $a$ and simulation run $j$ in condition $b$.) In \SchSp's p-value plot, each point corresponds to a single p-value in the meta-analysis (simulation run); the y-axis is the p-value itself and the x-axis is the descending rank of the p-value.](../out/samples_schsp.png){#fig:samples_schsp}

![Examples of Simonsohn et al.'s p-curve, drawn at random from the simulation results [@SimonsohnPcurveKeyFiledrawer2014].  Rows and colors correspond to conditions or real effects ($\delta$), from zero (0) to moderate-strong (0.6) and a mixed condition $\delta = \{0.0, 0.6\}$.  Columns correspond to indices for the simulation runs that produced these results, and are not meaningful.  (In particular, there is no relationship between simulation run $j$ in condition $a$ and simulation run $j$ in condition $b$.) Simonsohn et al.'s p-curve is restricted to p-values below the conventional 0.05 threshold; empty plots correspond to cases in which no p-values were below the threshold.  These p-values are binned at the thresholds $0.01, 0.02, 0.03, 0.04, 0.05$, and each point corresponds to the number of p-values in the given bin.  This plot is equivalent to a histogram.](../out/samples_simonsohn.png){#fig:samples_simonsohn}


# Linearity analyses #

The distribution of linearity calls for all three methods across conditions are shown in @fig:linearity and table \ref{tab:linearity}. 

![Distributions of outcomes for the linearity tests across conditions.  Each panel corresponds to one method of testing linearity:  AIC comparison of linear vs. quadratic regression, F-test of quadratic vs. linear regression, and a Kolmogorov-Smirnov test against the uniform distribution.  Yellow regions indicate the share of cases (simulation runs) in which the test concludes Young's p-value plot is non-linear; blue regions indicate the share of cases in which the test concludes the plot is linear.](../out/linearity.png){#fig:linearity}

\input{../out/linearity.tex}


# Severity analysis #

\input{../out/severity.tex}

# Likelihood analysis #

\input{../out/likelihood.tex}
