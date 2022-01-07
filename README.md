[![build](https://github.com/dhicks/p_curve/actions/workflows/build.yml/badge.svg)](https://github.com/dhicks/p_curve/actions/workflows/build.yml)

# Analysis for "The p-value plot does not provide evidence against air pollution hazards"
## Daniel J. Hicks
## University of California, Merced
## <hicks.daniel.j@gmail.com>

This package includes the complete code and manuscript for the paper "The p-value plot does not provide evidence against air pollution hazards" by Daniel J. Hicks.  

After cloning the repository, `scripts/setup.R` will install the necessary dependencies using `renv`.  Alternatively, manually install the dependencies listed in `DESCRIPTION`, or if GNU Make is available call `make setup`.  The analysis, output files, and `html` version of the analysis script can be recreated from the command line by calling `Rscript -e "rmarkdown::render('run_metas.R')"` within the `scripts` directory.  Alternavely, if GNU Make is available, call either `make script` from the command line within the top folder of the repository or `make knit` from within the `scripts` folder.   

This reproduction of the study is conducted automatically using GitHub Actions and the result is automatically published at <https://dhicks.github.io/p_curve/>. 

Note that reproducing the analysis creates a top-level `out` folder for plots and tables if this folder does not already exist, and quietly overwrites existing versions of the output files in this folder.  

The workhorse functions to conduct and analyze the simulation are all included in the package `p.curve`, which is included in the repository and installed using `renv::install()`.  Documentation for all of these functions can be constructed by calling `devtools::document(file.path('..', 'p.curve'))` from within the `scripts` folder, and then queried in the usual way, e.g., `?many_metas`. 

