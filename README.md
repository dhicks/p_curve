![Travis build status](https://travis-ci.org/dhicks/p_curve.svg?branch=master)

# Analysis for "Young's p-value plot does not provide evidence against air pollution hazards"
## Daniel J. Hicks
## University of California, Merced
## <hicks.daniel.j@gmail.com>

This package includes the complete code and manuscript for the paper "Young's p-value plot does not provide evidence against air pollution hazards" by Daniel J. Hicks.  

After cloning the repository and installing the R dependencies listed in `DESCRIPTION`, the analysis, output files, and `html` version of the analysis script can be recreated from the command line by calling `Rscript -e "rmarkdown::render('run_metas.R')"` within the `scripts` directory.  Alternavely, if `make` is installed, call either `make script` from the command line within the top folder of the repository or `make knit` from within the `scripts` folder.   

This reproduction of the study is conducted automatically using [Travis-CI](https://travis-ci.org/github/dhicks/p_curve) and the result is automatically published at <https://dhicks.github.io/p_curve/>. 

Note that reproducing the analysis creates a top-level `out` folder for plots and tables if this folder does not already exist, and quietly overwrites existing versions of the output files in this folder.  

The workhorse functions to conduct and analyze the simulation are all found in the `p.curve` folder, which is loaded as a local package using `devtools`.  Documentation for all of these functions can be constructed by calling `devtools::document(file.path('..', 'p.curve'))` from within the `scripts` folder, and then queried in the usual way, e.g., `?many_metas`. 
