if (!require(renv)) {
    install.packages('renv')
}

## Restore the other dependencies
renv::restore(exclude = 'p.curve')

## Recreate documentation for the local package
# devtools::document(file.path('..', 'p.curve'))
## Install the local package
renv::install(file.path('..', 'p.curve'))
