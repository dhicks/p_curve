P.CURVE := $(shell find ../p.curve -name "*.R")

all: knit

knit: run_metas.html
run_metas.html: run_metas.R \
                $(P.CURVE)
	Rscript -e "rmarkdown::render('run_metas.R')"

