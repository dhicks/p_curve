## As of 2022-01-13 need to use dev version of pandoc
## <https://github.com/jgm/pandoc/commit/4fdbb30a97d0b10d64888783d803d1654da1975d>
## But this breaks other stuff?  Runs in a regular command line after setting up PATH
# PATH=~/.local/bin:$PATH

PAPER := paper
OUT := out
SCRIPT := scripts

all: script summary paper zip title

setup: 
	cd scripts; Rscript setup.R

script:
	$(MAKE) -C $(SCRIPT)


supplement: $(PAPER)/supplement.pdf
$(PAPER)/supplement.pdf: $(PAPER)/header.yaml $(PAPER)/supplement.md $(PAPER)/Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	cd $(PAPER); pandoc header.yaml supplement.md -o supplement.pdf --filter=/usr/local/bin/pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve


## Figures
$(PAPER)/fig_%.png: $(OUT)/fig_%.png
	cp -a $^ $(PAPER)

## Grayscale figures
bw: $(wildcard $(PAPER)/fig_*_bw.png)
$(PAPER)/fig_%_bw.png: $(PAPER)/fig_%.png
	convert $^ -set colorspace gray -quiet $@ 

paper: $(PAPER)/paper.pdf
$(PAPER)/paper.pdf: $(PAPER)/header.yaml $(PAPER)/paper.md \
                    $(PAPER)/Young.bib \
                    $(wildcard $(PAPER)/fig_*.png)
	cd $(PAPER); pandoc header.yaml paper.md -o paper.pdf --filter=/usr/local/bin/pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve
	# cd $(PAPER); Rscript -e "rmarkdown::render('paper.md')"
	
summary: $(PAPER)/summary.pdf
$(PAPER)/summary.pdf: $(PAPER)/header.yaml $(PAPER)/summary.md \
                    $(PAPER)/Young.bib \
                    $(wildcard $(PAPER)/fig_*.png)
	cd $(PAPER); pandoc header.yaml summary.md -o summary.pdf --filter=pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve
	# cd $(PAPER); Rscript -e "rmarkdown::render('summary.md')"
	
title: $(PAPER)/title.pdf
$(PAPER)/title.pdf: $(PAPER)/title.md
	cd $(PAPER); pandoc title.md -o title.pdf
	
diff: $(PAPER)/diff.pdf
$(PAPER)/diff.pdf: $(PAPER)/paper_20201211.md $(PAPER)/paper.md
	pandiff $(PAPER)/paper_20201211.md $(PAPER)/paper.md -s -o $@ --pdf-engine=lualatex


tex: $(PAPER)/summary.tex
$(PAPER)/summary.tex: summary
	cd $(PAPER); pandoc header.yaml summary.md -o summary.tex  --filter=/usr/local/bin/pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve

zip: $(PAPER)/paper.zip
$(PAPER)/paper.zip: $(PAPER)/paper.pdf \
                    $(PAPER)/summary.tex \
                    $(wildcard $(PAPER)/fig_*.png) \
                    $(PAPER)/summary.pdf \
                    $(PAPER)/Young.bib
	zip $@ $^


docx: $(PAPER)/summary.docx
$(PAPER)/summary.docx: $(PAPER)/header.yaml $(PAPER)/summary.md $(PAPER)/Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	cd $(PAPER); pandoc header.yaml summary.md -o summary.docx --filter=/usr/local/bin/pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve
	
wc: 
	pdftotext $(PAPER)/paper.pdf - | wc -w



