PAPER := paper
OUT := out
SCRIPT := scripts

all: script supplement paper zip

script:
	$(MAKE) -C $(SCRIPT)


supplement: $(PAPER)/supplement.pdf
$(PAPER)/supplement.pdf: $(PAPER)/header.yaml $(PAPER)/supplement.md $(PAPER)/Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	cd $(PAPER); pandoc header.yaml supplement.md -o supplement.pdf --filter=pandoc-crossref --filter=pandoc-citeproc --pdf-engine=lualatex --wrap=preserve


## Figures
$(PAPER)/fig_%.png: $(OUT)/fig_%.png
	cp -a $^ $(PAPER)


paper: $(PAPER)/paper.pdf
$(PAPER)/paper.pdf: $(PAPER)/header.yaml $(PAPER)/paper.md $(PAPER)/Young.bib \
                    $(wildcard $(PAPER)/fig_*.png)
	# cd $(PAPER); pandoc header.yaml paper.md -o paper.pdf --filter=pandoc-crossref --filter=pandoc-citeproc --pdf-engine=lualatex --wrap=preserve
	cd $(PAPER); Rscript -e "rmarkdown::render('paper.md')"


zip: $(PAPER)/paper.zip
$(PAPER)/paper.zip: $(PAPER)/paper.pdf \
                    $(PAPER)/paper.tex \
                    $(wildcard $(PAPER)/fig_*.png)
	zip $@ $^


docx: $(PAPER)/paper.docx
$(PAPER)/paper.docx: $(PAPER)/header.yaml $(PAPER)/paper.md $(PAPER)/Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	cd $(PAPER); pandoc header.yaml paper.md -o paper.docx --filter=pandoc-crossref --filter=pandoc-citeproc --pdf-engine=lualatex --wrap=preserve
	
wc: 
	pdftotext $(PAPER)/paper.pdf - | wc -w



