PAPER := paper
OUT := out
SCRIPT := scripts

all: script supplement paper zip diff

script:
	$(MAKE) -C $(SCRIPT)


supplement: $(PAPER)/supplement.pdf
$(PAPER)/supplement.pdf: $(PAPER)/header.yaml $(PAPER)/supplement.md $(PAPER)/Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	cd $(PAPER); pandoc header.yaml supplement.md -o supplement.pdf --filter=pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve


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
	cd $(PAPER); pandoc header.yaml paper.md -o paper.pdf --filter=pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve
	# cd $(PAPER); Rscript -e "rmarkdown::render('paper.md')"
	
diff: diff.pdf
diff.pdf: paper_20210208.md paper.md
	pandiff paper_20210208.md paper.md -s -o diff.pdf


tex: $(PAPER)/paper.tex
$(PAPER)/paper.tex: paper
	cd $(PAPER); pandoc header.yaml paper.md -o paper.tex  --filter=pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve

zip: $(PAPER)/paper.zip
$(PAPER)/paper.zip: $(PAPER)/paper.pdf \
                    $(PAPER)/paper.tex \
                    $(wildcard $(PAPER)/fig_*.png) \
                    $(PAPER)/supplement.pdf \
                    $(PAPER)/Young.bib
	zip $@ $^


docx: $(PAPER)/paper.docx
$(PAPER)/paper.docx: $(PAPER)/header.yaml $(PAPER)/paper.md $(PAPER)/Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	cd $(PAPER); pandoc header.yaml paper.md -o paper.docx --filter=pandoc-crossref --citeproc --pdf-engine=lualatex --wrap=preserve
	
wc: 
	pdftotext $(PAPER)/paper.pdf - | wc -w



