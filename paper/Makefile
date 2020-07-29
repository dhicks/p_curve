OUT := ../out

all: supplement paper wc

supplement: supplement.pdf
supplement.pdf: header.yaml supplement.md Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	pandoc header.yaml supplement.md -o supplement.pdf --filter=pandoc-crossref --filter=pandoc-citeproc --pdf-engine=lualatex --wrap=preserve

paper: paper.pdf
paper.pdf: header.yaml paper.md Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	pandoc header.yaml paper.md -o paper.pdf --filter=pandoc-crossref --filter=pandoc-citeproc --pdf-engine=lualatex --wrap=preserve
	
wc: 
	pdftotext paper.pdf - | wc -w


