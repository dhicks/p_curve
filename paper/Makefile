OUT := ../out

all: paper wc

paper: paper.pdf
paper.pdf: header.yaml paper.md Young.bib \
           $(wildcard $(OUT)/*.png) \
           $(wildcard $(OUT)/*.tex)
	pandoc header.yaml paper.md -o paper.pdf --filter=pandoc-crossref --filter=pandoc-citeproc --pdf-engine=lualatex --wrap=preserve
	
wc: 
	pdftotext paper.pdf - | wc -w


