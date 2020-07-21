all: paper wc

paper: paper.pdf
paper.pdf: header.yaml paper.md Young.bib
	pandoc header.yaml paper.md -o paper.pdf --filter=pandoc-crossref --filter=pandoc-citeproc --pdf-engine=lualatex --wrap=preserve
	
wc: 
	pdftotext paper.pdf - | wc -w


