TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all:
	$(PDFTEX) sna_abstract
	$(BIB) sna_abstract
	$(PDFTEX) sna_abstract
	$(PDFTEX) sna_abstract
	$(PDFTEX) sna_abstract

clean:
	-rm -f *~
	-rm -f *.log
	-rm -f *.aux
	-rm -f *.dvi
	-rm -f *.rel
	-rm -f *.bbl
	-rm -f *.blg
	-rm -f *.out
	-rm -f *.toc *.lof *.lot
	-rm -f *.bak
	-rm -f *.pdftex_t
	-rm -f TAGS
	-rm -f *.pdf
