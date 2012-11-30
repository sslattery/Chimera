TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: stencil.pdf
	$(PDFTEX) slattery_wilson_evans
	$(BIB) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans
	$(PDFTEX) slattery_wilson_evans

stencil.pdf: stencil.fig
	fig2dev -L pdftex -p stencil.pdf -m 0.4 stencil.fig stencil.pdf
	fig2dev -L pdftex_t -p stencil.pdf -m 0.4 stencil.fig stencil.pdftex_t

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
