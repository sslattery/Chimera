TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: msod_example.pdf stencil.pdf
	$(PDFTEX) mandc_2013

msod_example.pdf: msod_example.fig
	fig2dev -L pdftex -p msod_example.pdf -m 0.4 msod_example.fig msod_example.pdf
	fig2dev -L pdftex_t -p msod_example.pdf -m 0.4 msod_example.fig msod_example.pdftex_t

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
	-rm -f *.nav
	-rm -f *.snm
	-rm -f TAGS
