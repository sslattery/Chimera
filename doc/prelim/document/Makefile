TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

# sample.pdf: sample.fig 
# 	fig2dev -L pdftex -p sample.pdf -m 0.4 sample.fig sample.pdf
# 	fig2dev -L pdftex_t -p sample.pdf -m 0.4 sample.fig sample.pdftex_t

all: 
#	sample.pdf
	$(PDFTEX) dissertation
	$(BIB) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation

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
#	-rm -f sample.pdf 3D_sample.pdftex_t
	-rm -f TAGS

realclean: clean
	-rm -f dissertation.pdf
