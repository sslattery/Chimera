TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

chap_intro = chapters/introduction
chap_lp = chapters/linear_problem
chap_pmcsa = chapters/parallel_mcsa
chap_nl = chapters/nonlinear_problem
chap_fanm = chapters/fanm
chap_chm = chapters/chimera_codebase
chap_prop = chapters/research_proposal

all:
	$(PDFTEX) dissertation
	$(BIB) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation

orthogonal_residual.pdf:
	fig2dev -L pdftex -p $(chap_lp)/orthogonal_residual.pdf -m 0.4 $(chap_lp)/orthogonal_residual.fig $(chap_lp)/orthogonal_residual.pdf
	fig2dev -L pdftex_t -p $(chap_lp)/orthogonal_residual.pdf -m 0.4 $(chap_lp)/orthogonal_residual.fig $(chap_lp)/orthogonal_residual.pdftex_t

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
	-rm -f $(chap_lp)/orthogonal_residual.pdf $(chap_lp)/orthogonal_residual.pdftex_t
	-rm -f TAGS

realclean: clean
	-rm -f dissertation.pdf
