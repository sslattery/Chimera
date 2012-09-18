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

all: orthogonal_residual.pdf partitioned_matrix.pdf matvec_proc_1.pdf
	$(PDFTEX) dissertation
	$(BIB) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation

# Linear problem figures
orthogonal_residual.pdf: $(chap_lp)/orthogonal_residual.fig
	fig2dev -L pdftex -p $(chap_lp)/orthogonal_residual.pdf -m 0.4 $(chap_lp)/orthogonal_residual.fig $(chap_lp)/orthogonal_residual.pdf
	fig2dev -L pdftex_t -p $(chap_lp)/orthogonal_residual.pdf -m 0.4 $(chap_lp)/orthogonal_residual.fig $(chap_lp)/orthogonal_residual.pdftex_t

partitioned_matrix.pdf: $(chap_lp)/partitioned_matrix.fig
	fig2dev -L pdftex -p $(chap_lp)/partitioned_matrix.pdf -m 0.4 $(chap_lp)/partitioned_matrix.fig $(chap_lp)/partitioned_matrix.pdf
	fig2dev -L pdftex_t -p $(chap_lp)/partitioned_matrix.pdf -m 0.4 $(chap_lp)/partitioned_matrix.fig $(chap_lp)/partitioned_matrix.pdftex_t

matvec_proc_1.pdf: $(chap_lp)/matvec_proc_1.fig
	fig2dev -L pdftex -p $(chap_lp)/matvec_proc_1.pdf -m 0.4 $(chap_lp)/matvec_proc_1.fig $(chap_lp)/matvec_proc_1.pdf
	fig2dev -L pdftex_t -p $(chap_lp)/matvec_proc_1.pdf -m 0.4 $(chap_lp)/matvec_proc_1.fig $(chap_lp)/matvec_proc_1.pdftex_t

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
	-rm -f $(chap_lp)/*.pdf $(chap_lp)/*.pdftex_t $(chap_lp)/*.bak
	-rm -f TAGS

realclean: clean
	-rm -f dissertation.pdf
