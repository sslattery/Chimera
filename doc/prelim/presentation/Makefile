TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: dnb_example.pdf orthogonal_residual.pdf partitioned_matrix.pdf matvec_proc_1.pdf gentile_example.pdf msod_example.pdf procassini_example.pdf heat_eq_setup.pdf natural_convection_cavity.pdf lid_driven_cavity.pdf backward_facing_step.pdf
	$(PDFTEX) prelim_presentation
	$(BIB) prelim_presentation
	$(PDFTEX) prelim_presentation
	$(PDFTEX) prelim_presentation
	$(PDFTEX) prelim_presentation

# Introduction figures
dnb_example.pdf: dnb_example.fig
	fig2dev -L pdftex -p dnb_example.pdf -m 0.4 dnb_example.fig dnb_example.pdf
	fig2dev -L pdftex_t -p dnb_example.pdf -m 0.4 dnb_example.fig dnb_example.pdftex_t

# Linear problem figures
orthogonal_residual.pdf: orthogonal_residual.fig
	fig2dev -L pdftex -p orthogonal_residual.pdf -m 0.4 orthogonal_residual.fig orthogonal_residual.pdf
	fig2dev -L pdftex_t -p orthogonal_residual.pdf -m 0.4 orthogonal_residual.fig orthogonal_residual.pdftex_t

partitioned_matrix.pdf: partitioned_matrix.fig
	fig2dev -L pdftex -p partitioned_matrix.pdf -m 0.4 partitioned_matrix.fig partitioned_matrix.pdf
	fig2dev -L pdftex_t -p partitioned_matrix.pdf -m 0.4 partitioned_matrix.fig partitioned_matrix.pdftex_t

matvec_proc_1.pdf: matvec_proc_1.fig
	fig2dev -L pdftex -p matvec_proc_1.pdf -m 0.4 matvec_proc_1.fig matvec_proc_1.pdf
	fig2dev -L pdftex_t -p matvec_proc_1.pdf -m 0.4 matvec_proc_1.fig matvec_proc_1.pdftex_t

# Parallel mc figures
gentile_example.pdf: gentile_example.fig
	fig2dev -L pdftex -p gentile_example.pdf -m 0.4 gentile_example.fig gentile_example.pdf
	fig2dev -L pdftex_t -p gentile_example.pdf -m 0.4 gentile_example.fig gentile_example.pdftex_t

msod_example.pdf: msod_example.fig
	fig2dev -L pdftex -p msod_example.pdf -m 0.4 msod_example.fig msod_example.pdf
	fig2dev -L pdftex_t -p msod_example.pdf -m 0.4 msod_example.fig msod_example.pdftex_t

procassini_example.pdf: procassini_example.fig
	fig2dev -L pdftex -p procassini_example.pdf -m 0.4 procassini_example.fig procassini_example.pdf
	fig2dev -L pdftex_t -p procassini_example.pdf -m 0.4 procassini_example.fig procassini_example.pdftex_t

# Research proposal figures
heat_eq_setup.pdf: heat_eq_setup.fig
	fig2dev -L pdftex -p heat_eq_setup.pdf -m 0.4 heat_eq_setup.fig heat_eq_setup.pdf
	fig2dev -L pdftex_t -p heat_eq_setup.pdf -m 0.4 heat_eq_setup.fig heat_eq_setup.pdftex_t

natural_convection_cavity.pdf: natural_convection_cavity.fig
	fig2dev -L pdftex -p natural_convection_cavity.pdf -m 0.4 natural_convection_cavity.fig natural_convection_cavity.pdf
	fig2dev -L pdftex_t -p natural_convection_cavity.pdf -m 0.4 natural_convection_cavity.fig natural_convection_cavity.pdftex_t

lid_driven_cavity.pdf: lid_driven_cavity.fig
	fig2dev -L pdftex -p lid_driven_cavity.pdf -m 0.4 lid_driven_cavity.fig lid_driven_cavity.pdf
	fig2dev -L pdftex_t -p lid_driven_cavity.pdf -m 0.4 lid_driven_cavity.fig lid_driven_cavity.pdftex_t

backward_facing_step.pdf: backward_facing_step.fig
	fig2dev -L pdftex -p backward_facing_step.pdf -m 0.4 backward_facing_step.fig backward_facing_step.pdf
	fig2dev -L pdftex_t -p backward_facing_step.pdf -m 0.4 backward_facing_step.fig backward_facing_step.pdftex_t


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
	-rm -f *.pdf
	-rm -f *.pdftex_t
	-rm -f TAGS
