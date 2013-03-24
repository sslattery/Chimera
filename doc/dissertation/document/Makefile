TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

chap_intro = chapters/introduction
chap_lp = chapters/linear_problem
chap_mcb = chapters/mc_background
chap_pmc = chapters/parallel_mc
chap_mct = chapters/mc_theory
chap_spne = chapters/spn_equations
chap_spns = chapters/spn_solutions
chap_nl = chapters/nonlinear_problem
chap_ns = chapters/navier_stokes_solutions
chap_conc = chapters/conclusion

all: dnb_example.pdf orthogonal_residual.pdf partitioned_matrix.pdf matvec_proc_1.pdf gentile_example.pdf msod_example.pdf procassini_example.pdf stencil.pdf natural_convection_cavity.pdf lid_driven_cavity.pdf backward_facing_step.pdf
	$(PDFTEX) dissertation
	$(BIB) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation

# Introduction figures
dnb_example.pdf: $(chap_intro)/dnb_example.fig
	fig2dev -L pdftex -p $(chap_intro)/dnb_example.pdf -m 0.4 $(chap_intro)/dnb_example.fig $(chap_intro)/dnb_example.pdf
	fig2dev -L pdftex_t -p $(chap_intro)/dnb_example.pdf -m 0.4 $(chap_intro)/dnb_example.fig $(chap_intro)/dnb_example.pdftex_t

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

# Parallel mc figures
gentile_example.pdf: $(chap_pmc)/gentile_example.fig
	fig2dev -L pdftex -p $(chap_pmc)/gentile_example.pdf -m 0.4 $(chap_pmc)/gentile_example.fig $(chap_pmc)/gentile_example.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/gentile_example.pdf -m 0.4 $(chap_pmc)/gentile_example.fig $(chap_pmc)/gentile_example.pdftex_t

msod_example.pdf: $(chap_pmc)/msod_example.fig
	fig2dev -L pdftex -p $(chap_pmc)/msod_example.pdf -m 0.4 $(chap_pmc)/msod_example.fig $(chap_pmc)/msod_example.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/msod_example.pdf -m 0.4 $(chap_pmc)/msod_example.fig $(chap_pmc)/msod_example.pdftex_t

procassini_example.pdf: $(chap_pmc)/procassini_example.fig
	fig2dev -L pdftex -p $(chap_pmc)/procassini_example.pdf -m 0.4 $(chap_pmc)/procassini_example.fig $(chap_pmc)/procassini_example.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/procassini_example.pdf -m 0.4 $(chap_pmc)/procassini_example.fig $(chap_pmc)/procassini_example.pdftex_t

stencil.pdf: $(chap_mct)/stencil.fig
	fig2dev -L pdftex -p $(chap_mct)/stencil.pdf -m 0.4 $(chap_mct)/stencil.fig $(chap_mct)/stencil.pdf
	fig2dev -L pdftex_t -p $(chap_mct)/stencil.pdf -m 0.4 $(chap_mct)/stencil.fig $(chap_mct)/stencil.pdftex_t

# Navier stokes figures
natural_convection_cavity.pdf: $(chap_ns)/natural_convection_cavity.fig
	fig2dev -L pdftex -p $(chap_ns)/natural_convection_cavity.pdf -m 0.4 $(chap_ns)/natural_convection_cavity.fig $(chap_ns)/natural_convection_cavity.pdf
	fig2dev -L pdftex_t -p $(chap_ns)/natural_convection_cavity.pdf -m 0.4 $(chap_ns)/natural_convection_cavity.fig $(chap_ns)/natural_convection_cavity.pdftex_t

lid_driven_cavity.pdf: $(chap_ns)/lid_driven_cavity.fig
	fig2dev -L pdftex -p $(chap_ns)/lid_driven_cavity.pdf -m 0.4 $(chap_ns)/lid_driven_cavity.fig $(chap_ns)/lid_driven_cavity.pdf
	fig2dev -L pdftex_t -p $(chap_ns)/lid_driven_cavity.pdf -m 0.4 $(chap_ns)/lid_driven_cavity.fig $(chap_ns)/lid_driven_cavity.pdftex_t

backward_facing_step.pdf: $(chap_ns)/backward_facing_step.fig
	fig2dev -L pdftex -p $(chap_ns)/backward_facing_step.pdf -m 0.4 $(chap_ns)/backward_facing_step.fig $(chap_ns)/backward_facing_step.pdf
	fig2dev -L pdftex_t -p $(chap_ns)/backward_facing_step.pdf -m 0.4 $(chap_ns)/backward_facing_step.fig $(chap_ns)/backward_facing_step.pdftex_t


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
	-rm -f $(chap_intro)/*.pdftex_t $(chap_intro)/*.bak
	-rm -f $(chap_lp)/*.pdftex_t $(chap_lp)/*.bak
	-rm -f $(chap_mcb)/*.pdftex_t $(chap_mcb)/*.bak
	-rm -f $(chap_mct)/*.pdftex_t $(chap_mct)/*.bak
	-rm -f $(chap_pmc)/*.pdftex_t $(chap_pmc)/*.bak
	-rm -f $(chap_spne)/*.pdftex_t $(chap_spne)/*.bak
	-rm -f $(chap_spns)/*.pdftex_t $(chap_spns)/*.bak
	-rm -f $(chap_nl)/*.pdftex_t $(chap_nl)/*.bak
	-rm -f $(chap_ns)/*.pdftex_t $(chap_ns)/*.bak
	-rm -f $(chap_conc)/*.pdftex_t $(chap_conc)/*.bak
	-rm -f TAGS

realclean: clean
	-rm -f dissertation.pdf
