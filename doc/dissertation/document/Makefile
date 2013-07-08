TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

chap_intro = chapters/introduction
chap_lp = chapters/linear_problem
chap_mcb = chapters/mc_background
chap_pmc = chapters/parallel_mc
chap_spne = chapters/spn_equations
chap_nl = chapters/nonlinear_problem
chap_conc = chapters/conclusion
apdx = backmatter

all: dnb_example.pdf orthogonal_residual.pdf partitioned_matrix.pdf matvec_proc_1.pdf heat_eq_setup.pdf block_jacobi.pdf ddmc_example.pdf ddnu_example.pdf gentile_example.pdf msod_example.pdf domain_to_domain.pdf binary_comm_tree.pdf master_comm_tree.pdf msod_construction.pdf msod_tally.pdf procassini_example.pdf stencil.pdf stencil_graph.pdf natural_convection_cavity.pdf lid_driven_cavity.pdf backward_facing_step.pdf mesh_cell.pdf
	$(PDFTEX) dissertation
	$(BIB) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation
	$(PDFTEX) dissertation

tex:
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

# mc background
heat_eq_setup.pdf: $(chap_mcb)/heat_eq_setup.fig
	fig2dev -L pdftex -p $(chap_mcb)/heat_eq_setup.pdf -m 0.4 $(chap_mcb)/heat_eq_setup.fig $(chap_mcb)/heat_eq_setup.pdf
	fig2dev -L pdftex_t -p $(chap_mcb)/heat_eq_setup.pdf -m 0.4 $(chap_mcb)/heat_eq_setup.fig $(chap_mcb)/heat_eq_setup.pdftex_t

# spn figures
block_jacobi.pdf: $(chap_spne)/block_jacobi.fig
	fig2dev -L pdftex -p $(chap_spne)/block_jacobi.pdf -m 0.4 $(chap_spne)/block_jacobi.fig $(chap_spne)/block_jacobi.pdf
	fig2dev -L pdftex_t -p $(chap_spne)/block_jacobi.pdf -m 0.4 $(chap_spne)/block_jacobi.fig $(chap_spne)/block_jacobi.pdftex_t

stencil.pdf: $(chap_spne)/stencil.fig
	fig2dev -L pdftex -p $(chap_spne)/stencil.pdf -m 0.4 $(chap_spne)/stencil.fig $(chap_spne)/stencil.pdf
	fig2dev -L pdftex_t -p $(chap_spne)/stencil.pdf -m 0.4 $(chap_spne)/stencil.fig $(chap_spne)/stencil.pdftex_t

# Parallel mc figures
ddmc_example.pdf: $(chap_pmc)/ddmc_example.fig
	fig2dev -L pdftex -p $(chap_pmc)/ddmc_example.pdf -m 0.4 $(chap_pmc)/ddmc_example.fig $(chap_pmc)/ddmc_example.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/ddmc_example.pdf -m 0.4 $(chap_pmc)/ddmc_example.fig $(chap_pmc)/ddmc_example.pdftex_t

ddnu_example.pdf: $(chap_pmc)/ddnu_example.fig
	fig2dev -L pdftex -p $(chap_pmc)/ddnu_example.pdf -m 0.4 $(chap_pmc)/ddnu_example.fig $(chap_pmc)/ddnu_example.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/ddnu_example.pdf -m 0.4 $(chap_pmc)/ddnu_example.fig $(chap_pmc)/ddnu_example.pdftex_t

gentile_example.pdf: $(chap_pmc)/gentile_example.fig
	fig2dev -L pdftex -p $(chap_pmc)/gentile_example.pdf -m 0.4 $(chap_pmc)/gentile_example.fig $(chap_pmc)/gentile_example.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/gentile_example.pdf -m 0.4 $(chap_pmc)/gentile_example.fig $(chap_pmc)/gentile_example.pdftex_t

msod_example.pdf: $(chap_pmc)/msod_example.fig
	fig2dev -L pdftex -p $(chap_pmc)/msod_example.pdf -m 0.4 $(chap_pmc)/msod_example.fig $(chap_pmc)/msod_example.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/msod_example.pdf -m 0.4 $(chap_pmc)/msod_example.fig $(chap_pmc)/msod_example.pdftex_t

procassini_example.pdf: $(chap_pmc)/procassini_example.fig
	fig2dev -L pdftex -p $(chap_pmc)/procassini_example.pdf -m 0.4 $(chap_pmc)/procassini_example.fig $(chap_pmc)/procassini_example.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/procassini_example.pdf -m 0.4 $(chap_pmc)/procassini_example.fig $(chap_pmc)/procassini_example.pdftex_t

domain_to_domain.pdf: $(chap_pmc)/domain_to_domain.fig
	fig2dev -L pdftex -p $(chap_pmc)/domain_to_domain.pdf -m 0.4 $(chap_pmc)/domain_to_domain.fig $(chap_pmc)/domain_to_domain.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/domain_to_domain.pdf -m 0.4 $(chap_pmc)/domain_to_domain.fig $(chap_pmc)/domain_to_domain.pdftex_t

binary_comm_tree.pdf: $(chap_pmc)/binary_comm_tree.fig
	fig2dev -L pdftex -p $(chap_pmc)/binary_comm_tree.pdf -m 0.4 $(chap_pmc)/binary_comm_tree.fig $(chap_pmc)/binary_comm_tree.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/binary_comm_tree.pdf -m 0.4 $(chap_pmc)/binary_comm_tree.fig $(chap_pmc)/binary_comm_tree.pdftex_t

master_comm_tree.pdf: $(chap_pmc)/master_comm_tree.fig
	fig2dev -L pdftex -p $(chap_pmc)/master_comm_tree.pdf -m 0.4 $(chap_pmc)/master_comm_tree.fig $(chap_pmc)/master_comm_tree.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/master_comm_tree.pdf -m 0.4 $(chap_pmc)/master_comm_tree.fig $(chap_pmc)/master_comm_tree.pdftex_t

msod_construction.pdf: $(chap_pmc)/msod_construction.fig
	fig2dev -L pdftex -p $(chap_pmc)/msod_construction.pdf -m 0.4 $(chap_pmc)/msod_construction.fig $(chap_pmc)/msod_construction.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/msod_construction.pdf -m 0.4 $(chap_pmc)/msod_construction.fig $(chap_pmc)/msod_construction.pdftex_t

msod_tally.pdf: $(chap_pmc)/msod_tally.fig
	fig2dev -L pdftex -p $(chap_pmc)/msod_tally.pdf -m 0.4 $(chap_pmc)/msod_tally.fig $(chap_pmc)/msod_tally.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/msod_tally.pdf -m 0.4 $(chap_pmc)/msod_tally.fig $(chap_pmc)/msod_tally.pdftex_t

stencil_graph.pdf: $(chap_pmc)/stencil_graph.fig
	fig2dev -L pdftex -p $(chap_pmc)/stencil_graph.pdf -m 0.4 $(chap_pmc)/stencil_graph.fig $(chap_pmc)/stencil_graph.pdf
	fig2dev -L pdftex_t -p $(chap_pmc)/stencil_graph.pdf -m 0.4 $(chap_pmc)/stencil_graph.fig $(chap_pmc)/stencil_graph.pdftex_t

# Navier stokes figures
natural_convection_cavity.pdf: $(chap_nl)/natural_convection_cavity.fig
	fig2dev -L pdftex -p $(chap_nl)/natural_convection_cavity.pdf -m 0.4 $(chap_nl)/natural_convection_cavity.fig $(chap_nl)/natural_convection_cavity.pdf
	fig2dev -L pdftex_t -p $(chap_nl)/natural_convection_cavity.pdf -m 0.4 $(chap_nl)/natural_convection_cavity.fig $(chap_nl)/natural_convection_cavity.pdftex_t

lid_driven_cavity.pdf: $(chap_nl)/lid_driven_cavity.fig
	fig2dev -L pdftex -p $(chap_nl)/lid_driven_cavity.pdf -m 0.4 $(chap_nl)/lid_driven_cavity.fig $(chap_nl)/lid_driven_cavity.pdf
	fig2dev -L pdftex_t -p $(chap_nl)/lid_driven_cavity.pdf -m 0.4 $(chap_nl)/lid_driven_cavity.fig $(chap_nl)/lid_driven_cavity.pdftex_t

backward_facing_step.pdf: $(chap_nl)/backward_facing_step.fig
	fig2dev -L pdftex -p $(chap_nl)/backward_facing_step.pdf -m 0.4 $(chap_nl)/backward_facing_step.fig $(chap_nl)/backward_facing_step.pdf
	fig2dev -L pdftex_t -p $(chap_nl)/backward_facing_step.pdf -m 0.4 $(chap_nl)/backward_facing_step.fig $(chap_nl)/backward_facing_step.pdftex_t

# Backmatter
mesh_cell.pdf: $(apdx)/mesh_cell.fig
	fig2dev -L pdftex -p $(apdx)/mesh_cell.pdf -m 0.4 $(apdx)/mesh_cell.fig $(apdx)/mesh_cell.pdf
	fig2dev -L pdftex_t -p $(apdx)/mesh_cell.pdf -m 0.4 $(apdx)/mesh_cell.fig $(apdx)/mesh_cell.pdftex_t

#cleanup
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
	-rm -f $(chap_pmc)/*.pdftex_t $(chap_pmc)/*.bak
	-rm -f $(chap_spne)/*.pdftex_t $(chap_spne)/*.bak
	-rm -f $(chap_nl)/*.pdftex_t $(chap_nl)/*.bak
	-rm -f $(chap_nl)/*.pdftex_t $(chap_nl)/*.bak
	-rm -f $(chap_conc)/*.pdftex_t $(chap_conc)/*.bak
	-rm -f TAGS

realclean: clean
	-rm -f dissertation.pdf
