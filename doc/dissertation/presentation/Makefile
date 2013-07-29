TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: heat_eq_setup.pdf ddmc_example.pdf ddnu_example.pdf msod_example.pdf domain_to_domain.pdf binary_comm_tree.pdf msod_construction.pdf msod_tally.pdf
	$(PDFTEX) defense.tex

heat_eq_setup.pdf: heat_eq_setup.fig
	fig2dev -L pdftex -p heat_eq_setup.pdf -m 0.4 heat_eq_setup.fig heat_eq_setup.pdf
	fig2dev -L pdftex_t -p heat_eq_setup.pdf -m 0.4 heat_eq_setup.fig heat_eq_setup.pdftex_t

ddmc_example.pdf: ddmc_example.fig
	fig2dev -L pdftex -p ddmc_example.pdf -m 0.4 ddmc_example.fig ddmc_example.pdf
	fig2dev -L pdftex_t -p ddmc_example.pdf -m 0.4 ddmc_example.fig ddmc_example.pdftex_t

ddnu_example.pdf: ddnu_example.fig
	fig2dev -L pdftex -p ddnu_example.pdf -m 0.4 ddnu_example.fig ddnu_example.pdf
	fig2dev -L pdftex_t -p ddnu_example.pdf -m 0.4 ddnu_example.fig ddnu_example.pdftex_t

msod_example.pdf: msod_example.fig
	fig2dev -L pdftex -p msod_example.pdf -m 0.4 msod_example.fig msod_example.pdf
	fig2dev -L pdftex_t -p msod_example.pdf -m 0.4 msod_example.fig msod_example.pdftex_t

domain_to_domain.pdf: domain_to_domain.fig
	fig2dev -L pdftex -p domain_to_domain.pdf -m 0.4 domain_to_domain.fig domain_to_domain.pdf
	fig2dev -L pdftex_t -p domain_to_domain.pdf -m 0.4 domain_to_domain.fig domain_to_domain.pdftex_t

binary_comm_tree.pdf: binary_comm_tree.fig
	fig2dev -L pdftex -p binary_comm_tree.pdf -m 0.4 binary_comm_tree.fig binary_comm_tree.pdf
	fig2dev -L pdftex_t -p binary_comm_tree.pdf -m 0.4 binary_comm_tree.fig binary_comm_tree.pdftex_t

msod_construction.pdf: msod_construction.fig
	fig2dev -L pdftex -p msod_construction.pdf -m 0.4 msod_construction.fig msod_construction.pdf
	fig2dev -L pdftex_t -p msod_construction.pdf -m 0.4 msod_construction.fig msod_construction.pdftex_t

msod_tally.pdf: msod_tally.fig
	fig2dev -L pdftex -p msod_tally.pdf -m 0.4 msod_tally.fig msod_tally.pdf
	fig2dev -L pdftex_t -p msod_tally.pdf -m 0.4 msod_tally.fig msod_tally.pdftex_t
