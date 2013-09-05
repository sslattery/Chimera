# Name of manuscript
manuscript = slattery_evans_wilson

# PdfLaTeX compilation options
latexopt = -halt-on-error -file-line-error

#=================================================================
# Generate PDF of manuscript using PdfLaTeX
#=================================================================

all:  domain_to_domain.pdf binary_comm_tree.pdf msod_example.pdf $(manuscript).pdf

$(manuscript).pdf: $(manuscript).tex references.bib
	pdflatex $(latexopt) $(manuscript)
	bibtex -terse $(manuscript)
	pdflatex $(latexopt) $(manuscript)
	pdflatex $(latexopt) $(manuscript)

domain_to_domain.pdf: domain_to_domain.fig
	fig2dev -L pdftex -p domain_to_domain.pdf -m 0.4 domain_to_domain.fig domain_to_domain.pdf
	fig2dev -L pdftex_t -p domain_to_domain.pdf -m 0.4 domain_to_domain.fig domain_to_domain.pdftex_t

binary_comm_tree.pdf: binary_comm_tree.fig
	fig2dev -L pdftex -p binary_comm_tree.pdf -m 0.4 binary_comm_tree.fig binary_comm_tree.pdf
	fig2dev -L pdftex_t -p binary_comm_tree.pdf -m 0.4 binary_comm_tree.fig binary_comm_tree.pdftex_t

msod_example.pdf: msod_example.fig
	fig2dev -L pdftex -p msod_example.pdf -m 0.4 msod_example.fig msod_example.pdf
	fig2dev -L pdftex_t -p msod_example.pdf -m 0.4 msod_example.fig msod_example.pdftex_t

#=================================================================
# Other
#=================================================================

clean:
	@rm -f *.aux *.bbl *.blg *.log *.out *.spl $(manuscript).pdf

.PHONY: all clean
