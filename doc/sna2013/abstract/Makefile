# Name of manuscript
manuscript = slattery_wilson_evans

#=================================================================
# Generate PDF of manuscript using PdfLaTeX
#=================================================================

all: $(manuscript).pdf

$(manuscript).pdf: $(manuscript).tex references.bib
	pdflatex $(manuscript)
	bibtex -terse $(manuscript)
	pdflatex $(manuscript)
	pdflatex $(manuscript)

#=================================================================
# Other
#=================================================================

clean:
	@rm -f *.aux *.bbl *.blg *.log *.out *.spl $(manuscript).pdf

.PHONY: all clean
