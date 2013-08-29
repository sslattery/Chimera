# Name of manuscript
manuscript = slattery_evans_wilson

# PdfLaTeX compilation options
latexopt = -halt-on-error -file-line-error

# List of images to include
images = fig1.png fig2.png

#=================================================================
# Generate PDF of manuscript using PdfLaTeX
#=================================================================

all: $(manuscript).pdf

$(manuscript).pdf: $(manuscript).tex $(images) references.bib
	pdflatex $(latexopt) $(manuscript)
	bibtex -terse $(manuscript)
	pdflatex $(latexopt) $(manuscript)
	pdflatex $(latexopt) $(manuscript)

#=================================================================
# Other
#=================================================================

clean:
	@rm -f *.aux *.bbl *.blg *.log *.out *.spl $(manuscript).pdf

.PHONY: all clean
