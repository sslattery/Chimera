TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

all: jcomp2012.tex 3D_mesh_segment.pdf mesh_cell.pdf
	$(PDFTEX) jcomp2012
	$(BIB) jcomp2012
	$(PDFTEX) jcomp2012
	$(PDFTEX) jcomp2012

#mesh_cell.pdf: mesh_cell.fig 
#fig2dev -L pdftex -p mesh_cell.pdf -m 0.4 mesh_cell.fig mesh_cell.pdf
#fig2dev -L pdftex_t -p mesh_cell.pdf -m 0.4 mesh_cell.fig mesh_cell.tex

#3D_mesh_segment.pdf: 3D_mesh_segment.fig 
#fig2dev -L pdftex -p 3D_mesh_segment.pdf -m 0.4 3D_mesh_segment.fig 3D_mesh_segment.pdf
#fig2dev -L pdftex_t -p 3D_mesh_segment.pdf -m 0.4 3D_mesh_segment.fig 3D_mesh_segment.tex

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
	-rm -f TAGS

realclean: clean
	-rm -f jcomp2012.pdf
