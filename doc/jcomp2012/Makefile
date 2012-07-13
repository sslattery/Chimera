TEX    := latex
BIB    := bibtex
PDFTEX := pdflatex
PDF    := dvipdf

mesh_cell.pdf: mesh_cell.fig 
	fig2dev -L pdftex -p mesh_cell.pdf -m 0.4 mesh_cell.fig mesh_cell.pdf
	fig2dev -L pdftex_t -p mesh_cell.pdf -m 0.4 mesh_cell.fig mesh_cell.pdftex_t

3D_mesh_segment.pdf: 3D_mesh_segment.fig 
	fig2dev -L pdftex -p 3D_mesh_segment.pdf -m 0.4 3D_mesh_segment.fig 3D_mesh_segment.pdf
	fig2dev -L pdftex_t -p 3D_mesh_segment.pdf -m 0.4 3D_mesh_segment.fig 3D_mesh_segment.pdftex_t

all: 
	$(TEX) jcomp2012
	$(BIB) jcomp2012
	$(TEX) jcomp2012
	$(TEX) jcomp2012
	$(PDF) jcomp2012
