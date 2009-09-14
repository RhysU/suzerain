## Process this file with automake to produce Makefile.in

TEXSOURCES = derivation.tex timestepping.tex smr91.tex skewadjoint.tex

all-local: dvi-local pdf-local
dvi-local: $(TEXSOURCES:.tex=.dvi)
pdf-local: $(TEXSOURCES:.tex=.pdf)

CLEANFILES = *.aux *.bbl *.blg *.dvi *.log *.pdf

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
