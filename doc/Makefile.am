## Process this file with automake to produce Makefile.in

TEXSOURCES = derivation.tex

all-local: dvi-local pdf-local
dvi-local: $(TEXSOURCES:.tex=.dvi)
pdf-local: $(TEXSOURCES:.tex=.pdf)

CLEANFILES = *.aux *.dvi *.log *.pdf

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
