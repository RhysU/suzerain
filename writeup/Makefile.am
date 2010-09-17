## Process this file with automake to produce Makefile.in

# Building writeups in parallel tends to break LaTeX-related scripting
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

TEXSOURCES  = derivation.tex
TEXSOURCES += skewadjoint.tex
TEXSOURCES += smr91.tex
TEXSOURCES += timestepping.tex

EXTRA_DIST = derivation.bib skewadjoint.bib

all-local: dvi-local pdf-local
dvi-local: $(TEXSOURCES:.tex=.dvi)
pdf-local: $(TEXSOURCES:.tex=.pdf)

CLEANFILES = *.aux *.bbl *.blg *.dvi *.log *.pdf

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
