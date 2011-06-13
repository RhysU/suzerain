## Process this file with automake to produce Makefile.in

# Building writeups in parallel tends to break LaTeX-related scripting
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

TEXSOURCES  = channel_treatment.tex
TEXSOURCES += derivation.tex
TEXSOURCES += L2.tex
TEXSOURCES += skewadjoint.tex
TEXSOURCES += smr91.tex

EXTRA_DIST = references.bib

all-local: dvi-local pdf-local
dvi-local: $(TEXSOURCES:.tex=.dvi)
pdf-local: $(TEXSOURCES:.tex=.pdf)

CLEANFILES = *.aux *.bbl *.blg *.dvi *.log *.pdf *.ps

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
