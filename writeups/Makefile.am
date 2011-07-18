## Process this file with automake to produce Makefile.in

# Building writeups in parallel tends to break LaTeX-related scripting
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

# LaTeX and BiBTeX source files
TEXSOURCES  = channel_treatment.tex
TEXSOURCES += derivation.tex
TEXSOURCES += L2.tex
TEXSOURCES += skewadjoint.tex
TEXSOURCES += smr91.tex
EXTRA_DIST  = references.bib

# Build and install documentation PDFs
PDFS     = $(TEXSOURCES:.tex=.pdf)
doc_DATA = $(PDFS)
all-local:         pdf-local
pdf-local:         $(PDFS)
install-pdf-local: pdf-local

# Files to clean up based on common TeX-based tools
CLEANFILES = *.aux *.bbl *.blg *.dvi *.log *.pdf *.ps

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
