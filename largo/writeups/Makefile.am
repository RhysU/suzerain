## Process this file with automake to produce Makefile.in

# Building writeups in parallel tends to break LaTeX-related scripting
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

# LaTeX and BiBTeX source files
TEXSOURCES  = skewadjoint.tex
EXTRA_DIST  = references.bib

# Build and install writeups
all-local:      pdf-local
pdf-local:      $(PDFS)
PDFS          = $(TEXSOURCES:.tex=.pdf)
writeupsdir   = $(pkgdatadir)/writeups
writeups_DATA = $(PDFS)

# Files to clean up based on common TeX-based tools
CLEANFILES = *.aux *.bbl *.blg *.dvi *.log *.pdf *.ps

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
