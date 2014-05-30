## Process this file with automake to produce Makefile.in

# LaTeX and BiBTeX source files
TEXSOURCES  = skewadjoint.tex
EXTRA_DIST  = references.bib
EXTRA_DIST += workbooks   # Entire directory contents

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
