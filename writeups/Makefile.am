## Process this file with automake to produce Makefile.in

# Building writeups in parallel tends to break LaTeX-related scripting
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

# LaTeX and BiBTeX source files
TEXSOURCES  = channel_treatment.tex
TEXSOURCES += covariance.tex
TEXSOURCES += filter.tex
TEXSOURCES += L2.tex
TEXSOURCES += lssdimex.tex
TEXSOURCES += perfectgas.tex
TEXSOURCES += skewadjoint.tex
TEXSOURCES += smr91.tex
TEXSOURCES += reactinggas.tex
TEXSOURCES += reacting_linearization.tex

# Images, References, etc.
EXTRA_DIST  = inhomogeneity.gz
EXTRA_DIST += inhomogeneity.eps
EXTRA_DIST += references.bib
EXTRA_DIST += transfer_function.eps

# Build and install writeups
all-local:      pdf-local
pdf-local:      $(PDFS)
PDFS          = $(TEXSOURCES:.tex=.pdf)
writeupsdir   = $(pkgdatadir)/writeups
writeups_DATA = $(PDFS)

# Files to clean up based on common TeX-based tools
CLEANFILES  = *.aux *.bbl *.blg *.dvi *.log *.out *.pdf *.ps *.toc
CLEANFILES += *-generated.pdf

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
