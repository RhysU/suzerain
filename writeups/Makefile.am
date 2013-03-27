## Process this file with automake to produce Makefile.in

# Building writeups in parallel tends to break LaTeX-related scripting
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

# LaTeX and BiBTeX source files
TEXSOURCES  = acoustic_linearization.tex
TEXSOURCES += channel_treatment.tex
TEXSOURCES += covariance.tex
TEXSOURCES += filter.tex
TEXSOURCES += L2.tex
TEXSOURCES += lssdimex.tex
TEXSOURCES += perfectgas.tex
TEXSOURCES += reactinggas.tex
TEXSOURCES += reacting_linearization.tex
TEXSOURCES += skewadjoint.tex
TEXSOURCES += smr91.tex
TEXSOURCES += viscous_linearization.tex

# Notebooks, Images, References, etc.
EXTRA_DIST  = notebooks/Acoustics_WallNormal_3D.nb
EXTRA_DIST += notebooks/Acoustics_WallNormal_MassFraction_3D.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_MassFraction_3D.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_1D_Pressure.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_3D.nb
EXTRA_DIST += inhomogeneity.gz
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
