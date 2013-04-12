## Process this file with automake to produce Makefile.in

# Building writeups in parallel tends to break LaTeX-related scripting
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

# LaTeX and BiBTeX source files

# Standalone LaTeX documents
TEXSOURCES  = channel_treatment.tex
TEXSOURCES += covariance.tex
TEXSOURCES += L2.tex
TEXSOURCES += lssdimex.tex
TEXSOURCES += perfectgas.tex
TEXSOURCES += skewadjoint.tex
TEXSOURCES += smr91.tex

# Reacting doc broken out amongst many srcs.
reacting_TEXSOURCES  = reacting_main.tex
reacting_TEXSOURCES += reacting_inputs/acoustic_decision_pt.tex
reacting_TEXSOURCES += reacting_inputs/acoustic_linearization.tex
reacting_TEXSOURCES += reacting_inputs/filter.tex
reacting_TEXSOURCES += reacting_inputs/full_linear_operator.tex
reacting_TEXSOURCES += reacting_inputs/reactinggas.tex
reacting_TEXSOURCES += reacting_inputs/reacting_linearization.tex
reacting_TEXSOURCES += reacting_inputs/viscous_linearization.tex
reacting_TEXSOURCES += reacting_inputs/transfer_function.eps

# Explicit dependencies to get reacting_main to rebuild when necessary
reacting_main.pdf: $(reacting_TEXSOURCES)

# Notebooks, Images, References, etc.
EXTRA_DIST  = notebooks/Acoustics_WallNormal_3D.nb
EXTRA_DIST += notebooks/Acoustics_WallNormal_3D_reduced.nb
EXTRA_DIST += notebooks/Acoustics_WallNormal_MassFraction_3D.nb
EXTRA_DIST += notebooks/Acoustic_WallNormal_EOS_3d.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_1D_Pressure.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_3D.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_MassFraction_3D.nb
EXTRA_DIST += notebooks/Giles_BC_Dimensional.nb
EXTRA_DIST += notebooks/Giles_BC_Nondimensional.nb
EXTRA_DIST += notebooks/SMR91_Linearization.nb
EXTRA_DIST += inhomogeneity.eps
EXTRA_DIST += inhomogeneity.gz
EXTRA_DIST += references.bib


# Build and install writeups
all-local:      pdf-local
pdf-local:      $(PDFS)
PDFS          = $(TEXSOURCES:.tex=.pdf) reacting_main.pdf
writeupsdir   = $(pkgdatadir)/writeups
writeups_DATA = $(PDFS)

# Files to clean up based on common TeX-based tools
CLEANFILES  = *.aux *.bbl *.blg *.dvi *.log *.out *.pdf *.ps *.toc
CLEANFILES += *-generated.pdf

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
