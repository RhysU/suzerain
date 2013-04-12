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

# Reacting doc broken out amongst many srcs.  Use EXTRA_DIST for
# everything that gets \input in reacting_main.tex s.t. we don't try
# to build them standalone (which will fail of course).
#
# FIXME: This approach breaks dependency tracking.  Thus, a mod to one
# of the files that is \input does *not* result in a rebuild of
# reacting_main.
TEXSOURCES += reacting_main.tex
EXTRA_DIST  = reacting_inputs/acoustic_decision_pt.tex
EXTRA_DIST += reacting_inputs/acoustic_linearization.tex
EXTRA_DIST += reacting_inputs/filter.tex
EXTRA_DIST += reacting_inputs/full_linear_operator.tex
EXTRA_DIST += reacting_inputs/reactinggas.tex
EXTRA_DIST += reacting_inputs/reacting_linearization.tex
EXTRA_DIST += reacting_inputs/viscous_linearization.tex
EXTRA_DIST += reacting_inputs/transfer_function.eps

# Notebooks, Images, References, etc.
EXTRA_DIST += notebooks/Acoustics_WallNormal_3D.nb
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
PDFS          = $(TEXSOURCES:.tex=.pdf)
writeupsdir   = $(pkgdatadir)/writeups
writeups_DATA = $(PDFS)

# Files to clean up based on common TeX-based tools
CLEANFILES  = *.aux *.bbl *.blg *.dvi *.log *.out *.pdf *.ps *.toc
CLEANFILES += *-generated.pdf

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
