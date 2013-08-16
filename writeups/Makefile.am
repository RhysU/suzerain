## Process this file with automake to produce Makefile.in

# Building writeups in parallel tends to break LaTeX-related scripting
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

# LaTeX and BiBTeX source files

# Standalone LaTeX documents
TEXSOURCES  = baseflow.tex
TEXSOURCES += channel_treatment.tex
TEXSOURCES += covariance.tex
TEXSOURCES += L2.tex
TEXSOURCES += lssdimex.tex
TEXSOURCES += perfectgas.tex
TEXSOURCES += plate_treatment.tex
TEXSOURCES += skewadjoint.tex
TEXSOURCES += smr91.tex

# Reacting doc broken out amongst many srcs.
reacting_TEXSOURCES  = reacting_main.tex
reacting_TEXSOURCES += reacting_inputs/boundary_conditions.tex
reacting_TEXSOURCES += reacting_inputs/discretization.tex
reacting_TEXSOURCES += reacting_inputs/filter.tex
reacting_TEXSOURCES += reacting_inputs/implicit_operator.tex
reacting_TEXSOURCES += reacting_inputs/reacting_flow_eqns.tex
reacting_TEXSOURCES += reacting_inputs/reacting_linearization.tex
reacting_TEXSOURCES += reacting_inputs/statistics.tex
reacting_TEXSOURCES += reacting_inputs/transfer_function.eps

# Explicit dependencies where necessary
reacting_main.dvi reacting_main.pdf: $(reacting_TEXSOURCES)
perfectgas.dvi perfectgas.pdf: perfectgas.tex
perfectgas.dvi perfectgas.pdf: inhomogeneity1.eps inhomogeneity2.eps
baseflow.dvi   baseflow.pdf:   nozzle_subsonic.eps nozzle_supersonic.eps

# Notebooks, Images, Source Data, References, etc.
EXTRA_DIST  = notebooks/Acoustics_WallNormal_3D.nb
EXTRA_DIST += notebooks/Acoustics_WallNormal_3D_reduced.nb
EXTRA_DIST += notebooks/Acoustics_WallNormal_MassFraction_3D.nb
EXTRA_DIST += notebooks/Acoustic_WallNormal_EOS_3d.nb
EXTRA_DIST += notebooks/Bump_Function_Powers.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_1D_Pressure.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_3D.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_3D_PerfectGas.nb
EXTRA_DIST += notebooks/Euler_Eigensystem_MassFraction_3D.nb
EXTRA_DIST += notebooks/FPG_BaseFlow_Problem.nb
EXTRA_DIST += notebooks/Giles_BC_Dimensional.nb
EXTRA_DIST += notebooks/Giles_BC_Nondimensional.nb
EXTRA_DIST += notebooks/nozzle1.m
EXTRA_DIST += notebooks/nozzle_baseflow.m
EXTRA_DIST += notebooks/nozzle_edge.m
EXTRA_DIST += notebooks/nozzle.m
EXTRA_DIST += notebooks/SMR91_Linearization.nb
EXTRA_DIST += notebooks/Verification_Giles_Reacting.nb
EXTRA_DIST += nozzle_schematic.eps nozzle_subsonic.eps nozzle_supersonic.eps
EXTRA_DIST += inhomogeneity1.gz inhomogeneity1.eps
EXTRA_DIST += inhomogeneity2.gz inhomogeneity2.eps
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
