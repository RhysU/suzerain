## Process this file with automake to produce Makefile.in

# LaTeX and BiBTeX source files

# Standalone LaTeX documents
TEXSOURCES  = covariance.tex
TEXSOURCES += L2.tex
TEXSOURCES += lssdimex.tex
TEXSOURCES += perfectgas.tex
TEXSOURCES += radialflow.tex
TEXSOURCES += skewadjoint.tex
TEXSOURCES += smr91.tex
TEXSOURCES += treatment_channel.tex
TEXSOURCES += treatment_plate.tex
TEXSOURCES += twopoint.tex

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
baseflow.dvi   radialflow.pdf: nozzle_subsonic.eps nozzle_supersonic.eps

# Images, Source Data, References, etc.
EXTRA_DIST  = astrommurray1011.eps helm.dot help.eps help.svg
EXTRA_DIST += cev_summary1.eps cev_summary2.eps cev_summary_fpg.eps
EXTRA_DIST += nozzle_schematic.eps nozzle_subsonic.eps nozzle_supersonic.eps
EXTRA_DIST += inhomogeneity1.gz inhomogeneity1.eps
EXTRA_DIST += inhomogeneity2.gz inhomogeneity2.eps
EXTRA_DIST += symplanenorm.eps
EXTRA_DIST += references.bib
EXTRA_DIST += writeup.sh

# Build and install writeups
all-local:       pdf-local
pdf-local:       $(PDFS)
PDFS           = $(TEXSOURCES:.tex=.pdf) reacting_main.pdf
PDFS          += viscous_nondim.pdf thicknesses.pdf
writeupsdir    = $(pkgdatadir)/writeups
writeups_DATA  = $(PDFS)

# Files to clean up based on common TeX-based tools
CLEANFILES  = *.aux *.bbl *.blg *.dvi *.log *.out *.pdf *.ps *.toc
CLEANFILES += *-generated.pdf

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
