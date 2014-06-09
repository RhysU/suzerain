## Process this file with automake to produce Makefile.in

postprocdir       = $(pkgdatadir)/postproc
postproc_SCRIPTS  = channel.m
postproc_SCRIPTS += datatimes.sh
postproc_SCRIPTS += frame.sh
postproc_SCRIPTS += gb2ge.m
postproc_SCRIPTS += gb.py
postproc_SCRIPTS += gplot
postproc_SCRIPTS += orderest
postproc_SCRIPTS += perfect.py
postproc_SCRIPTS += plan_suzerain.py
postproc_SCRIPTS += plot_blparam.py
postproc_SCRIPTS += plot_stats.py
postproc_SCRIPTS += plot_status.py
postproc_SCRIPTS += reacting_summary.py
postproc_SCRIPTS += turbulence.bl
postproc_SCRIPTS += remove_hdf5.py
postproc_SCRIPTS += rename_hdf5.py
postproc_SCRIPTS += select_files.py
postproc_SCRIPTS += sgr2suz.py
postproc_SCRIPTS += spectra.py
postproc_SCRIPTS += stationary.bl
postproc_SCRIPTS += summary_surf.py

# To automatically synchronize uncertainty and postprocessing computations
# generate Python code for pointwise computations using preproc/perfect.decl
perfect.py: $(top_srcdir)/preproc/perfect.decl Makefile.am
	@echo Generating $@ from $<
	@echo "# Generated from $+ at $(shell date)"                          > "$@"
	@echo "def pointwise(gam, Ma, Re, Pr, y, bar, tilde, local, tke):"   >> "$@"
	@echo '    """'                                                      >> "$@"
	@echo '    Populate Bunch-like bar, tilde, local, tke from pointwise'>> "$@"
	@echo '    collocation point data in bar on gridpoints y per scalars'>> "$@"
	@echo '    gamma, Mach number, Reynolds number, and Prandtl number.' >> "$@"
	@echo '    See, e.g., https://pypi.python.org/pypi/bunch re: Bunch.' >> "$@"
	@echo '    """'                                                      >> "$@"
	@$(SED) -e 's/\/\//# /g'       \
	        -e 's/(y)/   /g'       \
	        -e 's/bar_/bar./g'     \
	        -e 's/tilde_/tilde./g' \
	        -e 's/local_/local./g' \
	        -e 's/tke_/tke./g'     \
	        -e 's/^/    /' "$<" >> "$@"

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
