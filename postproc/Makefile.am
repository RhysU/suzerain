## Process this file with automake to produce Makefile.in

postprocdir       = $(pkgdatadir)/postproc
postproc_DATA     = perfect.decl
postproc_SCRIPTS  = channel.m
postproc_SCRIPTS += datatimes.sh
postproc_SCRIPTS += frame.sh
postproc_SCRIPTS += gb2ge.m
postproc_SCRIPTS += gb.py
postproc_SCRIPTS += gplot
postproc_SCRIPTS += orderest
postproc_SCRIPTS += perfect_decl.py
postproc_SCRIPTS += perfect.py
postproc_SCRIPTS += plan_suzerain.py
postproc_SCRIPTS += plot_blparam.py
postproc_SCRIPTS += plot_stats.py
postproc_SCRIPTS += plot_status.py
postproc_SCRIPTS += propagation.py
postproc_SCRIPTS += reacting_summary.py
postproc_SCRIPTS += relaminarize.bl
postproc_SCRIPTS += remove_hdf5.py
postproc_SCRIPTS += rename_hdf5.py
postproc_SCRIPTS += select_files.py
postproc_SCRIPTS += sgr2suz.py
postproc_SCRIPTS += spectra.py
postproc_SCRIPTS += stationary.bl
postproc_SCRIPTS += summary_surf.py
postproc_SCRIPTS += turbulence.bl

# To automatically synchronize uncertainty and postprocessing computations
# generate Python assignments for pointwise computations using perfect.decl.
# The resulting file is COMMITTED as SymPy prerequisites are touchy.
$(abs_srcdir)/perfect_decl.raw: perfect.decl propagation.py
	@echo "# Generating $@ from $< (sensitive to SymPy version)" | tee "$@.tmp"
	md5sum $+ | sed -e 's/^/# Based upon MD5 /' >> "$@.tmp"
	"$(srcdir)/propagation.py" -d "$<" dec   >> "$@.tmp"
	mv "$@.tmp" "$@"

# Convert the perfect_decl.raw output into valid Python code
CLEANFILES = perfect_decl.py *.pyc
perfect_decl.py: $(abs_srcdir)/perfect_decl.raw Makefile.am
	@echo "# Generating $@ from $+" | tee "$@"
	@md5sum $+ | sed -e 's/^/# Based upon MD5 /'                         >> "$@"
	@echo "def pointwise(gam, Ma, Re, Pr, y, bar, tilde, local, tke):"   >> "$@"
	@echo '    """'                                                      >> "$@"
	@echo '    Populate Bunch-like bar, tilde, local, tke from pointwise'>> "$@"
	@echo '    collocation point data in bar on gridpoints y per scalars'>> "$@"
	@echo '    gamma, Mach number, Reynolds number, and Prandtl number.' >> "$@"
	@echo '    See, e.g., https://pypi.python.org/pypi/bunch re: Bunch.' >> "$@"
	@echo '    """'                                                      >> "$@"
	@$(SED) -e 's/\/\//# /'        \
	        -e 's/(y)//g'          \
	        -e 's/bar_/bar./g'     \
	        -e 's/tilde_/tilde./g' \
	        -e 's/local_/local./g' \
	        -e 's/tke_/tke./g'     \
	        -e 's/^./    &/' "$<"  >> "$@"

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
