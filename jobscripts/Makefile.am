## Process this file with automake to produce Makefile.in

JOBSCRIPTS = qsub.lonestar qsub.ranger

# Installation, distribution, and clean up details
jobscriptsdir      = $(pkgdatadir)/jobscripts
jobscripts_SCRIPTS = $(JOBSCRIPTS)
EXTRA_DIST         = qsub.tacc.in
CLEANFILES         = $(JOBSCRIPTS)

# Job scripts are generated following the sed-based example at
# http://www.gnu.org/software/automake/manual/automake.html#Scripts
do_subst = $(SED) -e 's,[@]bindir[@],$(bindir),g'         \
                  -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

$(JOBSCRIPTS) : Makefile

qsub.lonestar : $(srcdir)/qsub.tacc.in
	$(do_subst) "$(srcdir)/$<" | $(SED) -e 's,[@]corespernode[@],12,g' > "$@"
	chmod +x "$@"

qsub.ranger : $(srcdir)/qsub.tacc.in
	$(do_subst) "$(srcdir)/$<" | $(SED) -e 's,[@]corespernode[@],16,g' > "$@"
	chmod +x "$@"

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
