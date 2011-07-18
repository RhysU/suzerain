## Process this file with automake to produce Makefile.in

# Append new system job scripts onto this list
JOBSCRIPTS = qsub.lonestar

# Installation, distribution, and clean up details
jobscriptsdir      = $(pkgdatadir)/jobscripts
jobscripts_SCRIPTS = $(JOBSCRIPTS)
EXTRA_DIST         = $(JOBSCRIPTS:=.in)
CLEANFILES         = $(JOBSCRIPTS)

# Job scripts are generated following the sed example at
# http://www.gnu.org/software/automake/manual/automake.html#Scripts
do_subst = $(SED) -e 's,[@]bindir[@],$(bindir),g'         \
                  -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

$(JOBSCRIPTS) : % : $(srcdir)/%.in Makefile
	$(do_subst) < "$(srcdir)/$@.in" > "$@"
	chmod +x "$@"

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
