## Process this file with automake to produce Makefile.in
dist_check_SCRIPTS  = # Append below
TESTS               = # Append below

JOBSCRIPTS  = qsub.lonestar qsub.ranger
UTILITIES   = qlast
UTILITIES  += channel.perfect  mms_case.perfect
UTILITIES  += channel.reacting mms_case.reacting

# Installation, distribution, and clean up details
jobscriptsdir      = $(pkgdatadir)/jobscripts
jobscripts_SCRIPTS = $(JOBSCRIPTS) $(UTILITIES)
EXTRA_DIST         = qsub.tacc.in
EXTRA_DIST        += channel.perfect.in  mms_case.perfect.in
EXTRA_DIST        += channel.reacting.in mms_case.reacting.in
CLEANFILES         = $(JOBSCRIPTS)

# Job scripts are generated following the sed-based example at
# http://www.gnu.org/software/automake/manual/automake.html#Scripts
do_subst = $(SED) -e 's,[@]bindir[@],$(bindir),g'         \
                  -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

%.perfect : %.perfect.in Makefile
	$(do_subst) "$<" > "$@"
	chmod +x "$@"

%.reacting : %.reacting.in Makefile
	$(do_subst) "$<" > "$@"
	chmod +x "$@"

qsub.lonestar : qsub.tacc.in Makefile
	$(do_subst) "$<" | $(SED) -e 's,[@]corespernode[@],12,g' > "$@"
	chmod +x "$@"

qsub.ranger : qsub.tacc.in Makefile
	$(do_subst) "$<" | $(SED) -e 's,[@]corespernode[@],16,g' > "$@"
	chmod +x "$@"

# Ensure bash-based job scripts have valid syntax
TESTS              += test_bash_syntax.sh
dist_check_SCRIPTS += test_bash_syntax.sh

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
