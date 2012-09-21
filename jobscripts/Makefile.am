## Process this file with automake to produce Makefile.in
dist_check_SCRIPTS  = # Append below
TESTS               = # Append below

JOBSCRIPTS = qsub.lonestar qsub.ranger mms_case.perfect
UTILITIES  = qlast

# Installation, distribution, and clean up details
jobscriptsdir      = $(pkgdatadir)/jobscripts
jobscripts_SCRIPTS = $(JOBSCRIPTS) $(UTILITIES)
EXTRA_DIST         = qsub.tacc.in
CLEANFILES         = $(JOBSCRIPTS)

# Job scripts are generated following the sed-based example at
# http://www.gnu.org/software/automake/manual/automake.html#Scripts
do_subst = $(SED) -e 's,[@]bindir[@],$(bindir),g'         \
                  -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

$(JOBSCRIPTS) : Makefile

mms_case.perfect : mms_case.perfect.in
	$(do_subst) "$<" > "$@"
	chmod +x "$@"

qsub.lonestar : qsub.tacc.in
	$(do_subst) "$<" | $(SED) -e 's,[@]corespernode[@],12,g' > "$@"
	chmod +x "$@"

qsub.ranger : qsub.tacc.in
	$(do_subst) "$<" | $(SED) -e 's,[@]corespernode[@],16,g' > "$@"
	chmod +x "$@"

# Ensure bash-based job scripts have valid syntax
TESTS              += test_bash_syntax.sh
dist_check_SCRIPTS += test_bash_syntax.sh

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
