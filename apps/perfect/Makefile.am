## Process this file with automake to produce Makefile.in

include ../Common.am

##
## Channel-based applications
##

bin_PROGRAMS         += channel_init
channel_init_LDADD    = ../libchannel.la
channel_init_SOURCES  = channel_init.cpp
channel_init_SOURCES += channel_init_svnrev.c
BUILT_SOURCES        += channel_init_svnrev.h
channel_init_svnrev.h: $(channel_init_SOURCES)                       \
                       $(__top_builddir__apps_libchannel_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

bin_PROGRAMS    += channel
channel_LDADD    = ../libchannel.la
channel_SOURCES  = channel_main.cpp
channel_SOURCES += channel_treatment.hpp
channel_SOURCES += nonlinear_fwd.hpp nonlinear.hpp
channel_SOURCES += explicit_op.cpp explicit_op.hpp
channel_SOURCES += hybrid_op.cpp   hybrid_op.hpp
channel_SOURCES += ../timers.hpp
channel_SOURCES += channel_main_svnrev.c
BUILT_SOURCES   += channel_main_svnrev.h
channel_main_svnrev.h: $(channel_SOURCES)                            \
                       $(__top_builddir__apps_libchannel_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

bin_PROGRAMS         += channel_mean
channel_mean_LDADD    = ../libchannel.la
channel_mean_SOURCES  = channel_mean.cpp
channel_mean_SOURCES += channel_mean_svnrev.c
BUILT_SOURCES        += channel_mean_svnrev.h
channel_mean_svnrev.h: $(channel_mean_SOURCES)                       \
                       $(__top_builddir__apps_libchannel_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

##
## Tests for channel applications
##
## Identical test logic varying only in the operator choices (e.g. --explicit
## versus implicit) are broken out into separate, templated tests to increase
## the parallelism available during 'make check'.
##

dist_check_SCRIPTS += test_channel_setup.sh     # Common initialization

# Restart handling using coefficient-based restart files
TESTS              += test_channel_restart_explicit.sh
TESTS              += test_channel_restart_implicit.sh
dist_check_SCRIPTS += test_channel_restart_template.sh
dist_check_SCRIPTS += test_channel_restart_explicit.sh
dist_check_SCRIPTS += test_channel_restart_implicit.sh

# Serial/parallel equivalence for degenerate parallel decompositions
# Not important from a utility perspective, but great for hitting edge cases
TESTS              += test_channel_dparallel.sh
dist_check_SCRIPTS += test_channel_dparallel.sh

# Serial/parallel equivalence for normal parallel decompositions
TESTS              += test_channel_parallel_explicit.sh
TESTS              += test_channel_parallel_implicit.sh
dist_check_SCRIPTS += test_channel_parallel_template.sh
dist_check_SCRIPTS += test_channel_parallel_explicit.sh
dist_check_SCRIPTS += test_channel_parallel_implicit.sh

# Read canned, legacy fields to test backwards compatibility
TESTS              += test_channel_fields.sh
dist_check_SCRIPTS += test_channel_fields.sh

# Adding fluctuations at startup
TESTS              += test_channel_fluct.sh
dist_check_SCRIPTS += test_channel_fluct.sh

# Restart handling using collocation-value-based restart files
TESTS              += test_channel_physical_explicit.sh
TESTS              += test_channel_physical_implicit.sh
dist_check_SCRIPTS += test_channel_physical_template.sh
dist_check_SCRIPTS += test_channel_physical_explicit.sh
dist_check_SCRIPTS += test_channel_physical_implicit.sh

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
