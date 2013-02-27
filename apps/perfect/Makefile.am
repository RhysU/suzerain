## Process this file with automake to produce Makefile.in

include ../Common.am

###
### Logic common to many different binaries within this subdirectory
###

noinst_LTLIBRARIES    += libperfect.la
libperfect_la_SOURCES  = perfect.hpp               perfect.cpp
libperfect_la_SOURCES += scenario_definition.hpp   scenario_definition.cpp
libperfect_la_SOURCES += mean_quantities.hpp       mean_quantities.cpp
libperfect_la_SOURCES += manufactured_solution.hpp manufactured_solution.cpp
libperfect_la_SOURCES += nsctpl_rholut_fwd.hpp     nsctpl_rholut.hpp
libperfect_la_SOURCES += driver.hpp                driver.cpp

##
## Applications
##

bin_PROGRAMS         += perfect_init
perfect_init_SOURCES  = main_init.cpp
perfect_init_SOURCES += main_init_svnrev.c
BUILT_SOURCES        += main_init_svnrev.h
perfect_init_LDADD    = ./libperfect.la ../../suzerain/libsuzerain.la
main_init_svnrev.h: $(perfect_init_SOURCES) $(libperfect_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

bin_PROGRAMS            += perfect_advance
perfect_advance_SOURCES  = main_advance.cpp
perfect_advance_SOURCES += nonlinear_operator_fwd.hpp nonlinear_operator.hpp
perfect_advance_SOURCES += explicit_operator.cpp explicit_operator.hpp
perfect_advance_SOURCES += hybrid_operator.cpp   hybrid_operator.hpp
perfect_advance_SOURCES += main_advance_svnrev.c
BUILT_SOURCES           += main_advance_svnrev.h
perfect_advance_LDADD    = ./libperfect.la ../../suzerain/libsuzerain.la
main_advance_svnrev.h: $(perfect_advance_SOURCES) $(libperfect_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

##
## Old, non-suzerain::support::application_base channel application(s)
## TODO Eliminate these as part of closing ticket #2736
##

bin_PROGRAMS         += channel_mean
channel_mean_SOURCES  = channel_mean.cpp
channel_mean_SOURCES += channel_mean_svnrev.c
BUILT_SOURCES        += channel_mean_svnrev.h
channel_mean_LDADD    = ./libperfect.la ../../suzerain/libsuzerain.la
channel_mean_svnrev.h: $(channel_mean_SOURCES) $(libperfect_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

##
## Tests for perfect gas applications
##
## Identical test logic varying only in the operator choices (e.g. --explicit
## versus implicit) are broken out into separate, templated tests to increase
## the parallelism available during 'make check'.
##

dist_check_SCRIPTS += test_setup.sh     # Common initialization

# Restart handling using coefficient-based restart files
TESTS              += test_restart_explicit.sh
TESTS              += test_restart_implicit.sh
dist_check_SCRIPTS += test_restart_template.sh
dist_check_SCRIPTS += test_restart_explicit.sh
dist_check_SCRIPTS += test_restart_implicit.sh

# Serial/parallel equivalence for degenerate parallel decompositions
# Not important from a utility perspective, but great for hitting edge cases
TESTS              += test_dparallel.sh
dist_check_SCRIPTS += test_dparallel.sh

# Serial/parallel equivalence for normal parallel decompositions
TESTS              += test_parallel_explicit.sh
TESTS              += test_parallel_implicit.sh
dist_check_SCRIPTS += test_parallel_template.sh
dist_check_SCRIPTS += test_parallel_explicit.sh
dist_check_SCRIPTS += test_parallel_implicit.sh

# Read canned, legacy fields to test backwards compatibility
TESTS              += test_fields.sh
dist_check_SCRIPTS += test_fields.sh

# Adding fluctuations at startup
TESTS              += test_fluct.sh
dist_check_SCRIPTS += test_fluct.sh

# Restart handling using collocation-value-based restart files
TESTS              += test_physical_explicit.sh
TESTS              += test_physical_implicit.sh
dist_check_SCRIPTS += test_physical_template.sh
dist_check_SCRIPTS += test_physical_explicit.sh
dist_check_SCRIPTS += test_physical_implicit.sh

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
