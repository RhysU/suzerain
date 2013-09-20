## Process this file with automake to produce Makefile.in

include ../Common.am

###
### A convenience library for logic used across this subdirectory
###

noinst_LTLIBRARIES    += libperfect.la
libperfect_la_SOURCES  = common_block.cpp
libperfect_la_SOURCES += common_block.hpp
libperfect_la_SOURCES += driver.cpp
libperfect_la_SOURCES += driver.hpp
libperfect_la_SOURCES += hybrid_operator.cpp
libperfect_la_SOURCES += hybrid_operator.hpp
libperfect_la_SOURCES += isothermal_mass_operator.cpp
libperfect_la_SOURCES += isothermal_mass_operator.hpp
libperfect_la_SOURCES += linearize_type.hpp
libperfect_la_SOURCES += manufactured_solution.cpp
libperfect_la_SOURCES += manufactured_solution.hpp
libperfect_la_SOURCES += navier_stokes.cpp
libperfect_la_SOURCES += navier_stokes.hpp
libperfect_la_SOURCES += nonlinear_operator.cpp
libperfect_la_SOURCES += nonlinear_operator.hpp
libperfect_la_SOURCES += nonreflecting_treatment.cpp
libperfect_la_SOURCES += nonreflecting_treatment.hpp
libperfect_la_SOURCES += nsctpl_rholut_fwd.hpp
libperfect_la_SOURCES += nsctpl_rholut.hpp
libperfect_la_SOURCES += perfect.cpp
libperfect_la_SOURCES += perfect.hpp
libperfect_la_SOURCES += quantities.cpp
libperfect_la_SOURCES += quantities.hpp
libperfect_la_SOURCES += scenario_definition.cpp
libperfect_la_SOURCES += scenario_definition.hpp

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
perfect_advance_SOURCES += main_advance_svnrev.c
BUILT_SOURCES           += main_advance_svnrev.h
perfect_advance_LDADD    = ./libperfect.la ../../suzerain/libsuzerain.la
main_advance_svnrev.h: $(perfect_advance_SOURCES) $(libperfect_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

##
## Old, non-suzerain::support::application_base channel application(s)
## TODO Eliminate these as part of closing ticket #2736
##

bin_PROGRAMS         += perfect_mean
perfect_mean_SOURCES  = perfect_mean.cpp
perfect_mean_SOURCES += perfect_mean_svnrev.c
BUILT_SOURCES        += perfect_mean_svnrev.h
perfect_mean_LDADD    = ./libperfect.la ../../suzerain/libsuzerain.la
perfect_mean_svnrev.h: $(perfect_mean_SOURCES) $(libperfect_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

##
## Tests for perfect gas applications
##
## Identical test logic varying only in the operator choices (e.g. --explicit
## versus implicit) are broken out into separate, templated tests to increase
## the parallelism available during 'make check'.
##

dist_check_SCRIPTS += test_setup.sh # Common initialization

# Restart handling using coefficient-based restart files
TESTS              += test_restart_explicit.sh
TESTS              += test_restart_rhome_xyz.sh
TESTS              += test_restart_rhome_y.sh
dist_check_SCRIPTS += test_restart_template.sh
dist_check_SCRIPTS += test_restart_explicit.sh
dist_check_SCRIPTS += test_restart_rhome_xyz.sh
dist_check_SCRIPTS += test_restart_rhome_y.sh

# Serial/parallel equivalence for degenerate parallel decompositions
# Not important from a utility perspective, but great for hitting edge cases
TESTS              += test_dparallel.sh
dist_check_SCRIPTS += test_dparallel.sh

# Serial/parallel equivalence for normal parallel decompositions
TESTS              += test_parallel_explicit.sh
TESTS              += test_parallel_rhome_xyz.sh
TESTS              += test_parallel_rhome_y.sh
dist_check_SCRIPTS += test_parallel_template.sh
dist_check_SCRIPTS += test_parallel_explicit.sh
dist_check_SCRIPTS += test_parallel_rhome_xyz.sh
dist_check_SCRIPTS += test_parallel_rhome_y.sh

# Read canned, legacy fields to test backwards compatibility
TESTS              += test_fields.sh
dist_check_SCRIPTS += test_fields.sh

# Adding fluctuations at startup
TESTS              += test_fluct.sh
dist_check_SCRIPTS += test_fluct.sh

# Restart handling using collocation-value-based restart files
TESTS              += test_physical_explicit.sh
TESTS              += test_physical_rhome_xyz.sh
TESTS              += test_physical_rhome_y.sh
dist_check_SCRIPTS += test_physical_template.sh
dist_check_SCRIPTS += test_physical_explicit.sh
dist_check_SCRIPTS += test_physical_rhome_xyz.sh
dist_check_SCRIPTS += test_physical_rhome_y.sh

# Ensure various implicit solver types give roughly equivalent results
TESTS              += test_implicit_solvers.sh
dist_check_SCRIPTS += test_implicit_solvers.sh

# Debugging issues around the tricks in largo_state will be awful
# in the larger context of a flow solver.  Beat on it specifically.
TESTS                     += test_largo_state
EXTRA_PROGRAMS            += test_largo_state
test_largo_state_SOURCES   = test_largo_state.cpp
test_largo_state_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_largo_state_LDADD     = ../../suzerain/libsuzerain.la
test_largo_state_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_largo_state_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

##
## Miscellaneous manual testing scripts
##

dist_check_SCRIPTS += manual.nrbc/pulse_perfect_init.sh
dist_check_SCRIPTS += manual.nrbc/pulse_perfect.py
dist_check_SCRIPTS += manual.nrbc/ticket2399.sh

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
