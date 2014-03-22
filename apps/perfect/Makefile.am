## Process this file with automake to produce Makefile.in

include ../Common.am

###
### A installed library bundling logic used across this subdirectory.
### Installation reduces the footprint of using multiple binaries,
### and, more importantly, provides better RPATH semantics.
###

lib_LTLIBRARIES                += libsuzerain_perfect.la
libsuzerain_perfect_la_LDFLAGS  = $(AM_LDFLAGS) -release $(PACKAGE_VERSION)
libsuzerain_perfect_la_LIBADD   = ../../suzerain/libsuzerain.la
libsuzerain_perfect_la_SOURCES  = common_block.cpp
libsuzerain_perfect_la_SOURCES += common_block.hpp
libsuzerain_perfect_la_SOURCES += definition_scenario.cpp
libsuzerain_perfect_la_SOURCES += definition_scenario.hpp
libsuzerain_perfect_la_SOURCES += driver.cpp
libsuzerain_perfect_la_SOURCES += driver.hpp
libsuzerain_perfect_la_SOURCES += linearize_type.hpp
libsuzerain_perfect_la_SOURCES += manufactured_solution.cpp
libsuzerain_perfect_la_SOURCES += manufactured_solution.hpp
libsuzerain_perfect_la_SOURCES += navier_stokes.cpp
libsuzerain_perfect_la_SOURCES += navier_stokes.hpp
libsuzerain_perfect_la_SOURCES += nsctpl_rholut_fwd.hpp
libsuzerain_perfect_la_SOURCES += nsctpl_rholut.hpp
libsuzerain_perfect_la_SOURCES += operator_hybrid_isothermal.cpp
libsuzerain_perfect_la_SOURCES += operator_hybrid_isothermal.hpp
libsuzerain_perfect_la_SOURCES += operator_mass_isothermal.cpp
libsuzerain_perfect_la_SOURCES += operator_mass_isothermal.hpp
libsuzerain_perfect_la_SOURCES += operator_nonlinear.cpp
libsuzerain_perfect_la_SOURCES += operator_nonlinear.hpp
libsuzerain_perfect_la_SOURCES += perfect.cpp
libsuzerain_perfect_la_SOURCES += perfect.hpp
libsuzerain_perfect_la_SOURCES += slowgrowth_type.hpp
libsuzerain_perfect_la_SOURCES += treatment_nonreflecting.cpp
libsuzerain_perfect_la_SOURCES += treatment_nonreflecting.hpp

##
## Applications
##

bin_PROGRAMS             += perfect_initial
perfect_initial_SOURCES   = main_initial.cpp
perfect_initial_LDADD     = libsuzerain_perfect.la
perfect_initial_CPPFLAGS  = $(AM_CPPFLAGS) $(REVISIONSTR)

bin_PROGRAMS             += perfect_advance
perfect_advance_SOURCES   = main_advance.cpp
perfect_advance_LDADD     = libsuzerain_perfect.la
perfect_advance_CPPFLAGS  = $(AM_CPPFLAGS) $(REVISIONSTR)

bin_PROGRAMS             += perfect_summary
perfect_summary_SOURCES   = main_summary.cpp
perfect_summary_LDADD     = libsuzerain_perfect.la
perfect_summary_CPPFLAGS  = $(AM_CPPFLAGS) $(REVISIONSTR)

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

##
## Miscellaneous manual testing scripts
## Adding a directory to EXTRA_DIST recursively distributes the contents
##
EXTRA_DIST += manual.nrbc
EXTRA_DIST += manual.gold

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
