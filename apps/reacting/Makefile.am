## Process this file with automake to produce Makefile.in

include ../Common.am

##
## While there are some files here that will compile w/out Antioch,
## the overall capability does not make sense without Antioch support,
## so skip unless Antioch is detected at configure time.
##
if BUILD_REACTING

AM_CPPFLAGS += $(ANTIOCH_CFLAGS)

###
### A installed library bundling logic used across this subdirectory.
### Installation reduces the footprint of using multiple binaries,
### and, more importantly, provides better RPATH semantics.
###

lib_LTLIBRARIES                 += libsuzerain_reacting.la
libsuzerain_reacting_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libsuzerain_reacting_la_LIBADD   = ../../suzerain/libsuzerain.la
libsuzerain_reacting_la_LIBADD  += $(ANTIOCH_LIBS)
libsuzerain_reacting_la_SOURCES  = antioch_constitutive.cpp
libsuzerain_reacting_la_SOURCES += antioch_constitutive.hpp
libsuzerain_reacting_la_SOURCES += definition_channel.cpp
libsuzerain_reacting_la_SOURCES += definition_channel.hpp
libsuzerain_reacting_la_SOURCES += definition_filter.cpp
libsuzerain_reacting_la_SOURCES += definition_filter.hpp
libsuzerain_reacting_la_SOURCES += driver.cpp
libsuzerain_reacting_la_SOURCES += driver.hpp
libsuzerain_reacting_la_SOURCES += manufactured_solution.cpp
libsuzerain_reacting_la_SOURCES += manufactured_solution.hpp
libsuzerain_reacting_la_SOURCES += nsctpl_fwd.hpp
libsuzerain_reacting_la_SOURCES += nsctpl.hpp
libsuzerain_reacting_la_SOURCES += operator_explicit.cpp
libsuzerain_reacting_la_SOURCES += operator_explicit.hpp
libsuzerain_reacting_la_SOURCES += operator_hybrid.cpp
libsuzerain_reacting_la_SOURCES += operator_hybrid.hpp
libsuzerain_reacting_la_SOURCES += operator_nonlinear_fwd.hpp
libsuzerain_reacting_la_SOURCES += operator_nonlinear.hpp
libsuzerain_reacting_la_SOURCES += quantities.cpp
libsuzerain_reacting_la_SOURCES += quantities.hpp
libsuzerain_reacting_la_SOURCES += reacting.cpp
libsuzerain_reacting_la_SOURCES += reacting.hpp
libsuzerain_reacting_la_SOURCES += reacting_ndx.hpp
libsuzerain_reacting_la_SOURCES += treatment_nonreflecting.cpp
libsuzerain_reacting_la_SOURCES += treatment_nonreflecting.hpp

##
## Applications
##

bin_PROGRAMS              += reacting_init
reacting_init_SOURCES      = main_init.cpp
reacting_init_LDADD        = libsuzerain_reacting.la

bin_PROGRAMS              += reacting_advance
reacting_advance_SOURCES   = main_advance.cpp
reacting_advance_LDADD     = libsuzerain_reacting.la

##
## Tests for reacting flow applications
##
## Identical test logic varying only in the operator choices (e.g. --explicit
## versus implicit) are broken out into separate, templated tests to increase
## the parallelism available during 'make check'.
##

dist_check_SCRIPTS += test_setup.sh              # Common initialization...
dist_check_SCRIPTS += test_setup_multispecies.sh # ...including for multispecies

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
#TESTS              += test_fields.sh
#dist_check_SCRIPTS += test_fields.sh

# Adding fluctuations at startup
TESTS              += test_fluct.sh
dist_check_SCRIPTS += test_fluct.sh

# Restart handling using collocation-value-based restart files
TESTS              += test_physical_explicit.sh
TESTS              += test_physical_implicit.sh
dist_check_SCRIPTS += test_physical_template.sh
dist_check_SCRIPTS += test_physical_explicit.sh
dist_check_SCRIPTS += test_physical_implicit.sh

# Unit tests for antioch_constitutive class
TESTS                             += test_antioch_constitutive.sh
dist_check_SCRIPTS                += test_antioch_constitutive.sh

check_PROGRAMS                     = test_antioch_constitutive
test_antioch_constitutive_SOURCES  = test_antioch_constitutive.cpp
test_antioch_constitutive_LDADD    = libsuzerain_reacting.la

# Augment test environment for 'make check'.
# One may override by manipulating $TESTS_ENVIRONMENT using the command line.
AM_TESTS_ENVIRONMENT += export ANTIOCH_DATA_DIR="$(abs_top_srcdir)/apps/reacting/data";

# Multispecies restart
TESTS              += test_restart_multispecies_explicit.sh
TESTS              += test_restart_multispecies_implicit.sh
dist_check_SCRIPTS += test_restart_multispecies_template.sh
dist_check_SCRIPTS += test_restart_multispecies_explicit.sh
dist_check_SCRIPTS += test_restart_multispecies_implicit.sh

endif # if BUILD_REACTING

##
## Adding a directory to EXTRA_DIST recursively distributes the contents
##
EXTRA_DIST += data
EXTRA_DIST += scripts

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
