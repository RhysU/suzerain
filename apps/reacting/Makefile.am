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
### Logic common to many different binaries within this subdirectory
###

noinst_LTLIBRARIES     += libreacting.la
libreacting_la_SOURCES  = reacting_ndx.hpp
libreacting_la_SOURCES += reacting.hpp              reacting.cpp
libreacting_la_SOURCES += channel_definition.hpp    channel_definition.cpp
libreacting_la_SOURCES += filter_definition.hpp     filter_definition.cpp
libreacting_la_SOURCES += quantities.hpp            quantities.cpp
libreacting_la_SOURCES += manufactured_solution.hpp manufactured_solution.cpp
libreacting_la_SOURCES += nsctpl_fwd.hpp            nsctpl.hpp
libreacting_la_SOURCES += driver.hpp                driver.cpp
libreacting_la_SOURCES += antioch_constitutive.hpp
libreacting_la_SOURCES += antioch_constitutive.cpp
libreacting_la_LIBADD   = $(ANTIOCH_LIBS)

##
## Applications
##

bin_PROGRAMS          += reacting_init
reacting_init_SOURCES  = main_init.cpp
reacting_init_SOURCES += main_init_svnrev.c
BUILT_SOURCES         += main_init_svnrev.h
reacting_init_LDADD    = ./libreacting.la ../../suzerain/libsuzerain.la
main_init_svnrev.h: $(reacting_init_SOURCES) $(libreacting_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

bin_PROGRAMS             += reacting_advance
reacting_advance_SOURCES  = main_advance.cpp
reacting_advance_SOURCES += nonlinear_operator_fwd.hpp nonlinear_operator.hpp
reacting_advance_SOURCES += explicit_operator.cpp explicit_operator.hpp
reacting_advance_SOURCES += hybrid_operator.cpp hybrid_operator.hpp
reacting_advance_SOURCES += main_advance_svnrev.c
BUILT_SOURCES            += main_advance_svnrev.h
reacting_advance_LDADD    = ./libreacting.la ../../suzerain/libsuzerain.la
main_advance_svnrev.h: $(reacting_advance_SOURCES) $(libreacting_la_SOURCES)
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

##
## Tests for reacting flow applications
##
## Identical test logic varying only in the operator choices (e.g. --explicit
## versus implicit) are broken out into separate, templated tests to increase
## the parallelism available during 'make check'.
##

dist_check_SCRIPTS += test_setup.sh # Common initialization...
dist_check_SCRIPTS += test_setup_multispecies.sh # Common initialization (for multispecies)
AM_TESTS_FD_REDIRECT = 9>&2         # ...which uses fd 9 to get human attention

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
test_antioch_constitutive_LDADD    = ./libreacting.la ../../suzerain/libsuzerain.la

TESTS_ENVIRONMENT += ANTIOCH_DATA_DIR=$(abs_top_srcdir)/apps/reacting/data

# Multispecies restart
TESTS              += test_restart_multispecies_explicit.sh
dist_check_SCRIPTS += test_restart_multispecies_explicit.sh

endif # if BUILD_REACTING

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
