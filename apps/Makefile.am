## Process this file with automake to produce Makefile.in

include Common.am

# Build targets in current directory before any subdirectories
SUBDIRS = . perfect reacting

##
## Miscellaneous programs and/or performance drivers
##

noinst_PROGRAMS      += conservation
conservation_SOURCES  = conservation.cpp
conservation_LDADD    = ../suzerain/libsuzerain.la

noinst_PROGRAMS    += operators
operators_SOURCES   = operators.cpp
operators_LDADD     = ../suzerain/libsuzerain.la

bin_PROGRAMS               += pencil_grid_bench
pencil_grid_bench_LDADD     = ../suzerain/libsuzerain.la
pencil_grid_bench_SOURCES   = pencil_grid_bench.cpp
pencil_grid_bench_CPPFLAGS  = $(AM_CPPFLAGS) $(REVISIONFLAG)

##
## Common test infrastructure and tests
##

dist_check_SCRIPTS += test_infrastructure.sh

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
