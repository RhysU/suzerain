## Process this file with automake to produce Makefile.in

AM_CPPFLAGS        = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS       += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS       += $(BOOST_CPPFLAGS)
AM_LDFLAGS         = # Append below
bin_PROGRAMS       = # Append below
check_PROGRAMS     = # Append below
dist_check_SCRIPTS = # Append below
noinst_HEADERS     = # Append below
BUILT_SOURCES      = # Append below
CLEANFILES         = # Append below
TESTS              = # Append below
TESTS_ENVIRONMENT  = # Append below

noinst_HEADERS += test_tools.hpp
noinst_HEADERS += fixtures.hpp

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CC=@MPICC@
CXX=@MPICXX@

# C-based to avoid some indecipherable Ubuntu 11.10 linking issues
TESTS              += test_error
check_PROGRAMS     += test_error
test_error_SOURCES  = test_error.c
test_error_CFLAGS   = $(AM_CFLAGS) $(FFTW3_MPI_THREADS_CFLAGS)
test_error_LDADD    = ../underling/libunderling.la

TESTS                += test_fftw_1D
check_PROGRAMS       += test_fftw_1D
test_fftw_1D_SOURCES  = test_fftw_1D.cpp
test_fftw_1D_CXXFLAGS = $(AM_CXXFLAGS) $(FFTW3_MPI_THREADS_CFLAGS)
test_fftw_1D_LDADD    = ../underling/libunderling.la
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_underling.sh
dist_check_SCRIPTS     += test_underling.sh
check_PROGRAMS         += test_underling
test_underling_SOURCES  = test_underling.cpp
test_underling_CXXFLAGS = $(AM_CXXFLAGS) $(FFTW3_MPI_THREADS_CFLAGS)
test_underling_LDADD    = ../underling/libunderling.la
test_underling_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_underling_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                      += test_underling_fftw.sh
dist_check_SCRIPTS         += test_underling_fftw.sh
check_PROGRAMS             += test_underling_fftw
test_underling_fftw_SOURCES  = test_underling_fftw.cpp
test_underling_fftw_CXXFLAGS = $(AM_CXXFLAGS) $(FFTW3_MPI_THREADS_CFLAGS)
test_underling_fftw_LDADD    = ../underling/libunderling.la
test_underling_fftw_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_underling_fftw_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

#############################################################################
### MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS ###
#############################################################################

# Build test programs without running 'make check'
.PHONY: check-programs
check-programs: $(check_PROGRAMS)

# Extra things to clean up
clean-local:
	rm -f *.gcno *.gcda

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS = -i underling -i underling_c_binding
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
