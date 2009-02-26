## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(GSL_CFLAGS) $(LOG4CXX_CFLAGS) $(BOOST_CPPFLAGS)
AM_CPPFLAGS   += $(P3DFFT_CFLAGS) $(FFTW3_CFLAGS)
AM_LDFLAGS     = $(BOOST_LDFLAGS)
bin_PROGRAMS   =
check_PROGRAMS =

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@

# Progams
bin_PROGRAMS       += helloworld
helloworld_SOURCES  = helloworld.cc
helloworld_LDADD    = $(GSL_LIBS)

bin_PROGRAMS        += driver_sine
driver_sine_SOURCES  = driver_sine.cc p3dfft_d.h
driver_sine_LDADD    = $(GSL_LIBS) $(LOG4CXX_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB)
driver_sine_LDADD   += $(P3DFFT_LIBS) $(FFTW3_LIBS)

bin_PROGRAMS             += test_gsl_bspline
test_gsl_bspline_SOURCES  = test_gsl_bspline.cc
test_gsl_bspline_LDADD    = $(GSL_LIBS) $(LOG4CXX_LIBS)

# Unit test executables
check_PROGRAMS      += test_example
test_example_SOURCES = test_example.cc
test_example_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

TESTS_ENVIRONMENT = LOG4CXX_CONFIGURATION=$(srcdir)/log4cxx.properties
TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
