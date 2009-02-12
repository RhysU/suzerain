## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(GSL_CFLAGS) $(LOG4CXX_CFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS     = $(GSL_LIBS)   $(LOG4CXX_LIBS)   $(BOOST_LDFLAGS)
bin_PROGRAMS   =
check_PROGRAMS = 

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CC=@MPICC@

# Progams
bin_PROGRAMS       += helloworld
helloworld_SOURCES  = helloworld.cc

bin_PROGRAMS        += driver_sine
driver_sine_SOURCES  = driver_sine.c
driver_sine_LDADD    = $(FFTW3_LIBS) $(P3DFFT_LIBS)
driver_sine_CPPFLAGS = $(FFTW3_CFLAGS) $(P3DFFT_CFLAGS)

bin_PROGRAMS             += test_gsl_bspline
test_gsl_bspline_SOURCES  = test_gsl_bspline.cc

# Unit test executables
check_PROGRAMS      += test_example
test_example_SOURCES = test_example.cc
test_example_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
