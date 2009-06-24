## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(GSL_CFLAGS) $(LOG4CXX_CFLAGS) $(BOOST_CPPFLAGS)
AM_CPPFLAGS   += $(P3DFFT_CFLAGS) $(FFTW3_CFLAGS)
AM_LDFLAGS     =
bin_PROGRAMS   =
check_PROGRAMS =

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@
CC=@MPICC@

check_PROGRAMS           += test_gsl_bspline
test_gsl_bspline_SOURCES  = test_gsl_bspline.cpp
test_gsl_bspline_LDADD    = $(GSL_LIBS) $(LOG4CXX_LIBS)

check_PROGRAMS                 += test_gsl_integrate_qng
test_gsl_integrate_qng_SOURCES  = test_gsl_integrate_qng.cpp
test_gsl_integrate_qng_LDADD    = $(GSL_LIBS) $(LOG4CXX_LIBS)

check_PROGRAMS             += test_bspline_deriv
test_bspline_deriv_SOURCES  = test_bspline_deriv.cpp
test_bspline_deriv_LDADD    = $(GSL_LIBS) $(LOG4CXX_LIBS)

check_PROGRAMS           += test_pencil_grid
test_pencil_grid_SOURCES  = test_pencil_grid.cpp
test_pencil_grid_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS      += test_pencil
test_pencil_SOURCES  = test_pencil.cpp
test_pencil_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS      += test_complex
test_complex_SOURCES = test_complex.cpp
test_complex_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS         += test_underling
test_underling_SOURCES  = test_underling.cpp
test_underling_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_underling_LDADD   += ../suzerain/libunderling.la

check_PROGRAMS              += test_suzerain_error
test_suzerain_error_SOURCES  = test_suzerain_error.cpp
test_suzerain_error_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_suzerain_error_LDADD   += ../suzerain/libsuzerain.la

check_PROGRAMS                += test_bspline_operator
test_bspline_operator_SOURCES  = test_bspline_operator.cpp
test_bspline_operator_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_bspline_operator_LDADD   += ../suzerain/libsuzerain.la
test_bspline_operator_LDADD   += $(GSL_LIBS)
test_bspline_operator_LDADD   += $(LOG4CXX_LIBS)

TESTS_ENVIRONMENT = LOG4CXX_CONFIGURATION=$(srcdir)/log4cxx.properties
TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
