## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(MKL_CFLAGS) $(GSL_CFLAGS) $(LOG4CXX_CFLAGS)
AM_CPPFLAGS   += $(BOOST_CPPFLAGS) $(P3DFFT_CFLAGS) $(FFTW3_CFLAGS)
AM_LDFLAGS     =
bin_PROGRAMS   =
check_PROGRAMS =
noinst_HEADERS =

noinst_HEADERS += test_tools.hpp

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@
CC=@MPICC@

check_PROGRAMS          += test_blas_et_al
test_blas_et_al_SOURCES  = test_blas_et_al.c
test_blas_et_al_LDADD    = $(GSL_LIBS)
test_blas_et_al_LDADD   += ../suzerain/libsuzerain.la

check_PROGRAMS           += test_pencil_grid
test_pencil_grid_SOURCES  = test_pencil_grid.cpp
test_pencil_grid_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencil_grid_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS      += test_pencil
test_pencil_SOURCES  = test_pencil.cpp
test_pencil_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencil_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS      += test_complex
test_complex_SOURCES = test_complex.cpp
test_complex_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_complex_LDADD  += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS         += test_underling
test_underling_SOURCES  = test_underling.cpp
test_underling_LDADD    = ../suzerain/libunderling.la
test_underling_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_underling_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS     += test_error
test_error_SOURCES  = test_error.cpp
test_error_LDADD    = ../suzerain/libsuzerain.la
test_error_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_error_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS       += test_bspline
test_bspline_SOURCES  = test_bspline.cpp
test_bspline_LDADD    = ../suzerain/libsuzerain.la
test_bspline_LDADD   += $(GSL_LIBS)
test_bspline_LDADD   += $(LOG4CXX_LIBS)
test_bspline_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_bspline_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS            += test_bspline_eval
test_bspline_eval_SOURCES  = test_bspline_eval.cpp
test_bspline_eval_LDADD    = ../suzerain/libsuzerain.la
test_bspline_eval_LDADD   += $(LOG4CXX_LIBS)
test_bspline_eval_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_bspline_eval_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS_ENVIRONMENT = LOG4CXX_CONFIGURATION=$(srcdir)/log4cxx.properties
TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
