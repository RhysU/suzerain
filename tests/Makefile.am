## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(MKL_CFLAGS) $(GSL_CFLAGS) $(LOG4CXX_CFLAGS)
AM_CPPFLAGS   += $(BOOST_CPPFLAGS) $(P3DFFT_CFLAGS) $(FFTW3_CFLAGS)
AM_CPPFLAGS   += $(EIGEN_CFLAGS)
AM_LDFLAGS     = # Append below
bin_PROGRAMS   = # Append below
check_PROGRAMS = # Append below
noinst_HEADERS = # Append below
EXTRA_DIST     = # Append below
CLEANFILES     = # Append below

noinst_HEADERS += test_tools.hpp

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@
CC=@MPICC@

# Clean up from Intel Compiler's -prof-gen output
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
# Clean up from Intel Compiler's precompiled header output
CLEANFILES += *.pchi

check_PROGRAMS     += test_error
test_error_SOURCES  = test_error.cpp
test_error_LDADD    = ../suzerain/libsuzerain.la
test_error_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_error_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS          += test_blas_et_al
test_blas_et_al_SOURCES  = test_blas_et_al.c
test_blas_et_al_LDADD    = $(GSL_LIBS)
test_blas_et_al_LDADD   += ../suzerain/libsuzerain.la

check_PROGRAMS       += test_bspline
test_bspline_SOURCES  = test_bspline.cpp
test_bspline_LDADD    = ../suzerain/libsuzerain.la
test_bspline_LDADD   += $(GSL_LIBS)
test_bspline_LDADD   += $(LOG4CXX_LIBS)
test_bspline_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_bspline_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS      += test_complex
test_complex_SOURCES = test_complex.cpp
test_complex_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_complex_LDADD  += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS           += test_pencil_grid
test_pencil_grid_SOURCES  = test_pencil_grid.cpp
test_pencil_grid_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencil_grid_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS      += test_pencil
test_pencil_SOURCES  = test_pencil.cpp
test_pencil_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencil_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

check_PROGRAMS           += test_orthonormal
test_orthonormal_SOURCES  = test_orthonormal.cpp
test_orthonormal_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_orthonormal_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
EXTRA_DIST               += test_orthonormal.sage

check_PROGRAMS          += test_richardson
test_richardson_SOURCES  = test_richardson.c
test_richardson_LDADD    = $(GSL_LIBS)
test_richardson_LDADD   += ../suzerain/libsuzerain.la

check_PROGRAMS           += test_timestepper
test_timestepper_SOURCES  = test_timestepper.c
test_timestepper_LDADD    = $(GSL_LIBS)
test_timestepper_LDADD   += ../suzerain/libsuzerain.la

check_PROGRAMS       += test_fftw_1D
test_fftw_1D_SOURCES  = test_fftw_1D.cpp
test_fftw_1D_LDADD    = $(FFTW3_LIBS)
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)

check_PROGRAMS                   += test_fftw_multi_array_1D
test_fftw_multi_array_1D_SOURCES  = test_fftw_multi_array_1D.cpp
test_fftw_multi_array_1D_LDADD    = $(FFTW3_LIBS)
test_fftw_multi_array_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)

check_PROGRAMS                   += test_fftw_multi_array_2D
test_fftw_multi_array_2D_SOURCES  = test_fftw_multi_array_2D.cpp
test_fftw_multi_array_2D_LDADD    = $(FFTW3_LIBS)
test_fftw_multi_array_2D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)

check_PROGRAMS         += test_underling
test_underling_SOURCES  = test_underling.cpp
test_underling_LDADD    = ../suzerain/libunderling.la
test_underling_LDADD   += $(GSL_LIBS)
test_underling_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_underling_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS_ENVIRONMENT = LOG4CXX_CONFIGURATION=$(srcdir)/log4cxx.properties
TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
