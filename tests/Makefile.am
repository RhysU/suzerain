## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(MKL_CFLAGS) $(GSL_CFLAGS) $(LOG4CXX_CFLAGS)
AM_CPPFLAGS   += $(BOOST_CPPFLAGS) $(P3DFFT_CFLAGS)
AM_CPPFLAGS   += $(EIGEN_CFLAGS)
AM_LDFLAGS     = # Append below
bin_PROGRAMS   = # Append below
check_PROGRAMS = # Append below
noinst_HEADERS = # Append below
CLEANFILES     = # Append below
EXTRA_DIST     = # Append below
TESTS          = # Append below

noinst_HEADERS += test_tools.hpp
noinst_HEADERS += test_underling_tools.hpp

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@
CC=@MPICC@

# Clean up from Intel Compiler's -prof-gen output
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
# Clean up from Intel Compiler's precompiled header output
CLEANFILES += *.pchi

# Provide a test-ready logging configuration
TESTS_ENVIRONMENT = LOG4CXX_CONFIGURATION=$(srcdir)/log4cxx.properties

TESTS              += test_error
check_PROGRAMS     += test_error
test_error_SOURCES  = test_error.cpp
test_error_LDADD    = ../suzerain/libsuzerain.la
test_error_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_error_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS           += test_os
check_PROGRAMS  += test_os
test_os_SOURCES  = test_os.c
test_os_LDADD    = ../suzerain/libsuzerain.la

TESTS                   += test_blas_et_al
check_PROGRAMS          += test_blas_et_al
test_blas_et_al_SOURCES  = test_blas_et_al.c
test_blas_et_al_LDADD    = ../suzerain/libsuzerain.la
test_blas_et_al_LDADD   += $(GSL_LIBS)
test_blas_et_al_LDADD   += $(MKL_LIBS)

TESTS            += test_mpl
check_PROGRAMS   += test_mpl
test_mpl_SOURCES  = test_mpl.cpp
test_mpl_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_mpl_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_allocator
check_PROGRAMS         += test_allocator
test_allocator_SOURCES  = test_allocator.cpp
test_allocator_LDADD    = ../suzerain/libsuzerain.la
test_allocator_LDADD   += $(MKL_LIBS)
test_allocator_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_bspline
check_PROGRAMS       += test_bspline
test_bspline_SOURCES  = test_bspline.cpp
test_bspline_LDADD    = ../suzerain/libsuzerain.la
test_bspline_LDADD   += $(GSL_LIBS)
test_bspline_LDADD   += $(MKL_LIBS)
test_bspline_LDADD   += $(LOG4CXX_LIBS)
test_bspline_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_bspline_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS             += test_math
check_PROGRAMS    += test_math
test_math_SOURCES  = test_math.cpp
test_math_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_math_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS               += test_complex
check_PROGRAMS      += test_complex
test_complex_SOURCES = test_complex.cpp
test_complex_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_complex_LDADD  += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                       += test_blas_et_al_cpp
check_PROGRAMS              += test_blas_et_al_cpp
test_blas_et_al_cpp_SOURCES  = test_blas_et_al_cpp.cpp
test_blas_et_al_cpp_LDADD    = ../suzerain/libsuzerain.la
test_blas_et_al_cpp_LDADD   += $(MKL_LIBS)
test_blas_et_al_cpp_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_blas_et_al_cpp_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_iterator
check_PROGRAMS       += test_iterator
test_iterator_SOURCES = test_iterator.cpp
test_iterator_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_iterator_LDADD  += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_functional
check_PROGRAMS         += test_functional
test_functional_SOURCES = test_functional.cpp
test_functional_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_functional_LDADD  += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                 += test_storage
check_PROGRAMS        += test_storage
test_storage_SOURCES   = test_storage.cpp
test_storage_LDADD     = ../suzerain/libsuzerain.la
test_storage_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_storage_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS               += test_state
check_PROGRAMS      += test_state
test_state_SOURCES   = test_state.cpp
test_state_LDADD     = ../suzerain/libsuzerain.la
test_state_LDADD    += $(MKL_LIBS)
test_state_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_state_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS               += test_pencil
check_PROGRAMS      += test_pencil
test_pencil_SOURCES  = test_pencil.cpp
test_pencil_LDADD    = ../suzerain/libsuzerain.la
test_pencil_LDADD   += $(P3DFFT_LIBS)
test_pencil_LDADD   += $(FFTW3_LIBS)
test_pencil_LDADD   += $(GSL_LIBS)
test_pencil_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencil_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                    += test_orthonormal
check_PROGRAMS           += test_orthonormal
test_orthonormal_SOURCES  = test_orthonormal.cpp
test_orthonormal_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_orthonormal_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
EXTRA_DIST               += test_orthonormal.sage

TESTS                   += test_richardson
check_PROGRAMS          += test_richardson
test_richardson_SOURCES  = test_richardson.c
test_richardson_LDADD    = ../suzerain/libsuzerain.la
test_richardson_LDADD   += $(GSL_LIBS)

TESTS                    += test_timestepper
check_PROGRAMS           += test_timestepper
test_timestepper_SOURCES  = test_timestepper.cpp
test_timestepper_LDADD    = ../suzerain/libsuzerain.la
test_timestepper_LDADD   += $(GSL_LIBS)
test_timestepper_LDADD   += $(MKL_LIBS)
test_timestepper_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_timestepper_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                    += test_multi_array
check_PROGRAMS           += test_multi_array
test_multi_array_SOURCES  = test_multi_array.cpp
test_multi_array_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_multi_array_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS             += test_fftw
check_PROGRAMS    += test_fftw
test_fftw_SOURCES  = test_fftw.cpp
test_fftw_CPPFLAGS = $(AM_CPPFLAGS) $(FFTW3_CFLAGS)
test_fftw_LDADD    = ../suzerain/libsuzerain.la
test_fftw_LDADD   += $(FFTW3_LIBS)
test_fftw_LDADD   += $(BOOST_PROGRAM_OPTIONS_LIBS)
test_fftw_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_fftw_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_fftw_1D
check_PROGRAMS       += test_fftw_1D
test_fftw_1D_SOURCES  = test_fftw_1D.cpp
test_fftw_1d_CPPFLAGS = $(AM_CPPFLAGS) $(FFTW3_CFLAGS)
test_fftw_1D_LDADD    = $(FFTW3_LIBS)
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                          += test_pencilfft_utility
check_PROGRAMS                 += test_pencilfft_utility
test_pencilfft_utility_SOURCES  = test_pencilfft_utility.cpp
test_pencilfft_utility_CPPFLAGS = $(AM_CPPFLAGS) $(FFTW3_CFLAGS)
test_pencilfft_utility_LDADD    = $(FFTW3_LIBS)
test_pencilfft_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_1D
check_PROGRAMS            += test_pencilfft_1D
test_pencilfft_1D_SOURCES  = test_pencilfft_1D.cpp
test_pencilfft_1D_CPPFLAGS = $(AM_CPPFLAGS) $(FFTW3_CFLAGS)
test_pencilfft_1D_LDADD    = $(FFTW3_LIBS)
test_pencilfft_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_2D
check_PROGRAMS            += test_pencilfft_2D
test_pencilfft_2D_SOURCES  = test_pencilfft_2D.cpp
test_pencilfft_2D_CPPFLAGS = $(AM_CPPFLAGS) $(FFTW3_CFLAGS)
test_pencilfft_2D_LDADD    = $(FFTW3_LIBS)
test_pencilfft_2D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_2D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_3D
check_PROGRAMS            += test_pencilfft_3D
test_pencilfft_3D_SOURCES  = test_pencilfft_3D.cpp
test_pencilfft_3D_CPPFLAGS = $(AM_CPPFLAGS) $(FFTW3_CFLAGS)
test_pencilfft_3D_LDADD    = $(FFTW3_LIBS)
test_pencilfft_3D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_3D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_underling.sh          # NB Wrapper script
EXTRA_DIST             += test_underling.sh          # NB Wrapper script
check_PROGRAMS         += test_underling
test_underling_SOURCES  = test_underling.cpp
test_underling_CPPFLAGS = $(AM_CPPFLAGS) $(FFTW3_THREADS_CFLAGS)
test_underling_LDADD    = ../suzerain/libunderling.la
test_underling_LDADD   += $(FFTW3_MPI_THREADS_LIBS)
test_underling_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_underling_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                      += test_underling_fft.sh  # NB Wrapper script
EXTRA_DIST                 += test_underling_fft.sh  # NB Wrapper script
check_PROGRAMS             += test_underling_fft
test_underling_fft_SOURCES  = test_underling_fft.cpp
test_underling_fft_CPPFLAGS = $(AM_CPPFLAGS) $(FFTW3_THREADS_CFLAGS)
test_underling_fft_LDADD    = ../suzerain/libunderling.la
test_underling_fft_LDADD   += $(FFTW3_MPI_THREADS_LIBS)
test_underling_fft_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_underling_fft_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
