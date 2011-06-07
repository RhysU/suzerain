## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(P3DFFT_CFLAGS) $(FFTW3_THREADS_CFLAGS)
AM_CPPFLAGS   += $(GSL_CFLAGS) $(MKL_CFLAGS)
AM_CPPFLAGS   += $(BOOST_CPPFLAGS) $(LOG4CXX_CFLAGS) $(EIGEN_CFLAGS)
AM_LDFLAGS     = # Append below
bin_PROGRAMS   = # Append below
check_PROGRAMS = # Append below
noinst_HEADERS = # Append below
CLEANFILES     = # Append below
EXTRA_DIST     = # Append below
TESTS          = # Append below

noinst_HEADERS += test_tools.hpp

# MPI-based compilation and linking
CXX = @MPICXX@
CC  = @MPICC@

# Provide a test-ready logging configuration
TESTS_ENVIRONMENT = LOG4CXX_CONFIGURATION=$(srcdir)/log4cxx.properties

# Demand that all tests build for 'make all' target to be successful.
# Much of the C++ functionality is header-only; this ensures it builds.
all-local: $(check_PROGRAMS)

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

TESTS                  += test_htstretch
check_PROGRAMS         += test_htstretch
test_htstretch_SOURCES  = test_htstretch.c
test_htstretch_LDADD    = ../suzerain/libsuzerain.la

TESTS            += test_mpl
check_PROGRAMS   += test_mpl
test_mpl_SOURCES  = test_mpl.cpp
test_mpl_LDADD    = ../suzerain/libsuzerain.la
test_mpl_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_mpl_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_allocator
check_PROGRAMS         += test_allocator
test_allocator_SOURCES  = test_allocator.cpp
test_allocator_LDADD    = ../suzerain/libsuzerain.la
test_allocator_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_bspline
check_PROGRAMS       += test_bspline
test_bspline_SOURCES  = test_bspline.c
test_bspline_LDADD    = ../suzerain/libsuzerain.la

TESTS                   += test_bsplineop
check_PROGRAMS          += test_bsplineop
test_bsplineop_SOURCES   = test_bsplineop.cpp
test_bsplineop_LDADD     = ../suzerain/libsuzerain.la
test_bsplineop_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_bsplineop_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS             += test_math
check_PROGRAMS    += test_math
test_math_SOURCES  = test_math.cpp
test_math_LDADD    = ../suzerain/libsuzerain.la
test_math_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_math_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_utility
check_PROGRAMS       += test_utility
test_utility_SOURCES  = test_utility.cpp
test_utility_LDADD    = ../suzerain/libsuzerain.la
test_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_complex
check_PROGRAMS       += test_complex
test_complex_SOURCES  = test_complex.cpp
test_complex_LDADD    = ../suzerain/libsuzerain.la
test_complex_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_complex_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_traits
check_PROGRAMS       += test_traits
test_traits_SOURCES   = test_traits.cpp
test_traits_LDADD     = ../suzerain/libsuzerain.la
test_traits_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_traits_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                      += test_mpi_datatype
check_PROGRAMS             += test_mpi_datatype
test_mpi_datatype_SOURCES   = test_mpi_datatype.cpp
test_mpi_datatype_LDADD     = ../suzerain/libsuzerain.la
test_mpi_datatype_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_mpi_datatype_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                       += test_blas_et_al_cpp
check_PROGRAMS              += test_blas_et_al_cpp
test_blas_et_al_cpp_SOURCES  = test_blas_et_al_cpp.cpp
test_blas_et_al_cpp_LDADD    = ../suzerain/libsuzerain.la
test_blas_et_al_cpp_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_blas_et_al_cpp_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_iterator
check_PROGRAMS       += test_iterator
test_iterator_SOURCES = test_iterator.cpp
test_iterator_LDADD   = ../suzerain/libsuzerain.la
test_iterator_LDADD  += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_iterator_LDADD  += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_functional
check_PROGRAMS         += test_functional
test_functional_SOURCES = test_functional.cpp
test_functional_LDADD   = ../suzerain/libsuzerain.la
test_functional_LDADD  += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
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
test_state_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_state_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS               += test_pencil
check_PROGRAMS      += test_pencil
test_pencil_SOURCES  = test_pencil.cpp
test_pencil_LDADD    = ../suzerain/libsuzerain.la
test_pencil_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencil_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS              += test_rholt
check_PROGRAMS     += test_rholt
test_rholt_SOURCES  = test_rholt.cpp
test_rholt_LDADD    = ../suzerain/libsuzerain.la
test_rholt_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rholt_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
EXTRA_DIST         += test_rholt.sage

TESTS               += test_rholut
check_PROGRAMS      += test_rholut
test_rholut_SOURCES  = test_rholut.cpp
test_rholut_LDADD    = ../suzerain/libsuzerain.la
test_rholut_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rholut_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
EXTRA_DIST          += test_rholut.sage

TESTS                   += test_richardson
check_PROGRAMS          += test_richardson
test_richardson_SOURCES  = test_richardson.c
test_richardson_LDADD    = ../suzerain/libsuzerain.la

TESTS                       += test_timecontroller
check_PROGRAMS              += test_timecontroller
test_timecontroller_SOURCES  = test_timecontroller.cpp
test_timecontroller_LDADD    = ../suzerain/libsuzerain.la
test_timecontroller_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_timecontroller_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                    += test_timestepper
check_PROGRAMS           += test_timestepper
test_timestepper_SOURCES  = test_timestepper.cpp
test_timestepper_LDADD    = ../suzerain/libsuzerain.la
test_timestepper_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_timestepper_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                    += test_multi_array
check_PROGRAMS           += test_multi_array
test_multi_array_SOURCES  = test_multi_array.cpp
test_multi_array_LDADD    = ../suzerain/libsuzerain.la
test_multi_array_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_multi_array_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS             += test_fftw
check_PROGRAMS    += test_fftw
test_fftw_SOURCES  = test_fftw.cpp
test_fftw_LDADD    = ../suzerain/libsuzerain.la
test_fftw_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_fftw_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_fftw_1D
check_PROGRAMS       += test_fftw_1D
test_fftw_1D_SOURCES  = test_fftw_1D.cpp
test_fftw_1D_LDADD    = ../suzerain/libsuzerain.la
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                          += test_pencilfft_utility
check_PROGRAMS                 += test_pencilfft_utility
test_pencilfft_utility_SOURCES  = test_pencilfft_utility.cpp
test_pencilfft_utility_LDADD    = ../suzerain/libsuzerain.la
test_pencilfft_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_1D
check_PROGRAMS            += test_pencilfft_1D
test_pencilfft_1D_SOURCES  = test_pencilfft_1D.cpp
test_pencilfft_1D_LDADD    = ../suzerain/libsuzerain.la
test_pencilfft_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_2D
check_PROGRAMS            += test_pencilfft_2D
test_pencilfft_2D_SOURCES  = test_pencilfft_2D.cpp
test_pencilfft_2D_LDADD    = ../suzerain/libsuzerain.la
test_pencilfft_2D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_2D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_3D
check_PROGRAMS            += test_pencilfft_3D
test_pencilfft_3D_SOURCES  = test_pencilfft_3D.cpp
test_pencilfft_3D_LDADD    = ../suzerain/libsuzerain.la
test_pencilfft_3D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_3D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                 += test_inorder
check_PROGRAMS        += test_inorder
test_inorder_SOURCES   = test_inorder.cpp
test_inorder_LDADD     = ../suzerain/libsuzerain.la
test_inorder_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_inorder_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                 += test_diffwave
check_PROGRAMS        += test_diffwave
test_diffwave_SOURCES  = test_diffwave.cpp
test_diffwave_LDADD    = ../suzerain/libsuzerain.la
test_diffwave_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_diffwave_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                         += test_diffwave_p3dfft.sh  # NB Wrapper script
EXTRA_DIST                    += test_diffwave_p3dfft.sh  # NB Wrapper script
check_PROGRAMS                += test_diffwave_p3dfft
test_diffwave_p3dfft_SOURCES   = test_diffwave_p3dfft.cpp
test_diffwave_p3dfft_LDADD     = ../suzerain/libsuzerain.la
test_diffwave_p3dfft_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_diffwave_p3dfft_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

# Clean up from Intel Compiler's -prof-gen output and precomiled headers
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
CLEANFILES += *.pchi

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
