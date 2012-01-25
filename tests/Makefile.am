## Process this file with automake to produce Makefile.in

AM_CPPFLAGS         = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS        += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS        += $(UNDERLING_CFLAGS) $(P3DFFT_CFLAGS)
AM_CPPFLAGS        += $(FFTW3_THREADS_CFLAGS) $(GSL_CFLAGS) $(MKL_CFLAGS)
AM_CPPFLAGS        += $(BOOST_CPPFLAGS) $(LOG4CXX_CFLAGS) $(EIGEN_CFLAGS)
AM_LDFLAGS          = # Append below
bin_PROGRAMS        = # Append below
check_PROGRAMS      = # Append below
dist_check_SCRIPTS  = # Append below
noinst_HEADERS      = # Append below
CLEANFILES          = # Append below
TESTS               = # Append below

noinst_HEADERS += test_tools.hpp

# MPI-based compilation and linking
CXX = @MPICXX@
CC  = @MPICC@

# Specify test environment for 'make check'
TESTS_ENVIRONMENT  = TMPDIR=$(abs_builddir)
TESTS_ENVIRONMENT += LOG4CXX_CONFIGURATION=$(abs_srcdir)/log4cxx.properties

# Demand that all tests build for 'make all' target to be successful.
# Much of the C++ functionality is header-only; this ensures it builds.
all-local: $(check_PROGRAMS)

TESTS              += test_empty
check_PROGRAMS     += test_empty
test_empty_SOURCES  = test_empty.cpp
test_empty_LDADD    = ../suzerain/libsuzerain.la
test_empty_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_empty_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

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

TESTS               += test_parcel
check_PROGRAMS      += test_parcel
test_parcel_SOURCES  = test_parcel.cpp
test_parcel_LDADD    = ../suzerain/libsuzerain.la

TESTS                   += test_blas_et_al
check_PROGRAMS          += test_blas_et_al
test_blas_et_al_SOURCES  = test_blas_et_al.c
test_blas_et_al_LDADD    = ../suzerain/libsuzerain.la

TESTS             += test_gbmv
check_PROGRAMS    += test_gbmv
test_gbmv_SOURCES  = test_gbmv.cpp
test_gbmv_LDADD    = ../suzerain/libsuzerain.la

TESTS              += test_gbdmv
check_PROGRAMS     += test_gbdmv
test_gbdmv_SOURCES  = test_gbdmv.cpp
test_gbdmv_LDADD    = ../suzerain/libsuzerain.la

TESTS               += test_gbddmv
check_PROGRAMS      += test_gbddmv
test_gbddmv_SOURCES  = test_gbddmv.cpp
test_gbddmv_LDADD    = ../suzerain/libsuzerain.la

TESTS               += test_gbidmv
check_PROGRAMS      += test_gbidmv
test_gbidmv_SOURCES  = test_gbidmv.cpp
test_gbidmv_LDADD    = ../suzerain/libsuzerain.la

TESTS                += test_gbdddmv
check_PROGRAMS       += test_gbdddmv
test_gbdddmv_SOURCES  = test_gbdddmv.cpp
test_gbdddmv_LDADD    = ../suzerain/libsuzerain.la

TESTS                += test_gbiddmv
check_PROGRAMS       += test_gbiddmv
test_gbiddmv_SOURCES  = test_gbiddmv.cpp
test_gbiddmv_LDADD    = ../suzerain/libsuzerain.la

TESTS             += test_sbmv
check_PROGRAMS    += test_sbmv
test_sbmv_SOURCES  = test_sbmv.cpp
test_sbmv_LDADD    = ../suzerain/libsuzerain.la

TESTS               += test_bsmbsm
check_PROGRAMS      += test_bsmbsm
test_bsmbsm_SOURCES  = test_bsmbsm.c
test_bsmbsm_LDADD    = ../suzerain/libsuzerain.la

TESTS                  += test_htstretch
check_PROGRAMS         += test_htstretch
test_htstretch_SOURCES  = test_htstretch.c
test_htstretch_LDADD    = ../suzerain/libsuzerain.la

TESTS                  += test_RngStream
check_PROGRAMS         += test_RngStream
test_RngStream_SOURCES  = test_RngStream.cpp
test_RngStream_LDADD    = ../suzerain/libsuzerain.la
test_RngStream_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_RngStream_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS            += test_mpl
check_PROGRAMS   += test_mpl
test_mpl_SOURCES  = test_mpl.cpp
test_mpl_LDADD    = ../suzerain/libsuzerain.la
test_mpl_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_mpl_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

### Compilation footprint of test_allocator is large even at -O0
### so we break it up into test_allocator{0,1,2,3,4,5,6,7} below.
## TESTS                   += test_allocator
## check_PROGRAMS          += test_allocator
## test_allocator_SOURCES  = test_allocator.cpp
## test_allocator_LDADD    = ../suzerain/libsuzerain.la
## test_allocator_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
## test_allocator_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator0
check_PROGRAMS            += test_allocator0
test_allocator0_SOURCES    = test_allocator.cpp
test_allocator0_LDADD      = ../suzerain/libsuzerain.la
test_allocator0_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator0_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_allocator0_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=0

TESTS                     += test_allocator1
check_PROGRAMS            += test_allocator1
test_allocator1_SOURCES    = test_allocator.cpp
test_allocator1_LDADD      = ../suzerain/libsuzerain.la
test_allocator1_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator1_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_allocator1_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=1

TESTS                     += test_allocator2
check_PROGRAMS            += test_allocator2
test_allocator2_SOURCES    = test_allocator.cpp
test_allocator2_LDADD      = ../suzerain/libsuzerain.la
test_allocator2_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator2_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_allocator2_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=2

TESTS                     += test_allocator3
check_PROGRAMS            += test_allocator3
test_allocator3_SOURCES    = test_allocator.cpp
test_allocator3_LDADD      = ../suzerain/libsuzerain.la
test_allocator3_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator3_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_allocator3_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=3

TESTS                     += test_allocator4
check_PROGRAMS            += test_allocator4
test_allocator4_SOURCES    = test_allocator.cpp
test_allocator4_LDADD      = ../suzerain/libsuzerain.la
test_allocator4_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator4_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_allocator4_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=4

TESTS                     += test_allocator5
check_PROGRAMS            += test_allocator5
test_allocator5_SOURCES    = test_allocator.cpp
test_allocator5_LDADD      = ../suzerain/libsuzerain.la
test_allocator5_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator5_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_allocator5_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=5

TESTS                     += test_allocator6
check_PROGRAMS            += test_allocator6
test_allocator6_SOURCES    = test_allocator.cpp
test_allocator6_LDADD      = ../suzerain/libsuzerain.la
test_allocator6_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator6_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_allocator6_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=6

TESTS                     += test_allocator7
check_PROGRAMS            += test_allocator7
test_allocator7_SOURCES    = test_allocator.cpp
test_allocator7_LDADD      = ../suzerain/libsuzerain.la
test_allocator7_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator7_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
test_allocator7_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=7

TESTS                     += test_shared_range
check_PROGRAMS            += test_shared_range
test_shared_range_SOURCES  = test_shared_range.cpp
test_shared_range_LDADD    = ../suzerain/libsuzerain.la
test_shared_range_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_shared_range_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                        += test_coalescing_pool
check_PROGRAMS               += test_coalescing_pool
test_coalescing_pool_SOURCES  = test_coalescing_pool.cpp
test_coalescing_pool_LDADD    = ../suzerain/libsuzerain.la
test_coalescing_pool_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_coalescing_pool_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_exprparse
check_PROGRAMS         += test_exprparse
test_exprparse_SOURCES  = test_exprparse.cpp
test_exprparse_LDADD    = ../suzerain/libsuzerain.la
test_exprparse_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_exprparse_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

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
dist_check_SCRIPTS += test_rholt.sage

TESTS               += test_rholut
check_PROGRAMS      += test_rholut
test_rholut_SOURCES  = test_rholut.cpp
test_rholut_LDADD    = ../suzerain/libsuzerain.la
test_rholut_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rholut_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
dist_check_SCRIPTS  += test_rholut.sage

TESTS                   += test_richardson
check_PROGRAMS          += test_richardson
test_richardson_SOURCES  = test_richardson.c
test_richardson_LDADD    = ../suzerain/libsuzerain.la

TESTS                    += test_multi_array
check_PROGRAMS           += test_multi_array
test_multi_array_SOURCES  = test_multi_array.cpp
test_multi_array_LDADD    = ../suzerain/libsuzerain.la
test_multi_array_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_multi_array_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

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

TESTS                         += test_diffwave_p3dfft.sh
dist_check_SCRIPTS            += test_diffwave_p3dfft.sh
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
