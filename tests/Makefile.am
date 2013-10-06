## Process this file with automake to produce Makefile.in

AM_CPPFLAGS         = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS        += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS        += $(UNDERLING_CFLAGS) $(P3DFFT_CFLAGS) $(GRVY_CFLAGS)
AM_CPPFLAGS        += $(FFTW3_THREADS_CFLAGS) $(GSL_CFLAGS) $(MKL_CFLAGS)
AM_CPPFLAGS        += $(BOOST_CPPFLAGS) $(LOG4CXX_CFLAGS) $(EIGEN_CFLAGS)
AM_LDFLAGS          = # Append below
bin_PROGRAMS        = # Append below
check_PROGRAMS      = # Append below
EXTRA_PROGRAMS      = # Append below
dist_check_SCRIPTS  = # Append below
noinst_HEADERS      = # Append below
CLEANFILES          = # Append below
TESTS               = # Append below

noinst_HEADERS += test_tools.hpp

# MPI-based compilation and linking
CXX = @MPICXX@
CC  = @MPICC@

# Permit building EXTRA_PROGRAMS and check_PROGRAMS without running tests
extra: $(EXTRA_PROGRAMS) $(check_PROGRAMS)

# Specify test environment for 'make check'.
# One may override by manipulating $TESTS_ENVIRONMENT using the command line.
AM_TESTS_ENVIRONMENT  = export TMPDIR="$(abs_builddir)";
AM_TESTS_ENVIRONMENT += export LOG4CXX_CONFIGURATION="$(abs_srcdir)/log4cxx.properties";
AM_TESTS_FD_REDIRECT  = 9>&2 # Use descriptor 9 to get human attention

# Notice that all tests build for 'make all' target to be successful.
# Much of the C++ functionality is header-only; this ensures it builds.

# Suspect a problem similar to test_running_statistics if Intel complains
TESTS                  += test_exprparse
EXTRA_PROGRAMS         += test_exprparse
test_exprparse_SOURCES  = test_exprparse.cpp
test_exprparse_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_exprparse_LDADD    = ../suzerain/libexprparse.la
test_exprparse_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_exprparse_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS              += test_empty
EXTRA_PROGRAMS     += test_empty
test_empty_SOURCES  = test_empty.cpp
test_empty_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_empty_LDADD    = ../suzerain/libsuzerain.la
test_empty_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_empty_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS              += test_error
EXTRA_PROGRAMS     += test_error
test_error_SOURCES  = test_error.cpp
test_error_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_error_LDADD    = ../suzerain/libsuzerain.la
test_error_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_error_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

# TODO Logic is header-only, but Intel requires libsuzerain.la link on Lonestar
TESTS                           += test_running_statistics
EXTRA_PROGRAMS                  += test_running_statistics
test_running_statistics_SOURCES  = test_running_statistics.cpp
test_running_statistics_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_running_statistics_LDADD    = ../suzerain/libsuzerain.la
test_running_statistics_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_running_statistics_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                 += test_gbmatrix
EXTRA_PROGRAMS        += test_gbmatrix
test_gbmatrix_SOURCES  = test_gbmatrix.c
test_gbmatrix_LDADD    = ../suzerain/libsuzerain.la

TESTS           += test_os
EXTRA_PROGRAMS  += test_os
test_os_SOURCES  = test_os.c
test_os_LDADD    = ../suzerain/libsuzerain.la

TESTS               += test_parcel
EXTRA_PROGRAMS      += test_parcel
test_parcel_SOURCES  = test_parcel.c
test_parcel_LDADD    = ../suzerain/libsuzerain.la

TESTS                   += test_blas_et_al
EXTRA_PROGRAMS          += test_blas_et_al
test_blas_et_al_SOURCES  = test_blas_et_al.c
test_blas_et_al_LDADD    = ../suzerain/libsuzerain.la

TESTS                             += test_zgbsv_specification
EXTRA_PROGRAMS                    += test_zgbsv_specification
test_zgbsv_specification_SOURCES   = test_zgbsv_specification.cpp
test_zgbsv_specification_CPPFLAGS  = $(AM_CPPFLAGS)
test_zgbsv_specification_CPPFLAGS += -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_zgbsv_specification_LDADD     = ../suzerain/libsuzerain.la
test_zgbsv_specification_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_zgbsv_specification_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS             += test_gbmv
EXTRA_PROGRAMS    += test_gbmv
test_gbmv_SOURCES  = test_gbmv.cpp
test_gbmv_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_gbmv_LDADD    = ../suzerain/libsuzerain.la
test_gbmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_gbmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS              += test_gbdmv
EXTRA_PROGRAMS     += test_gbdmv
test_gbdmv_SOURCES  = test_gbdmv.cpp
test_gbdmv_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_gbdmv_LDADD    = ../suzerain/libsuzerain.la
test_gbdmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_gbdmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS               += test_gbddmv
EXTRA_PROGRAMS      += test_gbddmv
test_gbddmv_SOURCES  = test_gbddmv.cpp
test_gbddmv_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_gbddmv_LDADD    = ../suzerain/libsuzerain.la
test_gbddmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_gbddmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_gbdddmv
EXTRA_PROGRAMS       += test_gbdddmv
test_gbdddmv_SOURCES  = test_gbdddmv.cpp
test_gbdddmv_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_gbdddmv_LDADD    = ../suzerain/libsuzerain.la
test_gbdddmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_gbdddmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                 += test_gbddddmv
EXTRA_PROGRAMS        += test_gbddddmv
test_gbddddmv_SOURCES  = test_gbddddmv.cpp
test_gbddddmv_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_gbddddmv_LDADD    = ../suzerain/libsuzerain.la
test_gbddddmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_gbddddmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_gbdddddmv
EXTRA_PROGRAMS         += test_gbdddddmv
test_gbdddddmv_SOURCES  = test_gbdddddmv.cpp
test_gbdddddmv_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_gbdddddmv_LDADD    = ../suzerain/libsuzerain.la
test_gbdddddmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_gbdddddmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS             += test_sbmv
EXTRA_PROGRAMS    += test_sbmv
test_sbmv_SOURCES  = test_sbmv.cpp
test_sbmv_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_sbmv_LDADD    = ../suzerain/libsuzerain.la
test_sbmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_sbmv_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS               += test_bsmbsm
EXTRA_PROGRAMS      += test_bsmbsm
test_bsmbsm_SOURCES  = test_bsmbsm.cpp
test_bsmbsm_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_bsmbsm_LDADD    = ../suzerain/libsuzerain.la
test_bsmbsm_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_bsmbsm_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_htstretch
EXTRA_PROGRAMS         += test_htstretch
test_htstretch_SOURCES  = test_htstretch.c
test_htstretch_LDADD    = ../suzerain/libsuzerain.la

TESTS                      += test_radial_nozzle
EXTRA_PROGRAMS             += test_radial_nozzle
test_radial_nozzle_SOURCES  = test_radial_nozzle.c
test_radial_nozzle_LDADD    = ../suzerain/libsuzerain.la

TESTS                  += test_rngstream
EXTRA_PROGRAMS         += test_rngstream
test_rngstream_SOURCES  = test_rngstream.cpp
test_rngstream_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_rngstream_LDADD    = ../suzerain/libsuzerain.la
test_rngstream_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rngstream_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS            += test_mpl
EXTRA_PROGRAMS   += test_mpl
test_mpl_SOURCES  = test_mpl.cpp
test_mpl_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_mpl_LDADD    = ../suzerain/libsuzerain.la
test_mpl_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_mpl_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

### Compilation footprint of test_allocator is large even at -O0
### so we break it up into test_allocator{0,1,2,3,4,5,6,7} below.
## TESTS                   += test_allocator
## EXTRA_PROGRAMS          += test_allocator
## test_allocator_SOURCES   = test_allocator.cpp
## test_allocator_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
## test_allocator_LDADD     = ../suzerain/libsuzerain.la
## test_allocator_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
## test_allocator_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator0
EXTRA_PROGRAMS            += test_allocator0
test_allocator0_SOURCES    = test_allocator.cpp
test_allocator0_CPPFLAGS   = $(AM_CPPFLAGS) -DSTL_ALLOCATOR_TEST_PART=0
test_allocator0_CPPFLAGS  += -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_allocator0_LDADD      = ../suzerain/libsuzerain.la
test_allocator0_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator0_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator1
EXTRA_PROGRAMS            += test_allocator1
test_allocator1_SOURCES    = test_allocator.cpp
test_allocator1_CPPFLAGS   = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_allocator1_CPPFLAGS  += -DSTL_ALLOCATOR_TEST_PART=1
test_allocator1_LDADD      = ../suzerain/libsuzerain.la
test_allocator1_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator1_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator2
EXTRA_PROGRAMS            += test_allocator2
test_allocator2_SOURCES    = test_allocator.cpp
test_allocator2_CPPFLAGS   = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_allocator2_CPPFLAGS  += -DSTL_ALLOCATOR_TEST_PART=2
test_allocator2_LDADD      = ../suzerain/libsuzerain.la
test_allocator2_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator2_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator3
EXTRA_PROGRAMS            += test_allocator3
test_allocator3_SOURCES    = test_allocator.cpp
test_allocator3_CPPFLAGS   = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_allocator3_CPPFLAGS  += -DSTL_ALLOCATOR_TEST_PART=3
test_allocator3_LDADD      = ../suzerain/libsuzerain.la
test_allocator3_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator3_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator4
EXTRA_PROGRAMS            += test_allocator4
test_allocator4_SOURCES    = test_allocator.cpp
test_allocator4_CPPFLAGS   = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_allocator4_CPPFLAGS  += -DSTL_ALLOCATOR_TEST_PART=4
test_allocator4_LDADD      = ../suzerain/libsuzerain.la
test_allocator4_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator4_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator5
EXTRA_PROGRAMS            += test_allocator5
test_allocator5_SOURCES    = test_allocator.cpp
test_allocator5_CPPFLAGS   = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_allocator5_CPPFLAGS  += -DSTL_ALLOCATOR_TEST_PART=5
test_allocator5_LDADD      = ../suzerain/libsuzerain.la
test_allocator5_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator5_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator6
EXTRA_PROGRAMS            += test_allocator6
test_allocator6_SOURCES    = test_allocator.cpp
test_allocator6_CPPFLAGS   = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_allocator6_CPPFLAGS  += -DSTL_ALLOCATOR_TEST_PART=6
test_allocator6_LDADD      = ../suzerain/libsuzerain.la
test_allocator6_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator6_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_allocator7
EXTRA_PROGRAMS            += test_allocator7
test_allocator7_SOURCES    = test_allocator.cpp
test_allocator7_CPPFLAGS   = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_allocator7_CPPFLAGS  += -DSTL_ALLOCATOR_TEST_PART=7
test_allocator7_LDADD      = ../suzerain/libsuzerain.la
test_allocator7_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_allocator7_LDADD     += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_shared_range
EXTRA_PROGRAMS            += test_shared_range
test_shared_range_SOURCES  = test_shared_range.cpp
test_shared_range_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_shared_range_LDADD    = ../suzerain/libsuzerain.la
test_shared_range_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_shared_range_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                        += test_coalescing_pool
EXTRA_PROGRAMS               += test_coalescing_pool
test_coalescing_pool_SOURCES  = test_coalescing_pool.cpp
test_coalescing_pool_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_coalescing_pool_LDADD    = ../suzerain/libsuzerain.la
test_coalescing_pool_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_coalescing_pool_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_bspline
EXTRA_PROGRAMS       += test_bspline
test_bspline_SOURCES  = test_bspline.c
test_bspline_LDADD    = ../suzerain/libsuzerain.la

TESTS                   += test_bsplineop
EXTRA_PROGRAMS          += test_bsplineop
test_bsplineop_SOURCES   = test_bsplineop.cpp
test_bsplineop_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_bsplineop_LDADD     = ../suzerain/libsuzerain.la
test_bsplineop_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_bsplineop_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS            += test_bl
EXTRA_PROGRAMS   += test_bl
test_bl_SOURCES   = test_bl.cpp
test_bl_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_bl_LDADD     = ../suzerain/libsuzerain.la
test_bl_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_bl_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_filterop
EXTRA_PROGRAMS         += test_filterop
test_filterop_SOURCES   = test_filterop.cpp
test_filterop_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_filterop_LDADD     = ../suzerain/libsuzerain.la
test_filterop_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_filterop_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS             += test_math
EXTRA_PROGRAMS    += test_math
test_math_SOURCES  = test_math.cpp
test_math_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_math_LDADD    = ../suzerain/libsuzerain.la
test_math_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_math_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_utility
EXTRA_PROGRAMS       += test_utility
test_utility_SOURCES  = test_utility.cpp
test_utility_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_utility_LDADD    = ../suzerain/libsuzerain.la
test_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_complex
EXTRA_PROGRAMS       += test_complex
test_complex_SOURCES  = test_complex.cpp
test_complex_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_complex_LDADD    = ../suzerain/libsuzerain.la
test_complex_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_complex_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS              += test_kahan
EXTRA_PROGRAMS     += test_kahan
test_kahan_SOURCES  = test_kahan.cpp
test_kahan_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_kahan_LDADD    = ../suzerain/libsuzerain.la
test_kahan_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_kahan_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_traits
EXTRA_PROGRAMS       += test_traits
test_traits_SOURCES   = test_traits.cpp
test_traits_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_traits_LDADD     = ../suzerain/libsuzerain.la
test_traits_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_traits_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                      += test_mpi_datatype
EXTRA_PROGRAMS             += test_mpi_datatype
test_mpi_datatype_SOURCES   = test_mpi_datatype.cpp
test_mpi_datatype_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_mpi_datatype_LDADD     = ../suzerain/libsuzerain.la
test_mpi_datatype_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_mpi_datatype_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                       += test_blas_et_al_cpp
EXTRA_PROGRAMS              += test_blas_et_al_cpp
test_blas_et_al_cpp_SOURCES  = test_blas_et_al_cpp.cpp
test_blas_et_al_cpp_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_blas_et_al_cpp_LDADD    = ../suzerain/libsuzerain.la
test_blas_et_al_cpp_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_blas_et_al_cpp_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                 += test_iterator
EXTRA_PROGRAMS        += test_iterator
test_iterator_SOURCES  = test_iterator.cpp
test_iterator_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_iterator_LDADD    = ../suzerain/libsuzerain.la
test_iterator_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_iterator_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                   += test_functional
EXTRA_PROGRAMS          += test_functional
test_functional_SOURCES  = test_functional.cpp
test_functional_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_functional_LDADD    = ../suzerain/libsuzerain.la
test_functional_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_functional_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                 += test_storage
EXTRA_PROGRAMS        += test_storage
test_storage_SOURCES   = test_storage.cpp
test_storage_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_storage_LDADD     = ../suzerain/libsuzerain.la
test_storage_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_storage_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS               += test_state
EXTRA_PROGRAMS      += test_state
test_state_SOURCES   = test_state.cpp
test_state_CPPFLAGS  = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_state_LDADD     = ../suzerain/libsuzerain.la
test_state_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_state_LDADD    += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS               += test_pencil
EXTRA_PROGRAMS      += test_pencil
test_pencil_SOURCES  = test_pencil.cpp
test_pencil_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_pencil_LDADD    = ../suzerain/libsuzerain.la
test_pencil_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencil_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS              += test_rholt
EXTRA_PROGRAMS     += test_rholt
test_rholt_SOURCES  = test_rholt.cpp
test_rholt_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_rholt_LDADD    = ../suzerain/libsuzerain.la
test_rholt_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rholt_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
dist_check_SCRIPTS += test_rholt.sage

TESTS               += test_rholut
EXTRA_PROGRAMS      += test_rholut
test_rholut_SOURCES  = test_rholut.cpp
test_rholut_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_rholut_LDADD    = ../suzerain/libsuzerain.la
test_rholut_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rholut_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
dist_check_SCRIPTS  += test_rholut.sage

TESTS                      += test_rholut_imexop
EXTRA_PROGRAMS             += test_rholut_imexop
test_rholut_imexop_SOURCES  = test_rholut_imexop.cpp
test_rholut_imexop_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_rholut_imexop_LDADD    = ../suzerain/libsuzerain.la
test_rholut_imexop_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rholut_imexop_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                        += test_rholut_imexop00
EXTRA_PROGRAMS               += test_rholut_imexop00
test_rholut_imexop00_SOURCES  = test_rholut_imexop00.cpp
test_rholut_imexop00_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_rholut_imexop00_LDADD    = ../suzerain/libsuzerain.la
test_rholut_imexop00_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rholut_imexop00_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                        += test_reacting_imexop
EXTRA_PROGRAMS               += test_reacting_imexop
test_reacting_imexop_SOURCES  = test_reacting_imexop.cpp
test_reacting_imexop_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_reacting_imexop_LDADD    = ../suzerain/libsuzerain.la
test_reacting_imexop_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_reacting_imexop_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                  += test_richardson
EXTRA_PROGRAMS         += test_richardson
test_richardson_SOURCES = test_richardson.c
test_richardson_LDADD   = ../suzerain/libsuzerain.la

TESTS                    += test_multi_array
EXTRA_PROGRAMS           += test_multi_array
test_multi_array_SOURCES  = test_multi_array.cpp
test_multi_array_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_multi_array_LDADD    = ../suzerain/libsuzerain.la
test_multi_array_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_multi_array_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                       += test_timecontroller
EXTRA_PROGRAMS              += test_timecontroller
test_timecontroller_SOURCES  = test_timecontroller.cpp
test_timecontroller_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_timecontroller_LDADD    = ../suzerain/libsuzerain.la
test_timecontroller_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_timecontroller_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                   += test_lowstorage
EXTRA_PROGRAMS          += test_lowstorage
test_lowstorage_SOURCES  = test_lowstorage.cpp
test_lowstorage_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_lowstorage_LDADD    = ../suzerain/libsuzerain.la
test_lowstorage_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_lowstorage_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
dist_check_SCRIPTS      += test_lowstorage.sage

TESTS             += test_fftw
EXTRA_PROGRAMS    += test_fftw
test_fftw_SOURCES  = test_fftw.cpp
test_fftw_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_fftw_LDADD    = ../suzerain/libsuzerain.la
test_fftw_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_fftw_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_fftw_1D
EXTRA_PROGRAMS       += test_fftw_1D
test_fftw_1D_SOURCES  = test_fftw_1D.cpp
test_fftw_1D_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_fftw_1D_LDADD    = ../suzerain/libsuzerain.la
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_fftw_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                          += test_pencilfft_utility
EXTRA_PROGRAMS                 += test_pencilfft_utility
test_pencilfft_utility_SOURCES  = test_pencilfft_utility.cpp
test_pencilfft_utility_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_pencilfft_utility_LDADD    = ../suzerain/libsuzerain.la
test_pencilfft_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_utility_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_1D
EXTRA_PROGRAMS            += test_pencilfft_1D
test_pencilfft_1D_SOURCES  = test_pencilfft_1D.cpp
test_pencilfft_1D_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_pencilfft_1D_LDADD    = ../suzerain/libsuzerain.la
test_pencilfft_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_1D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_2D
EXTRA_PROGRAMS            += test_pencilfft_2D
test_pencilfft_2D_SOURCES  = test_pencilfft_2D.cpp
test_pencilfft_2D_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_pencilfft_2D_LDADD    = ../suzerain/libsuzerain.la
test_pencilfft_2D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_2D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                     += test_pencilfft_3D
EXTRA_PROGRAMS            += test_pencilfft_3D
test_pencilfft_3D_SOURCES  = test_pencilfft_3D.cpp
test_pencilfft_3D_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_pencilfft_3D_LDADD    = ../suzerain/libsuzerain.la
test_pencilfft_3D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_pencilfft_3D_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                += test_inorder
EXTRA_PROGRAMS       += test_inorder
test_inorder_SOURCES  = test_inorder.cpp
test_inorder_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_inorder_LDADD    = ../suzerain/libsuzerain.la
test_inorder_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_inorder_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS              += test_constraint_treatment
EXTRA_PROGRAMS     += test_constraint_treatment
test_constraint_treatment_SOURCES  = test_constraint_treatment.cpp
test_constraint_treatment_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_constraint_treatment_LDADD    = ../suzerain/libsuzerain.la
test_constraint_treatment_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_constraint_treatment_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS                 += test_diffwave
EXTRA_PROGRAMS        += test_diffwave
test_diffwave_SOURCES  = test_diffwave.cpp
test_diffwave_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_diffwave_LDADD    = ../suzerain/libsuzerain.la
test_diffwave_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_diffwave_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

# Filed under check_PROGRAMS so binary is built before test is executed
# There must be a better way to express this dependency, but this works
TESTS                        += test_diffwave_p3dfft.sh
dist_check_SCRIPTS           += test_diffwave_p3dfft.sh
check_PROGRAMS               += test_diffwave_p3dfft
test_diffwave_p3dfft_SOURCES  = test_diffwave_p3dfft.cpp
test_diffwave_p3dfft_CPPFLAGS = $(AM_CPPFLAGS) -DBOOST_TEST_DYN_LINK -DBOOST_TEST_NO_LIB
test_diffwave_p3dfft_LDADD    = ../suzerain/libsuzerain.la
test_diffwave_p3dfft_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_diffwave_p3dfft_LDADD   += $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

# Filed under check_PROGRAMS so binary is built before test is executed
# There must be a better way to express this dependency, but this works
TESTS              += test_l2xz.sh
dist_check_SCRIPTS += test_l2xz.sh
check_PROGRAMS     += test_l2xz
test_l2xz_SOURCES   = test_l2xz.cpp
test_l2xz_LDADD     = ../suzerain/libsuzerain.la

# Clean up from Intel Compiler's -prof-gen output
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi

# Clean up EXTRA_PROGRAMS per hints in "Simple Tests using ‘parallel-tests’"
CLEANFILES += $(EXTRA_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
