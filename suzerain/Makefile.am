## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(MKL_CFLAGS) $(GSL_CFLAGS) $(LOG4CXX_CFLAGS)
AM_CPPFLAGS   += $(BOOST_CPPFLAGS) $(P3DFFT_CFLAGS) $(FFTW3_CFLAGS)
AM_CPPFLAGS   += $(EIGEN_CFLAGS)
AM_LDFLAGS     = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
AM_LDFLAGS    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
bin_PROGRAMS       = # Append below
check_PROGRAMS     = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
pkginclude_HEADERS = # Append below
CLEANFILES         = # Append below

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@
CC=@MPICC@

# Clean up from Intel Compiler's -prof-gen output
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
# Clean up from Intel Compiler's precompiled header output
CLEANFILES += *.pchi

BUILT_SOURCES = suzerain-config.h

pkginclude_HEADERS += allocator.hpp
pkginclude_HEADERS += blas_et_al.h
pkginclude_HEADERS += blas_et_al.hpp
pkginclude_HEADERS += bspline_definition.hpp
pkginclude_HEADERS += bspline.h
pkginclude_HEADERS += bspline.hpp
pkginclude_HEADERS += common.h
pkginclude_HEADERS += common.hpp
pkginclude_HEADERS += complex.hpp
pkginclude_HEADERS += diffwave.h
pkginclude_HEADERS += diffwave.hpp
pkginclude_HEADERS += error.h
pkginclude_HEADERS += fftw.hpp
pkginclude_HEADERS += functional.hpp
pkginclude_HEADERS += function.h
pkginclude_HEADERS += gbmatrix.h
pkginclude_HEADERS += grid_definition.hpp
pkginclude_HEADERS += htstretch.h
pkginclude_HEADERS += iterator.hpp
pkginclude_HEADERS += math.hpp
pkginclude_HEADERS += mpi.h
pkginclude_HEADERS += mpi.hpp
pkginclude_HEADERS += mpi_datatype.hpp
pkginclude_HEADERS += mpl.hpp
pkginclude_HEADERS += multi_array.hpp
pkginclude_HEADERS += orthonormal.hpp
pkginclude_HEADERS += os.h
pkginclude_HEADERS += pencilfft.hpp
pkginclude_HEADERS += pencil_grid.hpp
pkginclude_HEADERS += pencil.hpp
pkginclude_HEADERS += problem.hpp
pkginclude_HEADERS += profile.hpp
pkginclude_HEADERS += program_options.hpp
pkginclude_HEADERS += restart_definition.hpp
pkginclude_HEADERS += richardson.h
pkginclude_HEADERS += scenario_definition.hpp
pkginclude_HEADERS += state.hpp
pkginclude_HEADERS += state_impl.hpp
pkginclude_HEADERS += storage.hpp
pkginclude_HEADERS += timecontroller.hpp
pkginclude_HEADERS += timestepper.hpp
pkginclude_HEADERS += traits.hpp
pkginclude_HEADERS += types.hpp
pkginclude_HEADERS += underling_fft.h
pkginclude_HEADERS += underling_fft.hpp
pkginclude_HEADERS += underling.h
pkginclude_HEADERS += underling.hpp
pkginclude_HEADERS += utility.hpp
pkginclude_HEADERS += validation.hpp

noinst_LTLIBRARIES     += libsuzerain.la
libsuzerain_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libsuzerain_la_SOURCES  = blas_et_al.c
libsuzerain_la_SOURCES += bspline.c
libsuzerain_la_SOURCES += diffwave.c
libsuzerain_la_SOURCES += error.c
libsuzerain_la_SOURCES += fftw.cpp
libsuzerain_la_SOURCES += gbmatrix.c
libsuzerain_la_SOURCES += htstretch.c
libsuzerain_la_SOURCES += mpi.cpp
libsuzerain_la_SOURCES += os.c
libsuzerain_la_SOURCES += pencil_grid.cpp
libsuzerain_la_SOURCES += profile.cpp
libsuzerain_la_SOURCES += richardson.c

noinst_LTLIBRARIES      += libunderling.la
libunderling_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libunderling_la_SOURCES  = error.c
libunderling_la_SOURCES += underling.c
libunderling_la_SOURCES += underling_fft.c

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
