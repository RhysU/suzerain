## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(MKL_CFLAGS) $(GSL_CFLAGS) $(LOG4CXX_CFLAGS)
AM_CPPFLAGS   += $(BOOST_CPPFLAGS) $(P3DFFT_CFLAGS) $(FFTW3_CFLAGS)
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

pkginclude_HEADERS += blas_et_al.h
pkginclude_HEADERS += blas_et_al.hpp
pkginclude_HEADERS += bspline.h
pkginclude_HEADERS += common.h
pkginclude_HEADERS += complex.hpp
pkginclude_HEADERS += error.h
pkginclude_HEADERS += function.h
pkginclude_HEADERS += functional.hpp
pkginclude_HEADERS += math.hpp
pkginclude_HEADERS += multi_array.hpp
pkginclude_HEADERS += orthonormal.hpp
pkginclude_HEADERS += pencilfft.hpp
pkginclude_HEADERS += pencil_grid.hpp
pkginclude_HEADERS += pencil.hpp
pkginclude_HEADERS += problem.hpp
pkginclude_HEADERS += program_options.hpp
pkginclude_HEADERS += richardson.h
pkginclude_HEADERS += suzerain-config.h
pkginclude_HEADERS += timestepper.hpp
pkginclude_HEADERS += underling.h

noinst_LTLIBRARIES     += libsuzerain.la
libsuzerain_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libsuzerain_la_SOURCES  = blas_et_al.c
libsuzerain_la_SOURCES += bspline.c
libsuzerain_la_SOURCES += error.c
libsuzerain_la_SOURCES += richardson.c

noinst_LTLIBRARIES      += libunderling.la
libunderling_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libunderling_la_SOURCES  = error.c
libunderling_la_SOURCES += underling.c

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
