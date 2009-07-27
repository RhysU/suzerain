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

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@
CC=@MPICC@

pkginclude_HEADERS += blas_et_al.h
pkginclude_HEADERS += bspline_eval.h
pkginclude_HEADERS += bspline.h
pkginclude_HEADERS += error.h
pkginclude_HEADERS += exceptions.hpp
pkginclude_HEADERS += function.h
pkginclude_HEADERS += pencil_grid.hpp
pkginclude_HEADERS += pencil.hpp
pkginclude_HEADERS += underling.h

noinst_LTLIBRARIES     += libsuzerain.la
libsuzerain_la_SOURCES  = error.c bspline.c bspline_eval.c blas_et_al.c

noinst_LTLIBRARIES      += libunderling.la
libunderling_la_SOURCES  = error.c underling.c

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
