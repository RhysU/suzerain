## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(GSL_CFLAGS) $(LOG4CXX_CFLAGS) $(BOOST_CPPFLAGS)
AM_CPPFLAGS   += $(P3DFFT_CFLAGS) $(FFTW3_CFLAGS)
AM_LDFLAGS     = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
AM_LDFLAGS    += $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
bin_PROGRAMS   =
check_PROGRAMS =

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@

pkginclude_HEADERS  = exceptions.hpp pencil.hpp pencil_grid.hpp
pkginclude_HEADERS += underling.h underling_error.h

noinst_LTLIBRARIES      = libunderling.la
libunderling_la_SOURCES = underling_error.c # underling.c

## check_PROGRAMS += test
## test_SOURCES    = test.c
## test_LDADD      = libunderling.la
## 
## TESTS                   = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
