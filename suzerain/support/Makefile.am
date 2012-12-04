## Process this file with automake to produce Makefile.in

AM_CPPFLAGS         = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS        += -I$(top_builddir)/lib -I$(top_srcdir)/lib -I$(builddir)
AM_CPPFLAGS        += $(UNDERLING_CFLAGS) $(P3DFFT_CFLAGS)
AM_CPPFLAGS        += $(FFTW3_THREADS_CFLAGS) $(MKL_CFLAGS) $(GSL_CFLAGS)
AM_CPPFLAGS        += $(LOG4CXX_CFLAGS) $(BOOST_CPPFLAGS)
AM_CPPFLAGS        += $(EIGEN_CFLAGS) $(GRVY_CFLAGS) $(ESIO_CFLAGS)
bin_PROGRAMS       = # Append below
check_PROGRAMS     = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
lib_LTLIBRARIES    = # Append below
pkginclude_HEADERS = # Append below
CLEANFILES         = # Append below

# MPI-based compilation and linking
CXX = @MPICXX@
CC  = @MPICC@

# Implementation-only headers
# noinst_HEADERS +=

# Public API headers
supportdir = $(pkgincludedir)/support
# support_HEADERS +=

lib_LTLIBRARIES                      += libsuzerain-support.la
libsuzerain_support_la_SOURCES        =
EXTRA_libsuzerain_support_la_SOURCES  =
libsuzerain_support_la_LDFLAGS        = -release $(PACKAGE_VERSION)
libsuzerain_support_la_LIBADD         = ../libsuzerain.la
libsuzerain_support_la_LIBADD        += $(ESIO_LIBS)

# Clean up from Intel Compiler's -prof-gen output and precompiled headers
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
CLEANFILES += *.pchi

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
