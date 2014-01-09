## Process this file with automake to produce Makefile.in

AM_CFLAGS          = # Append below
AM_CPPFLAGS        = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
AM_CPPFLAGS       += -I$(top_builddir)        -I$(top_srcdir)
AM_LDFLAGS         = # Append below
BUILT_SOURCES      = # Append below
check_PROGRAMS     = # Append below
CLEANFILES         = # Append below
include_HEADERS    = # Append below
lib_LTLIBRARIES    = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
EXTRA_DIST         = # Append below

# Use MPI compiler wrappers
CC=@MPICC@

# Make FFTW3 APIs visible during compilation
AM_CPPFLAGS += $(FFTW3_THREADS_CFLAGS)

# Make public C APIs visible in resulting libraries-- see visibility.h
AM_CPPFLAGS += -DUNDERLING_SHARED_EXPORTS

# Public underling C library
lib_LTLIBRARIES         += libunderling.la
libunderling_la_SOURCES  = error.c
libunderling_la_SOURCES += underling.c
libunderling_la_SOURCES += underling_fftw.c
libunderling_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libunderling_la_LDFLAGS += $(FFTW3_MPI_THREADS_LIBS)
libunderling_la_LDFLAGS += $(LTLIBINTL) # gnulib requirement
libunderling_la_LIBADD   = ../gnulib/libgnu.la

# Internal headers
noinst_HEADERS += common.h

# Config-time details saved using AX_PREFIX_CONFIG_H
BUILT_SOURCES += underling-config.h

# Headers shipped with the library
underling_includedir       = $(includedir)/underling
underling_include_HEADERS  = error.h
underling_include_HEADERS += periodic_function.hpp
underling_include_HEADERS += underling-config.h
underling_include_HEADERS += underling_fftw.h
underling_include_HEADERS += underling_fftw.hpp
underling_include_HEADERS += underling.h
underling_include_HEADERS += underling.hpp
underling_include_HEADERS += version.h
underling_include_HEADERS += visibility.h

# Extra things to clean up
CLEANFILES += *.gcno *.gcna

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS =
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
