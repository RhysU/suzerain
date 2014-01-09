## Process this file with automake to produce Makefile.in

# Note this Makefile.am lacks Gnulib-related -I options like
#    AM_CPPFLAGS     = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
#    AM_FCFLAGS      = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
# This example code deliberately lacks #include <config.h> dependence
# and therefore cannot rely on Gnulib per
# http://www.gnu.org/software/gnulib/manual/gnulib.html#Source-changes

AM_CFLAGS       =
AM_CPPFLAGS     = -I$(top_builddir) -I$(top_srcdir)
AM_LDFLAGS      = # Append below
noinst_PROGRAMS = # Append below
BUILT_DOCS      = # Append below
BUILT_SOURCES   = # Append below
CLEANFILES      = # Append below

# Use MPI compiler wrappers since HDF5 wrappers break autodependencies.
CC=@MPICC@

## Sample progams
#noinst_PROGRAMS        += concepts1
#concepts1_SOURCES       = concepts1.c
#concepts1_CFLAGS        = ../underling/libunderling.la $(FFTW3_MPI_THREADS_CFLAGS)
#concepts1_LDADD         = ../underling/libunderling.la $(FFTW3_MPI_THREADS_LIBS)

# Additional files to clean
CLEANFILES += *.gcno *.gcda

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS = -i underling
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
