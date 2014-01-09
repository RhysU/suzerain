## Process this file with automake to produce Makefile.in

AM_CFLAGS           = $(GRVY_CFLAGS) $(FFTW3_MPI_THREADS_CFLAGS)
AM_CPPFLAGS         = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
AM_CPPFLAGS        += -I$(top_builddir)        -I$(top_srcdir)
AM_LDFLAGS          = # Append below
bin_PROGRAMS        = # Append below
BUILT_DOCS          = # Append below
BUILT_SOURCES       = # Append below
CLEANFILES          = # Append below
dist_check_SCRIPTS  = # Append below
noinst_HEADERS      = mpi_argp.h
noinst_PROGRAMS     = # Append below
TESTS               = # Append below
TESTS_ENVIRONMENT   = # Append below

# Use MPI compiler wrappers
CC=@MPICC@

## C-based benchmarking tool
bin_PROGRAMS             += underling_bench
underling_bench_SOURCES   = underling_bench.c mpi_argp.c
underling_bench_LDADD     = ../underling/libunderling.la
underling_bench_LDADD    += ../gnulib/libgnu.la $(LTLIBINTL)
underling_bench_LDADD    += $(GRVY_LIBS)

## Tests built atop the benchmarking tool
dist_check_SCRIPTS += test_setup.sh   # Common initialization
TESTS              += test_bench.sh   # A varienty of
dist_check_SCRIPTS += test_bench.sh   # ...benchmark cases

# Additional files to clean
CLEANFILES += *.gcno *.gcda

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS = -i underling
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
