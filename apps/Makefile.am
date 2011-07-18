## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(P3DFFT_CFLAGS) $(FFTW3_THREADS_CFLAGS)
AM_CPPFLAGS   += $(MKL_CFLAGS) $(GSL_CFLAGS)
AM_CPPFLAGS   += $(LOG4CXX_CFLAGS) $(BOOST_CPPFLAGS)
AM_CPPFLAGS   += $(EIGEN_CFLAGS) $(GRVY_CFLAGS) $(ESIO_CFLAGS)
bin_PROGRAMS   = # Append below
bin_SCRIPTS    = # Append below
check_PROGRAMS = # Append below
check_SCRIPTS  = # Append below
CLEANFILES     = # Append below
TESTS          = # Append below

# MPI-based compilation and linking
CXX = @MPICXX@
CC  = @MPICC@

# Specify test environment for 'make check'
TESTS_ENVIRONMENT  = TMPDIR=$(abs_builddir)
TESTS_ENVIRONMENT += LOG4CXX_CONFIGURATION=$(abs_top_srcdir)/log4cxx.properties

# Test progams and/or performance drivers

bin_PROGRAMS                    += driver_suzerain_p3dfft
driver_suzerain_p3dfft_CPPFLAGS  = $(AM_CPPFLAGS)
driver_suzerain_p3dfft_LDFLAGS   = $(MPIP_LDFLAGS)
driver_suzerain_p3dfft_LDADD     = ../suzerain/libsuzerain.la
driver_suzerain_p3dfft_SOURCES   = driver_suzerain_p3dfft.cpp
driver_suzerain_p3dfft_SOURCES  += logger.cpp  logger.hpp
driver_suzerain_p3dfft_LDADD    += $(MPIP_LIBS) $(LOG4CXX_LIBS)

# Channel-based Applications

noinst_LTLIBRARIES     = libchannel.la
libchannel_la_SOURCES  = precision.hpp
libchannel_la_SOURCES += logger.cpp  logger.hpp
libchannel_la_SOURCES += channel.cpp channel.hpp
libchannel_la_SOURCES += nsctpl_rholut_fwd.hpp nsctpl_rholut.hpp
libchannel_la_LIBADD   = $(ESIO_LIBS) $(LOG4CXX_LIBS)
libchannel_la_LIBADD  += ../suzerain/libsuzerain.la

bin_PROGRAMS         += channel_init
channel_init_SOURCES  = channel_init.cpp
channel_init_LDADD    = libchannel.la

bin_PROGRAMS             += channel_explicit
channel_explicit_SOURCES  = channel_explicit.cpp
channel_explicit_SOURCES += explicit_op.cpp explicit_op.hpp
channel_explicit_LDADD    = libchannel.la

bin_PROGRAMS         += channel_mean
channel_mean_SOURCES  = channel_mean.cpp
channel_mean_LDADD    = libchannel.la

bin_PROGRAMS         += conservation
conservation_SOURCES  = conservation.cpp
conservation_LDADD    = libchannel.la

# Tests for channel applications

check_SCRIPTS += test_channel_setup.sh     # Common initialization

TESTS         += test_channel_restart.sh   # Restart handling
check_SCRIPTS += test_channel_restart.sh

TESTS         += test_channel_parallel.sh  # Serial/parallel equivalence tests
check_SCRIPTS += test_channel_parallel.sh

# Clean up from Intel Compiler's -prof-gen output and precomiled headers
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
CLEANFILES += *.pchi

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
