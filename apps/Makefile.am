## Process this file with automake to produce Makefile.in

AM_CPPFLAGS         = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS        += -I$(top_builddir)/lib -I$(top_srcdir)/lib -I$(builddir)
AM_CPPFLAGS        += $(UNDERLING_CFLAGS) $(P3DFFT_CFLAGS)
AM_CPPFLAGS        += $(FFTW3_THREADS_CFLAGS) $(MKL_CFLAGS) $(GSL_CFLAGS)
AM_CPPFLAGS        += $(LOG4CXX_CFLAGS) $(BOOST_CPPFLAGS)
AM_CPPFLAGS        += $(EIGEN_CFLAGS) $(GRVY_CFLAGS) $(ESIO_CFLAGS)
bin_PROGRAMS        = # Append below
noinst_PROGRAMS     = # Append below
dist_check_SCRIPTS  = # Append below
CLEANFILES          = # Append below
TESTS               = # Append below
lib_LTLIBRARIES     = # Append below
noinst_LTLIBRARIES  = # Append below
BUILT_SOURCES       = # Append below

# MPI-based compilation and linking
CXX = @MPICXX@
CC  = @MPICC@

##
## MPI-aware logging built atop log4cxx
##

noinst_LTLIBRARIES    += liblogging.la
liblogging_la_SOURCES  = logging.cpp logging.hpp
liblogging_la_LIBADD   = $(LOG4CXX_LIBS)

##
## Channel-based applications
##

lib_LTLIBRARIES       += libchannel.la
libchannel_la_SOURCES  = precision.hpp
libchannel_la_SOURCES += channel.cpp channel.hpp
libchannel_la_SOURCES += nsctpl_rholut_fwd.hpp nsctpl_rholut.hpp
libchannel_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libchannel_la_LIBADD   = $(ESIO_LIBS)
libchannel_la_LIBADD  += ../suzerain/libsuzerain.la
libchannel_la_LIBADD  += ./liblogging.la

bin_PROGRAMS         += channel_init
channel_init_LDADD    = libchannel.la
channel_init_SOURCES  = channel_init.cpp
channel_init_SOURCES += channel_init_svnrev.c
BUILT_SOURCES        += channel_init_svnrev.h
channel_init_svnrev.h: $(channel_init_SOURCES) $(libchannel_la_SOURCES)
	../build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

bin_PROGRAMS             += channel_explicit
channel_explicit_LDADD    = libchannel.la
channel_explicit_SOURCES  = channel_explicit.cpp
channel_explicit_SOURCES += nonlinear_fwd.hpp nonlinear.hpp
channel_explicit_SOURCES += explicit_op.cpp explicit_op.hpp
channel_explicit_SOURCES += hybrid_op.cpp   hybrid_op.hpp
channel_explicit_SOURCES += channel_explicit_svnrev.c
BUILT_SOURCES            += channel_explicit_svnrev.h
channel_explicit_svnrev.h: $(channel_explicit_SOURCES) $(libchannel_la_SOURCES)
	../build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

bin_PROGRAMS         += channel_mean
channel_mean_LDADD    = libchannel.la
channel_mean_SOURCES  = channel_mean.cpp
channel_mean_SOURCES += channel_mean_svnrev.c
BUILT_SOURCES        += channel_mean_svnrev.h
channel_mean_svnrev.h: $(channel_mean_SOURCES) $(libchannel_la_SOURCES)
	../build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?

##
## Miscellaneous programs and/or performance drivers
##

noinst_PROGRAMS      += conservation
conservation_SOURCES  = conservation.cpp
conservation_LDADD    = libchannel.la

if HAVE_P3DFFT
bin_PROGRAMS                += driver_suzerain_p3dfft
driver_suzerain_p3dfft_LDADD     = ../suzerain/libsuzerain.la
driver_suzerain_p3dfft_LDADD    += ./liblogging.la
driver_suzerain_p3dfft_LDADD    += $(MPIP_LIBS)
driver_suzerain_p3dfft_CPPFLAGS  = $(AM_CPPFLAGS)
driver_suzerain_p3dfft_LDFLAGS   = $(MPIP_LDFLAGS)
driver_suzerain_p3dfft_SOURCES   = driver_suzerain_p3dfft.cpp
driver_suzerain_p3dfft_SOURCES  += driver_suzerain_p3dfft_svnrev.c
BUILT_SOURCES                   += driver_suzerain_p3dfft_svnrev.h
driver_suzerain_p3dfft_svnrev.h: $(driver_suzerain_p3dfft_SOURCES) $(liblogging_la_SOURCES)
	../build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?
endif

##
## Tests for channel applications
##

# Specify test environment for 'make check'
TESTS_ENVIRONMENT  = TMPDIR=$(abs_builddir)
TESTS_ENVIRONMENT += FIELDSDIR=$(abs_top_srcdir)/fields

dist_check_SCRIPTS += test_channel_setup.sh     # Common initialization

TESTS              += test_channel_restart.sh   # Restart handling...
dist_check_SCRIPTS += test_channel_restart.sh   # ...using coefficients

TESTS              += test_channel_dparallel.sh # Serial/parallel equivalence...
dist_check_SCRIPTS += test_channel_dparallel.sh # ...for degenerate parallel

TESTS              += test_channel_parallel.sh  # Serial/parallel equivalence...
dist_check_SCRIPTS += test_channel_parallel.sh  # ...for regular parallel

TESTS              += test_channel_fields.sh    # Read canned, legacy fields
dist_check_SCRIPTS += test_channel_fields.sh

TESTS              += test_channel_fluct.sh     # Adding fluctuations
dist_check_SCRIPTS += test_channel_fluct.sh

TESTS              += test_channel_physical.sh  # Restart handling...
dist_check_SCRIPTS += test_channel_physical.sh  # ...using collocation values

# Clean up from Intel Compiler's -prof-gen output and precomiled headers
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
CLEANFILES += *.pchi

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
