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

bin_PROGRAMS    += channel
channel_LDADD    = libchannel.la
channel_SOURCES  = channel_main.cpp
channel_SOURCES += channel_treatment.hpp
channel_SOURCES += nonlinear_fwd.hpp nonlinear.hpp
channel_SOURCES += explicit_op.cpp explicit_op.hpp
channel_SOURCES += hybrid_op.cpp   hybrid_op.hpp
channel_SOURCES += timers.hpp
channel_SOURCES += channel_main_svnrev.c
BUILT_SOURCES   += channel_main_svnrev.h
channel_main_svnrev.h: $(channel_SOURCES) $(libchannel_la_SOURCES)
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
bin_PROGRAMS                    += driver_suzerain_p3dfft
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

# Clean up from Intel Compiler's -prof-gen output and precomiled headers
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
CLEANFILES += *.pchi

##
## Tests for channel applications
##
## Identical test logic varying only in the operator choices (e.g. --explicit
## versus implicit) are broken out into separate, templated tests to increase
## the parallelism available during 'make check'.
##

# Specify test environment for 'make check'
TESTS_ENVIRONMENT  = TMPDIR=$(abs_builddir)
TESTS_ENVIRONMENT += FIELDSDIR=$(abs_top_srcdir)/fields

dist_check_SCRIPTS += test_channel_setup.sh     # Common initialization

# Restart handling using coefficient-based restart files
TESTS              += test_channel_restart_explicit.sh
TESTS              += test_channel_restart_implicit.sh
dist_check_SCRIPTS += test_channel_restart_template.sh
dist_check_SCRIPTS += test_channel_restart_explicit.sh
dist_check_SCRIPTS += test_channel_restart_implicit.sh

# Serial/parallel equivalence for degenerate parallel decompositions
# Not important from a utility perspective, but great for hitting edge cases
TESTS              += test_channel_dparallel.sh
dist_check_SCRIPTS += test_channel_dparallel.sh

# Serial/parallel equivalence for normal parallel decompositions
TESTS              += test_channel_parallel_explicit.sh
TESTS              += test_channel_parallel_implicit.sh
dist_check_SCRIPTS += test_channel_parallel_template.sh
dist_check_SCRIPTS += test_channel_parallel_explicit.sh
dist_check_SCRIPTS += test_channel_parallel_implicit.sh

# Read canned, legacy fields to test backwards compatibility
TESTS              += test_channel_fields.sh
dist_check_SCRIPTS += test_channel_fields.sh

# Adding fluctuations at startup
TESTS              += test_channel_fluct.sh
dist_check_SCRIPTS += test_channel_fluct.sh

# Restart handling using collocation-value-based restart files
TESTS              += test_channel_physical_explicit.sh
TESTS              += test_channel_physical_implicit.sh
dist_check_SCRIPTS += test_channel_physical_template.sh
dist_check_SCRIPTS += test_channel_physical_explicit.sh
dist_check_SCRIPTS += test_channel_physical_implicit.sh

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
