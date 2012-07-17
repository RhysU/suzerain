## Process this file with automake to produce Makefile.in

include Common.am

# Build targets in current directory before any subdirectories
SUBDIRS = . perfect reacting

##
## MPI-aware logging built atop log4cxx
##

noinst_LTLIBRARIES    += liblogging.la
liblogging_la_SOURCES  = logging.cpp logging.hpp
liblogging_la_LIBADD   = $(LOG4CXX_LIBS)

##
## Support logic common to potentially many applications kept in subdirectories
##

noinst_LTLIBRARIES    += libsupport.la
libsupport_la_SOURCES  = precision.hpp
libsupport_la_SOURCES += support.cpp support.hpp
libsupport_la_SOURCES += nsctpl_rholut_fwd.hpp nsctpl_rholut.hpp
libsupport_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libsupport_la_LIBADD   = $(ESIO_LIBS)
libsupport_la_LIBADD  += ../suzerain/libsuzerain.la
libsupport_la_LIBADD  += ./liblogging.la

##
## Miscellaneous programs and/or performance drivers
##

noinst_PROGRAMS      += conservation
conservation_SOURCES  = conservation.cpp
conservation_LDADD    = libsupport.la

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
	$(top_builddir)/build-aux/svnrev -n -o$@ "-f$(PACKAGE_VERSION).#" -i $?
endif

##
## Common test infrastructure and tests
##

dist_check_SCRIPTS += test_infrastructure.sh

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
