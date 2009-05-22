## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign

#Build in these directories:
SUBDIRS = lib doc suzerain tests

#Distribute these directories:
DIST_SUBDIRS = lib doc suzerain tests

EXTRA_DIST= bootstrap

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh    \
           build-aux/mkinstalldirs \
           build-aux/missing

MOSTLYCLEANFILES = $(DX_CLEANFILES)

# Files to be deleted by 'make maintainer-clean'
MAINTAINERCLEANFILES =     \
	aclocal.m4             \
	aminclude.am           \
	autom4te.cache/*       \
	$(AUX_DIST)            \
	config.log             \
	config.status          \
	config.sub             \
	configure              \
	Makefile.in            \
	src/config.h.in        \
	src/Makefile.in        \
	src/stamp.h.in

# Additional tasks run with the 'all' target
all-local: libtool

# Make sure that libtool is up to date
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Required for AX_PROG_DOXYGEN
include doxygen.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
