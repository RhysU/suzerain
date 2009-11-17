## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS  = -I m4
AUTOMAKE_OPTIONS = foreign
CLEANFILES       = # Append below

# Build in these directories:
SUBDIRS = lib suzerain tests apps doc

# Distribute these items:
DIST_SUBDIRS = lib suzerain tests apps doc
EXTRA_DIST= bootstrap $(DX_CONFIG)

# Eliminate .svn directories in dist tarball
dist-hook:
	find $(distdir)/ -name '.svn' -type d -exec rm {} \;

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh    \
           build-aux/mkinstalldirs \
           build-aux/missing

# Clean up from Intel Compiler's -prof-gen output
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
clean-local:
	rm -rf lib/pgopti.spi lib/pgopti.spl lib/*.dyn lib/*.dpi

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
	suzerain/config.h.in   \
	suzerain/Makefile.in   \
	suzerain/stamp.h.in

# Additional tasks run with the 'all' target
all-local:  libtool

# Documentation hooks in with the 'html' target
html-local: doxygen-run

# Make sure that libtool is up to date
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Required for AX_PROG_DOXYGEN
include doxygen.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
