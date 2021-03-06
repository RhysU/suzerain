## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign

# Build in these directories:
if BUILD_WRITEUPS
    MAYBE_WRITEUPS_SUBDIR = writeups
endif
SUBDIRS = largo tests $(MAYBE_WRITEUPS_SUBDIR) doxygen

# Distribute these items:
DIST_SUBDIRS = largo tests writeups doxygen
EXTRA_DIST   = # bootstrap not to be distributed in tarballs
EXTRA_DIST  += $(DX_CONFIG)
EXTRA_DIST  += build-aux/sfmakedepend.pl
EXTRA_DIST  += build-aux/update_license.pl

dist-hook: dist-hook-scm dist-hook-docs
dist-hook-scm:
	: # Eliminate .svn and .git directories in dist tarball
	find $(distdir)/ -name '.svn' -type d -exec rm {} \;
	find $(distdir)/ -name '.git' -type d -exec rm {} \;
dist-hook-docs:
	: # Include generated documentation in the dist tarball
	cp -rv $(top_builddir)/docs $(top_distdir)

# Tools in the auxillary directory
AUX_DIST  = build-aux/install-sh
AUX_DIST += build-aux/mkinstalldirs
AUX_DIST += build-aux/missing

# Support for pkgconfig
pkgconfigdir   = $(libdir)/pkgconfig
pkgconfig_DATA = largo.pc

# Files to be deleted by make clean variants
MOSTLYCLEANFILES =  $(DX_CLEANFILES)
DISTCLEANFILES =           \
    aminclude.am           \
    _configs.sed           \
    *__genmod.f90          \
    *__genmod.mod
MAINTAINERCLEANFILES =     \
    aclocal.m4             \
    autom4te.cache/*       \
    $(AUX_DIST)            \
    config.log             \
    config.status          \
    config.sub             \
    configure              \
    Makefile.in            \
    largo/config.h.in      \
    largo/Makefile.in      \
    largo/stamp.h.in

# Run toplevel makefile serially; submakes may still be parallel
.NOTPARALLEL:

# Additional tasks run with the 'all' target
all-local:  libtool

# Documentation hooks in with the 'html' target
html-local: doxygen-run

# Target for all documentation
docs: doxygen-run

# Target to build test-related programs within tests subdirectory
check-programs: all
	cd tests && $(MAKE) $(AM_MAKEFLAGS) check-programs

# -------------------------------------------
# Optional support for code coverage analysis
# -------------------------------------------

if CODE_COVERAGE_ENABLED

lcov_dir=$(top_builddir)/docs/html/lcov

# General philosophy is to maintain code coverage for the
# base library as generated by "make check" tests.

lcov-report:
	@mkdir -p $(lcov_dir)

	$(top_srcdir)/lcov/lcov --compat-libtool --directory . --capture --output-file $(lcov_dir)/lcov.info
	$(top_srcdir)/lcov/lcov --list-full-path -l $(lcov_dir)/lcov.info | grep -v "`cd -P $(top_srcdir)/largo && pwd`" | cut -d\| -f1 > $(lcov_dir)/remove
	$(top_srcdir)/lcov/lcov -q -r $(lcov_dir)/lcov.info `cat $(lcov_dir)/remove` > $(lcov_dir)/lcov.cleaned.info
	@rm $(lcov_dir)/remove
	@mv $(lcov_dir)/lcov.cleaned.info $(lcov_dir)/lcov.info
	$(top_srcdir)/lcov/genhtml -t "largo" -o $(lcov_dir) $(lcov_dir)/lcov.info

lcov-reset:
	@rm -rf $(lcov_dir)
	@find . -name "*.gcda" -exec rm {} \;
	$(top_srcdir)/lcov/lcov --directory . --zerocounters

coverage: docs lcov-reset check lcov-report

endif

# Make sure that libtool is up to date
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Push license updates into appropriate subdirectories
.license.stamp:
	cd largo && $(MAKE) $(AM_MAKEFLAGS) .license.stamp
	cd tests && $(MAKE) $(AM_MAKEFLAGS) .license.stamp

# Required for AX_PROG_DOXYGEN
include doxygen/doxygen.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
