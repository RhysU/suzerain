## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS  = -I m4
AUTOMAKE_OPTIONS = foreign
BUILT_SOURCES    = # Append below
CLEANFILES       = # Append below
EXTRA_DIST       = # Append below

# Build toplevel directory in serial
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

# Build in these directories:
if BUILD_WRITEUPS
    MAYBE_WRITEUPS_SUBDIR = writeups
endif
if BUILD_TESTS
    MAYBE_TESTS_SUBDIR = tests
endif
SUBDIRS = lib largo build-aux suzerain $(MAYBE_TESTS_SUBDIR) fields apps \
          $(MAYBE_WRITEUPS_SUBDIR) jobscripts preproc postproc

# Distribute these directories:
DIST_SUBDIRS = lib largo build-aux suzerain tests fields apps \
               writeups jobscripts preproc postproc
EXTRA_DIST  += AUTHORS astylerc bootstrap COPYING $(DX_CONFIG)
EXTRA_DIST  += notebooks # Entire directory contents

# Eliminate .svn and/or .git directories in dist tarball
dist-hook:
	find $(distdir)/ -name '.svn' -type d -exec rm {} \;
	find $(distdir)/ -name '.git' -type d -exec rm {} \;

# Clean up from Intel Compiler's -prof-gen output
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
clean-local:
	rm -rf lib/pgopti.spi lib/pgopti.spl lib/*.dyn lib/*.dpi

MOSTLYCLEANFILES = $(DX_CLEANFILES)

DISTCLEANFILES = _configs.sed               \
                 aminclude.am               \
                 suzerain/suzerain-config.h

MAINTAINERCLEANFILES = aclocal.m4             \
                       autom4te.cache/*       \
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

# Target for all documentation
docs: doxygen-run

# Make sure that libtool is up to date
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Required for AX_PROG_DOXYGEN
include doxygen.am

##
## Version number management using build-aux/version-gen modifying ideas from
## https://lists.gnu.org/archive/html/autoconf/2009-12/msg00056.html
## and those embedded within build-aux/version-gen itself.  That mailing
## list suggests more robust tarball procedures worth eventually adopting.
##
EXTRA_DIST    += .version
BUILT_SOURCES += $(srcdir)/.version
$(srcdir)/.version: configure.ac
	echo $(VERSION) >$@-t && mv $@-t $@
dist-hook: dist-hook-version-dist
dist-hook-version-dist:
	echo $(VERSION) >$(distdir)/.version-dist

# -------------------------------------------
# Optional support for code coverage analysis
# -------------------------------------------
if CODE_COVERAGE_ENABLED
lcov_dir=$(top_builddir)/docs/html/lcov

lcov-report:
	@mkdir -p $(lcov_dir)
	$(top_srcdir)/lcov/lcov --compat-libtool --directory . --capture --output-file $(lcov_dir)/lcov.info
	$(top_srcdir)/lcov/lcov --list-full-path -l $(lcov_dir)/lcov.info \
	    | egrep -v "`cd -P $(top_srcdir) && pwd`/(suzerain|apps)" \
	    | cut -d\| -f1 > $(lcov_dir)/remove
	$(top_srcdir)/lcov/lcov -r $(lcov_dir)/lcov.info `cat $(lcov_dir)/remove` -o $(lcov_dir)/lcov.cleaned.info
	@rm $(lcov_dir)/remove
	@mv $(lcov_dir)/lcov.cleaned.info $(lcov_dir)/lcov.info
	$(top_srcdir)/lcov/genhtml -t "suzerain" -o $(lcov_dir) $(lcov_dir)/lcov.info

lcov-reset:
	@rm -rf $(lcov_dir)
	@find . -name "*.gcda" -exec rm {} \;
	$(top_srcdir)/lcov/lcov --directory . --zerocounters

coverage: lcov-reset check lcov-report
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
