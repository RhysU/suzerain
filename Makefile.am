## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS  = -I m4
AUTOMAKE_OPTIONS = foreign
CLEANFILES       = # Append below

# Build in these directories:
if BUILD_WRITEUP
    MAYBE_WRITEUP_SUBDIR = writeup
endif
SUBDIRS = lib suzerain tests apps $(MAYBE_WRITEUP_SUBDIR)

# Distribute these items:
DIST_SUBDIRS = lib suzerain apps tests writeup
EXTRA_DIST = bootstrap $(DX_CONFIG)

# Build toplevel directory in serial
# http://www.gnu.org/software/automake/manual/make/Special-Targets.html
.NOTPARALLEL:

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

# Target for all documentation
docs: doxygen-run

# Make sure that libtool is up to date
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Required for AX_PROG_DOXYGEN
include doxygen.am

# -------------------------------------------
# Optional support for code coverage analysis
# -------------------------------------------
# DISABLED: Pattern to remove some lines from the report
# $(top_srcdir)/lcov/lcov --list-full-path -l $(lcov_dir)/lcov.info | grep "grvy_hpct_compat.cpp" | cut -d\| -f1 >> $(lcov_dir)/remove
if CODE_COVERAGE_ENABLED
lcov_dir=$(top_builddir)/docs/html/lcov

lcov-report:
	@mkdir -p $(lcov_dir)
	$(top_srcdir)/lcov/lcov --compat-libtool --directory . --capture --output-file $(lcov_dir)/lcov.info
	$(top_srcdir)/lcov/lcov --list-full-path -l $(lcov_dir)/lcov.info | grep -v "`cd -P $(top_srcdir)/suzerain && pwd`" | cut -d\| -f1 > $(lcov_dir)/remove
	$(top_srcdir)/lcov/lcov -q -r $(lcov_dir)/lcov.info `cat $(lcov_dir)/remove` > $(lcov_dir)/lcov.cleaned.info
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

