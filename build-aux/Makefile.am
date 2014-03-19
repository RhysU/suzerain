## Process this file with automake to produce Makefile.in

AM_CPPFLAGS     = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS    += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_LDFLAGS      =
noinst_PROGRAMS = # Append below

AUX_DIST  = ar-lib
AUX_DIST += config.guess
AUX_DIST += config.sub
AUX_DIST += depcomp
AUX_DIST += install-sh
AUX_DIST += ltmain.sh
AUX_DIST += missing
AUX_DIST += mkinstalldirs
AUX_DIST += test-driver
AUX_DIST += tex-it

MAINTAINERCLEANFILES = $(AUX_DIST)

# Why isn't noinst_PROGRAMS automatically building svnrev on 'make dist'?
dist-hook:         svnrev

# FIXME Do some revision-specific magic accommodating git as well as SVN.
# SVNREV utility from http://www.compuphase.com/svnrev.htm
# Requires SQLite from https://www.sqlite.org/ to accommodate SVN 1.7+
# but (weirdly) svnrev.c directly #include-s sqllite3.c source code.
# Define HAVE_MREMAP=0 present to work around issue on Sysnet Desktops
noinst_PROGRAMS += svnrev
svnrev_SOURCES   = svnrev.c  svnrev.h
svnrev_CPPFLAGS  = $(AM_CPPFLAGS) -DHAVE_MREMAP=0
svnrev_LDADD     = @LIBADD_DL@

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
