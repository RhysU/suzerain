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
noinst_PROGRAMS += svnrev
svnrev_SOURCES   = svnrev.c svnrev.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
