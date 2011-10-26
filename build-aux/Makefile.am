## Process this file with automake to produce Makefile.in

AM_CPPFLAGS     = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS    += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_LDFLAGS      =
noinst_PROGRAMS = # Append below

AUX_DIST = install-sh    \
           mkinstalldirs \
           missing

MAINTAINERCLEANFILES = $(AUX_DIST)

noinst_PROGRAMS += svnrev
svnrev_SOURCES   = svnrev.c svnrev.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
