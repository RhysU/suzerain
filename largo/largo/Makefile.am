## Process this file with automake to produce Makefile.in

AM_CFLAGS          = # Append below
AM_CPPFLAGS        = -I$(top_builddir)        -I$(top_srcdir)
AM_FCFLAGS         = -I$(top_builddir)        -I$(top_srcdir)
AM_LDFLAGS         = # Append below
BUILT_SOURCES      = # Append below
CLEANFILES         = # Append below
include_HEADERS    = # Append below
lib_LTLIBRARIES    = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
moduledir          = @libdir@

# Fortran compilers tend to be brittle wrt module file output
.NOTPARALLEL:

# Public largo Fortran library and module files
lib_LTLIBRARIES           += liblargof.la
liblargof_la_LDFLAGS       = -release $(PACKAGE_VERSION)
liblargof_la_LIBTOOLFLAGS  = --tag=FC
liblargof_la_SOURCES       = largo_workspace.f90
liblargof_la_SOURCES      += largo_bl_spatial.f90
liblargof_la_SOURCES      += largo_bl_spatial_first_order.f90
liblargof_la_SOURCES      += largo_bl_spatial_mixed.f90
liblargof_la_SOURCES      += largo_bl_spatiotemporal.f90
liblargof_la_SOURCES      += largo_bl_spatiotemporal_consistent.f90
liblargof_la_SOURCES      += largo_bl_temporal_chemistry.f90
liblargof_la_SOURCES      += largo_bl_temporal.f90
liblargof_la_SOURCES      += largo_bl_temporal_consistent.f90
liblargof_la_SOURCES      += largo_bl_temporal_tconsistent.f90
liblargof_la_SOURCES      += largo.f90
module_DATA                = largo_workspace.mod
module_DATA               += largo_bl_spatial.mod
module_DATA               += largo_bl_spatial_first_order.mod
module_DATA               += largo_bl_spatial_mixed.mod
module_DATA               += largo_bl_spatiotemporal.mod
module_DATA               += largo_bl_spatiotemporal_consistent.mod
module_DATA               += largo_bl_temporal_chemistry.mod
module_DATA               += largo_bl_temporal_consistent.mod
module_DATA               += largo_bl_temporal_tconsistent.mod
module_DATA               += largo_bl_temporal.mod
module_DATA               += largo.mod

# Public largo C library consumes the Fortran version
lib_LTLIBRARIES     += liblargo.la
liblargo_la_SOURCES  = largo.c
liblargo_la_SOURCES += largo.h  # Generic, unfinished
#liblargo_la_SOURCES += largo_bl_spatial.h
#liblargo_la_SOURCES += largo_bl_spatial_first_order.h
#liblargo_la_SOURCES += largo_bl_spatial_mixed.h
#liblargo_la_SOURCES += largo_bl_temporal_chemistry.h
#liblargo_la_SOURCES += largo_bl_temporal_tconsistent.h
#liblargo_la_SOURCES += largo_bl_temporal.h
liblargo_la_LDFLAGS  = -release $(PACKAGE_VERSION)
liblargo_la_LIBADD   = liblargof.la

# Public largo C headers
pkginclude_HEADERS  = largo.h version.h
#pkginclude_HEADERS += largo_bl_spatial.h
#pkginclude_HEADERS += largo_bl_spatial_first_order.h
#pkginclude_HEADERS += largo_bl_spatial_mixed.h
#pkginclude_HEADERS += largo_bl_temporal_chemistry.h
#pkginclude_HEADERS += largo_bl_temporal_tconsistent.h
#pkginclude_HEADERS += largo_bl_temporal.h

# Extra things to clean up
CLEANFILES += *.mod *.gcno *.gcna

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS =
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
