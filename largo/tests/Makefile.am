## Process this file with automake to produce Makefile.in

AM_CFLAGS           = # Append below
AM_CPPFLAGS         = -I$(top_builddir)        -I$(top_srcdir)
AM_CPPFLAGS        += -I$(top_builddir)/largo
AM_FCFLAGS          = -I$(top_builddir)        -I$(top_srcdir)
AM_FCFLAGS         += -I$(top_builddir)/largo
AM_LDFLAGS          = # Append below
BUILT_SOURCES       = # Append below
check_PROGRAMS      = # Append below
CLEANFILES          = # Append below
EXTRA_DIST          = # Append below
TESTS               = # Append below
XFAIL_TESTS         = # Append below
TESTS_ENVIRONMENT   = # Append below

#############################################################################
## BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC ##
#############################################################################

## Basic C API sanity checks
TESTS              += basic
check_PROGRAMS     += basic
basic_SOURCES       = basic.c
basic_LDADD         = ../largo/liblargo.la

## Ensure Fortran testframework can return non-zero exit status
TESTS            += sanity_f
XFAIL_TESTS      += sanity_f
check_PROGRAMS   += sanity_f
sanity_f_SOURCES  = sanity_f.f90 testframework.f90
sanity_f_LDADD    = ../largo/liblargof.la

## Basic Fortran API sanity checks
TESTS           += basic_f
check_PROGRAMS  += basic_f
basic_f_SOURCES  = basic_f.f90 testframework.f90
basic_f_LDADD    = ../largo/liblargof.la

## Sensible behavior on unknown model requests
TESTS                 += unknown_model
check_PROGRAMS        += unknown_model
unknown_model_SOURCES  = unknown_model.c
unknown_model_LDADD    = ../largo/liblargo.la

## Temporal boundary layer, fortran
TESTS                 += bl_temporal_f
check_PROGRAMS        += bl_temporal_f
bl_temporal_f_SOURCES  = bl_temporal_f.f90 testframework.f90
bl_temporal_f_LDADD    = ../largo/liblargof.la

## Temporal boundary layer, baseflow, fortran
TESTS                          += bl_temporal_baseflow_f
check_PROGRAMS                 += bl_temporal_baseflow_f
bl_temporal_baseflow_f_SOURCES  = bl_temporal_baseflow_f.f90 testframework.f90
bl_temporal_baseflow_f_LDADD    = ../largo/liblargof.la

## Above case modified to provide a uniform baseflow and ensure non-NANs result
TESTS                                  += bl_temporal_baseflow_uniform_f
check_PROGRAMS                         += bl_temporal_baseflow_uniform_f
bl_temporal_baseflow_uniform_f_SOURCES  = bl_temporal_baseflow_f.f90 testframework.f90
bl_temporal_baseflow_uniform_f_FCFLAGS  = $(AM_FCFLAGS) $(FC_DEFINE)BASEFLOW_UNIFORM
bl_temporal_baseflow_uniform_f_LDADD    = ../largo/liblargof.la

## Temporal boundary layer, baseflow, fortran, rans
TESTS                               += bl_temporal_baseflow_rans_f
check_PROGRAMS                      += bl_temporal_baseflow_rans_f
bl_temporal_baseflow_rans_f_SOURCES  = bl_temporal_baseflow_rans_f.f90 testframework.f90
bl_temporal_baseflow_rans_f_LDADD    = ../largo/liblargof.la

## Temporal tensor consistent boundary layer, fortran
TESTS                             += bl_temporal_tconsistent_f
check_PROGRAMS                    += bl_temporal_tconsistent_f
bl_temporal_tconsistent_f_SOURCES  = bl_temporal_tconsistent_f.f90 testframework.f90
bl_temporal_tconsistent_f_LDADD    = ../largo/liblargof.la

## Temporal new consistent boundary layer, baseflow, fortran
TESTS                                     += bl_temporal_consistent_baseflow_f
check_PROGRAMS                            += bl_temporal_consistent_baseflow_f
bl_temporal_consistent_baseflow_f_SOURCES  = bl_temporal_consistent_baseflow_f.f90 testframework.f90
bl_temporal_consistent_baseflow_f_LDADD    = ../largo/liblargof.la

## Above case modified to provide a uniform baseflow and ensure non-NANs result
TESTS                                             += bl_temporal_consistent_baseflow_uniform_f
check_PROGRAMS                                    += bl_temporal_consistent_baseflow_uniform_f
bl_temporal_consistent_baseflow_uniform_f_SOURCES  = bl_temporal_consistent_baseflow_f.f90 testframework.f90
bl_temporal_consistent_baseflow_uniform_f_FCFLAGS  = $(AM_FCFLAGS) $(FC_DEFINE)BASEFLOW_UNIFORM
bl_temporal_consistent_baseflow_uniform_f_LDADD    = ../largo/liblargof.la

## Spatiotemporal boundary layer, baseflow, fortran
TESTS                                += bl_spatiotemporal_baseflow_f
check_PROGRAMS                       += bl_spatiotemporal_baseflow_f
bl_spatiotemporal_baseflow_f_SOURCES  = bl_spatiotemporal_baseflow_f.f90 testframework.f90
bl_spatiotemporal_baseflow_f_LDADD    = ../largo/liblargof.la

## Above case modified to provide a uniform baseflow and ensure non-NANs result
TESTS                                        += bl_spatiotemporal_baseflow_uniform_f
check_PROGRAMS                               += bl_spatiotemporal_baseflow_uniform_f
bl_spatiotemporal_baseflow_uniform_f_SOURCES  = bl_spatiotemporal_baseflow_f.f90 testframework.f90
bl_spatiotemporal_baseflow_uniform_f_FCFLAGS  = $(AM_FCFLAGS) $(FC_DEFINE)BASEFLOW_UNIFORM
bl_spatiotemporal_baseflow_uniform_f_LDADD    = ../largo/liblargof.la

## Spatiotemporal consistent boundary layer, baseflow, fortran
TESTS                                           += bl_spatiotemporal_consistent_baseflow_f
check_PROGRAMS                                  += bl_spatiotemporal_consistent_baseflow_f
bl_spatiotemporal_consistent_baseflow_f_SOURCES  = bl_spatiotemporal_consistent_baseflow_f.f90 testframework.f90
bl_spatiotemporal_consistent_baseflow_f_LDADD    = ../largo/liblargof.la

## Above case modified to provide a uniform baseflow and ensure non-NANs result
TESTS                                                   += bl_spatiotemporal_consistent_baseflow_uniform_f
check_PROGRAMS                                          += bl_spatiotemporal_consistent_baseflow_uniform_f
bl_spatiotemporal_consistent_baseflow_uniform_f_SOURCES  = bl_spatiotemporal_consistent_baseflow_f.f90 testframework.f90
bl_spatiotemporal_consistent_baseflow_uniform_f_FCFLAGS  = $(AM_FCFLAGS) $(FC_DEFINE)BASEFLOW_UNIFORM
bl_spatiotemporal_consistent_baseflow_uniform_f_LDADD    = ../largo/liblargof.la

## Spatiotemporal consistent boundary layer, baseflow, rans, fortran
TESTS                                                += bl_spatiotemporal_consistent_baseflow_rans_f
check_PROGRAMS                                       += bl_spatiotemporal_consistent_baseflow_rans_f
bl_spatiotemporal_consistent_baseflow_rans_f_SOURCES  = bl_spatiotemporal_consistent_baseflow_rans_f.f90 testframework.f90
bl_spatiotemporal_consistent_baseflow_rans_f_LDADD    = ../largo/liblargof.la

## Generic interface, temporal boundary layer, c
TESTS                         += generic_bl_temporal
check_PROGRAMS                += generic_bl_temporal
generic_bl_temporal_SOURCES    = generic_bl_temporal.c
generic_bl_temporal_LDADD      = ../largo/liblargo.la

## Generic interface, temporal boundary layer, baseflow, c
TESTS                                  += generic_bl_temporal_baseflow
check_PROGRAMS                         += generic_bl_temporal_baseflow
generic_bl_temporal_baseflow_SOURCES    = generic_bl_temporal_baseflow.c
generic_bl_temporal_baseflow_LDADD      = ../largo/liblargo.la

## Generic interface, temporal boundary layer, baseflow, rans, c
TESTS                                     += generic_bl_temporal_baseflow_rans
check_PROGRAMS                            += generic_bl_temporal_baseflow_rans
generic_bl_temporal_baseflow_rans_SOURCES  = generic_bl_temporal_baseflow_rans.c
generic_bl_temporal_baseflow_rans_LDADD    = ../largo/liblargof.la

## Generic interface, temporal tensor consistent boundary layer, c
TESTS                                   += generic_bl_temporal_tconsistent
check_PROGRAMS                          += generic_bl_temporal_tconsistent
generic_bl_temporal_tconsistent_SOURCES  = generic_bl_temporal_tconsistent.c
generic_bl_temporal_tconsistent_LDADD    = ../largo/liblargo.la

## Generic interface, temporal consistent boundary layer, baseflow, c
TESTS                                           += generic_bl_temporal_consistent_baseflow
check_PROGRAMS                                  += generic_bl_temporal_consistent_baseflow
generic_bl_temporal_consistent_baseflow_SOURCES  = generic_bl_temporal_consistent_baseflow.c
generic_bl_temporal_consistent_baseflow_LDADD    = ../largo/liblargo.la

## Generic interface, spatiotemporal boundary layer, baseflow, c
TESTS                                        += generic_bl_spatiotemporal_baseflow
check_PROGRAMS                               += generic_bl_spatiotemporal_baseflow
generic_bl_spatiotemporal_baseflow_SOURCES    = generic_bl_spatiotemporal_baseflow.c
generic_bl_spatiotemporal_baseflow_LDADD      = ../largo/liblargo.la

## Generic interface, spatiotemporal consistent boundary layer, baseflow, c
TESTS                                                 += generic_bl_spatiotemporal_consistent_baseflow
check_PROGRAMS                                        += generic_bl_spatiotemporal_consistent_baseflow
generic_bl_spatiotemporal_consistent_baseflow_SOURCES  = generic_bl_spatiotemporal_consistent_baseflow.c
generic_bl_spatiotemporal_consistent_baseflow_LDADD    = ../largo/liblargo.la

## Generic interface, temporal boundary layer, fortran
TESTS                         += generic_bl_temporal_f
check_PROGRAMS                += generic_bl_temporal_f
generic_bl_temporal_f_SOURCES  = generic_bl_temporal_f.f90 testframework.f90
generic_bl_temporal_f_LDADD    = ../largo/liblargof.la

## Generic interface, temporal boundary layer, baseflow, fortran
TESTS                                  += generic_bl_temporal_baseflow_f
check_PROGRAMS                         += generic_bl_temporal_baseflow_f
generic_bl_temporal_baseflow_f_SOURCES  = generic_bl_temporal_baseflow_f.f90 testframework.f90
generic_bl_temporal_baseflow_f_LDADD    = ../largo/liblargof.la

## Generic interface, temporal tensor consistent boundary layer, fortran
TESTS                                     += generic_bl_temporal_tconsistent_f
check_PROGRAMS                            += generic_bl_temporal_tconsistent_f
generic_bl_temporal_tconsistent_f_SOURCES  = generic_bl_temporal_tconsistent_f.f90 testframework.f90
generic_bl_temporal_tconsistent_f_LDADD    = ../largo/liblargof.la

## Generic interface, temporal consistent boundary layer, fortran
TESTS                                              += generic_bl_temporal_consistent_baseflow_f
check_PROGRAMS                                     += generic_bl_temporal_consistent_baseflow_f
generic_bl_temporal_consistent_baseflow_f_SOURCES   = generic_bl_temporal_consistent_baseflow_f.f90 testframework.f90
generic_bl_temporal_consistent_baseflow_f_LDADD     = ../largo/liblargof.la

## Generic interface, spatiotemporal boundary layer, baseflow, fortran
TESTS                                        += generic_bl_spatiotemporal_baseflow_f
check_PROGRAMS                               += generic_bl_spatiotemporal_baseflow_f
generic_bl_spatiotemporal_baseflow_f_SOURCES  = generic_bl_spatiotemporal_baseflow_f.f90 testframework.f90
generic_bl_spatiotemporal_baseflow_f_LDADD    = ../largo/liblargof.la

## Generic interface, spatiotemporal consistent boundary layer, baseflow, fortran
TESTS                                                   += generic_bl_spatiotemporal_consistent_baseflow_f
check_PROGRAMS                                          += generic_bl_spatiotemporal_consistent_baseflow_f
generic_bl_spatiotemporal_consistent_baseflow_f_SOURCES  = generic_bl_spatiotemporal_consistent_baseflow_f.f90 testframework.f90
generic_bl_spatiotemporal_consistent_baseflow_f_LDADD    = ../largo/liblargof.la

#############################################################################
### MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS ###
#############################################################################

# Build test programs without running 'make check'
.PHONY: check-programs
check-programs: $(check_PROGRAMS)

# Extra things to clean up
clean-local:
	rm -f *.mod *.gcno *.gcda

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS = -i largo
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
