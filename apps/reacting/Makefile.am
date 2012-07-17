## Process this file with automake to produce Makefile.in

include ../Common.am

# Build targets in current directory before any subdirectories
if HAVE_CANTERA
    MAYBE_CANTERA = cantera
endif
SUBDIRS = . $(MAYBE_CANTERA)

##
## Logic common to potentially many applications kept in subdirectories
##

##
## Miscellaneous programs and/or performance drivers
##

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
