# ===========================================================================
#         http://autoconf-archive.cryp.to/ax_trilinos_amesos.html
# ===========================================================================
#
# SYNOPSIS
#
#   Test for the Trilinos Amesos
#   (http://trilinos.sandia.gov/packages/amesos) library.
#
#   AX_TRILINOS_AMESOS([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   On success, adds "include Makefile.export.amesos" statements to
#   every Automake file containing @INC_AMINCLUDE@.  Requires that
#   Trilinos was configured with the --enable-export-makefiles option.
#   When ACTION-IF-NOT-FOUND is not specified, the default behavior is
#   for configure to fail.
#
# LAST MODIFICATION
#
#   2008-08-13
#
# COPYLEFT
#
#   Copyright (c) 2008 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_TRILINOS_AMESOS],[
    AC_REQUIRE([AX_TRILINOS_BASE])
    ax_trilinos_amesos=yes
    AC_HAVE_LIBRARY([amesos],[:],[ax_trilinos_amesos=no])
    AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([amesos.macros],[ax_trilinos_amesos=no])
    AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([amesos],[ax_trilinos_amesos=no])
    if test "$ax_trilinos_amesos" = yes; then
        : # NOP
		ifelse([$1],,,
            [$1])
    else
        : # NOP
		ifelse([$2],,AC_MSG_ERROR([Trilinos Amesos not found.]),
            [$2])
    fi
])
