# ===========================================================================
#         http://autoconf-archive.cryp.to/ax_trilinos_epetraext.html
# ===========================================================================
#
# SYNOPSIS
#
#   Test for the Trilinos EpetraExt
#   (http://trilinos.sandia.gov/packages/epetraext) library.
#
#   AX_TRILINOS_EPETRAEXT([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   On success, adds "include Makefile.export.epetraext" statements to
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

AC_DEFUN([AX_TRILINOS_EPETRAEXT],[
    AC_REQUIRE([AX_TRILINOS_BASE])
    ax_trilinos_epetraext=yes
    AC_HAVE_LIBRARY([epetraext],[:],[ax_trilinos_epetraext=no])
    AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([epetraext.macros],[ax_trilinos_epetraext=no])
    AX_ADD_AM_TRILINOS_MAKEFILE_EXPORT([epetraext],[ax_trilinos_epetraext=no])
    if test "$ax_trilinos_epetraext" = yes; then
        : # NOP
		ifelse([$1],,,
            [$1])
    else
        : # NOP
		ifelse([$2],,AC_MSG_ERROR([Trilinos EpetraExt not usable.]),
            [$2])
    fi
])
