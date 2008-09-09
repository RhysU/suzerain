# ===========================================================================
#        http://autoconf-archive.cryp.to/ax_prog_ctangle.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PROG_CTANGLE([ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Search for Levy/Knuth's CWEB ctangle implementation.  Set the
#   variable CTANGLE to the name of the application if found.
#   If ACTION-IF-NOT-FOUND is not specified and ctangle is not located,
#   configure fails.
#
# LAST MODIFICATION
#
#   2008-08-11
#
# COPYLEFT
#
#   Copyright (c) 2008 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   This library is free software; you can redistribute it and/or modify it
#   under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation; either version 2.1 of the License, or (at
#   your option) any later version.
#
#   This library is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
#   General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with this library. If not, see <http://www.gnu.org/licenses/>.

AC_DEFUN([AX_PROG_CTANGLE],[
AC_ARG_VAR(CTANGLE,[Command to run Levy/Knuth CWEB's ctangle])

if test -n "$CTANGLE"; then
    AC_CHECK_PROG([ax_prog_ctangle],[$CTANGLE],[$CTANGLE],[no])
else
    AC_CHECK_PROGS([ax_prog_ctangle],[ctangle],[no])
fi

if test "$ax_prog_ctangle" = "no"; then
	ifelse($#,0,[AC_MSG_ERROR([Levy/Knuth CWEB's ctangle not found.])],
        $1)
else
    if "$ax_prog_ctangle" --version | grep "Levy" &>/dev/null; then
        ax_prog_ctangle_type="Levy/Knuth"
    elif "$ax_prog_ctangle" | grep "friend" &>/dev/null; then
        AC_MSG_WARN([Found ctangle from van Leeuwen's CWEBx-- using compatibility flag, but typesetting may give errors.])
        ax_prog_ctangle="$ax_prog_ctangle +c"
        ax_prog_ctangle_type="van Leeuwen"
    else
        AC_MSG_WARN([Found unknown ctangle implementation-- proceeding, but typesetting may give errors])
        ax_prog_ctangle_type="unknown"
    fi

    AC_SUBST(CTANGLE,$ax_prog_ctangle)
fi
])
