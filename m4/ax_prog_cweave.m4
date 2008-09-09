# ===========================================================================
#        http://autoconf-archive.cryp.to/ax_prog_cweave.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PROG_CWEAVE([ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Search for Levy/Knuth's CWEB cweave implementation.  Set the
#   variable CWEAVE to the name of the application if found.
#   If ACTION-IF-NOT-FOUND is not specified and cweave is not located,
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

AC_DEFUN([AX_PROG_CWEAVE],[
AC_ARG_VAR(CWEAVE,[Command to run Levy/Knuth CWEB's cweave])

if test -n "$CWEAVE"; then
    AC_CHECK_PROG([ax_prog_cweave],[$CWEAVE],[$CWEAVE],[no])
else
    AC_CHECK_PROGS([ax_prog_cweave],[cweave],[no])
fi

if test "$ax_prog_cweave" = "no"; then
	ifelse($#,0,[AC_MSG_ERROR([Levy/Knuth CWEB's cweave not found.])],
        $1)
else
    if "$ax_prog_cweave" --version | grep "Levy" &>/dev/null; then
        ax_prog_cweave_type="Levy/Knuth"
    elif "$ax_prog_cweave" | grep "friend" &>/dev/null; then
        AC_MSG_WARN([Found cweave from van Leeuwen's CWEBx-- using compatibility flag, but compilation may give errors.])
        ax_prog_cweave="$ax_prog_cweave +c"
        ax_prog_cweave_type="van Leeuwen"
    else
        AC_MSG_WARN([Found unknown cweave implementation-- proceeding, but compilation may give errors.])
        ax_prog_cweave_type="unknown"
    fi

    AC_SUBST(CWEAVE,$ax_prog_cweave)
fi
])
