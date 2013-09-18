# ===========================================================================
#            http://autoconf-archive.cryp.to/acltx_prog_tex.html
# ===========================================================================
#
# SYNOPSIS
#
#   ACLTX_PROG_TEX([ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   If TEX=PROGRAM is specified, do nothing.  Otherwise find TeX and
#   set TEX to the name of the application.  If ACTION-IF-NOT-FOUND is
#   not specified, configure fails when TeX is not found.
#
# LAST MODIFICATION
#
#   2008-08-06
#
# COPYLEFT
#
#   Copyright (c) 2008 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Boretti Mathieu <boretti@eig.unige.ch>
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

AC_DEFUN([ACLTX_PROG_TEX],[
AC_ARG_VAR(TEX,[specify default TeX application])
if test "x$TEX" != "x" ; then
    AC_MSG_CHECKING([Checking for tex])
    AC_MSG_RESULT([$TEX (from parameter)])
else
    AC_CHECK_PROGS(TEX,[tex])
fi
if test "x$TEX" = "x" ;
then
	ifelse($#,0,[AC_MSG_ERROR([Unable to find TeX])],
        $1)
fi])
