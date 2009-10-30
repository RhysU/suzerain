# SYNOPSIS
#
#   AX_FORCEINLINE
#
# DESCRIPTION
#
#   Provides a test for the compiler support of forced inlining.  If usable,
#   #define FORCEINLINE to the appropriate force inline keyword.  Otherwise
#   #define FORCEINLINE to be 'inline'.
#
# LICENSE
#
#   Copyright (c) 2008 Alan Woodland <ajw05@aber.ac.uk>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_FORCEINLINE], [
  AC_LANG_PUSH([C])
  AC_MSG_CHECKING(for forced inline keyword)
  AC_CACHE_VAL(ac_cv_forceinline, [
    ax_forceinline_keywords="__forceinline inline none"
    for ax_forceinline_keyword in $ax_forceinline_keywords; do
       case $ax_forceinline_keyword in
          none) ac_cv_forceinline=none ; break ;;
      *)
             AC_TRY_COMPILE(
                [#include <stdlib.h>
                 ] $ax_forceinline_keyword [
                 static void
                 foo(void) {
                 exit(1);
                 }],
                 [],
                 [ac_cv_forceinline=$ax_forceinline_keyword ; break],
                 ac_cv_forceinline=none
             )
      esac
    done
])

  if test "$ac_cv_forceinline" = "none"; then
    ax_forceinline_keyword=
  else
    ax_forceinline_keyword=$ac_cv_forceinline
  fi
  AC_DEFINE_UNQUOTED([FORCEINLINE],$ax_forceinline_keyword,
    [The most forceful inline keyword known by the compiler])
  AC_MSG_RESULT($ac_cv_forceinline)
  AC_LANG_POP([C])
])
