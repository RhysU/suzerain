# SYNOPSIS
#
#   AX_ASM_INLINE
#
# DESCRIPTION
#
#   Provides a test for the compiler support of inline assembly instructions.
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

AC_DEFUN([AX_ASM_INLINE], [
  AC_LANG_PUSH([C])
  AC_MSG_CHECKING(for inline assembly style)
  AC_CACHE_VAL(ac_cv_asm_inline, [
    ax_asm_inline_keywords="__asm__ __asm none"
    for ax_asm_inline_keyword in $ax_asm_inline_keywords; do
       case $ax_asm_inline_keyword in
          none) ac_cv_asm_inline=none ; break ;;
      *)
             AC_TRY_COMPILE(
                [#include <stdlib.h>
                 static void
                 foo(void) {
                 ] $ax_asm_inline_keyword [("");
                 exit(1);
                 }],
                 [],
                 [ac_cv_asm_inline=$ax_asm_inline_keyword ; break],
                 ac_cv_asm_inline=none
             )
      esac
    done
])

  if test "$ac_cv_asm_inline" != "none"; then
    AC_DEFINE_UNQUOTED([ASM_INLINE], $ac_cv_asm_inline, [If the compiler supports inline assembly define it to that keyword here])
  fi
  AC_MSG_RESULT($ac_cv_asm_inline)
  AC_LANG_POP([C])
])
