# SYNOPSIS
#
#   AX_ASM_INLINE_COMMENT
#
# DESCRIPTION
#
#   If possible, create the ASM_INLINE_COMMENT(COMMENT) macro which generates
#   an inline assembly comment like "__FILE__:__LINE__ COMMENT".  Useful in
#   conjunction with compile-to-assembly compiler options as a way to check the
#   assembly output within a critical code region.
#
# LICENSE
#
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

AC_DEFUN([AX_ASM_INLINE_COMMENT], [
  AC_REQUIRE([AX_ASM_INLINE])
  AC_LANG_PUSH([C])
  AC_MSG_CHECKING(for inline assembly comment style)
  AC_CACHE_VAL(ac_cv_asm_inline_comment, [
    if test "$ac_cv_asm_inline" != "none"; then
        ax_asm_inline_comment_keywords="# ; ;; none"
    else
        ax_asm_inline_comment_keywords="none"
    fi
    for ax_asm_inline_comment_keyword in $ax_asm_inline_comment_keywords; do
       case $ax_asm_inline_comment_keyword in
          none) ac_cv_asm_inline_comment=none ; break ;;
      *)
         AC_TRY_COMPILE(
            [#include <stdlib.h>
             static void
             foo(void) {
             ]$ax_asm_inline_keyword[("]$ax_asm_inline_comment_keyword[!?!$$!?!");
             exit(1);
             }],
             [],
             [ac_cv_asm_inline_comment=$ax_asm_inline_comment_keyword ; break],
             ac_cv_asm_inline_comment=none
         )
      esac
    done
])

  if test "$ac_cv_asm_inline_comment" != "none"; then
    ax_asm_inline_comment_macro="$ac_cv_asm_inline(\"$ac_cv_asm_inline_comment \"__FILE__##\":\"##ASM_INLINE_COMMENT_HELPER_TOSTRING(__LINE__)##\" \"## #comment)"
  else
    ax_asm_inline_comment_macro="/* unsupported */"
  fi
  AC_DEFINE_UNQUOTED([ASM_INLINE_COMMENT(comment)],
            [$ax_asm_inline_comment_macro],
            [If supported, define a macro to use inline assembly comments])
  AC_DEFINE([ASM_INLINE_COMMENT_HELPER_TOSTRING(x)],
            [ASM_INLINE_COMMENT_HELPER_STRINGIFY(x)],
            [Helps ASM_INLINE_COMMENT convert __LINE__ to a (char *)])
  AC_DEFINE([ASM_INLINE_COMMENT_HELPER_STRINGIFY(x)],
            [#x],
            [Helps ASM_INLINE_COMMENT convert __LINE__ to a (char *)])
  AC_MSG_RESULT($ac_cv_asm_inline_comment)
  AC_LANG_POP([C])
])
