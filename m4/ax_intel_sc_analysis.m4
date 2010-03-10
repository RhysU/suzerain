# SYNOPSIS
#
#   AX_INTEL_SC_ANALYSIS([Default analysis level in {0,1,2,3}])
#
# DESCRIPTION
#
#   Provide a configure-time option to execute Intel's source code analysis
#   (http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/cpp/lin/compiler_c/bldaps_cls/common/bldaps_svover.htm)
#   at build time.  An analysis level of 0 disables this diagnostic.  An
#   analysis level of {1,2,3} corresponds to '-diag-enable s{1,2,3}'.  If not
#   provided, the analysis level defaults to 0.  If enabled without any value
#   given, the analysis level used is 2.  When enabled, CFLAGS, CXXFLAGS, and
#   LDFLAGS are modified appropriately.
#
#   Requires macros: AX_COMPILER_VENDOR, AX_CHECK_COMPILER_FLAGS
#
# LAST MODIFICATION
#
#   2010-03-10
#
# COPYLEFT
#
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
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
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_INTEL_SC_ANALYSIS],
[

dnl Process the default level
m4_define([ax_intel_sc_analysis_default],[m4_tolower(m4_normalize(ifelse([$1],,[0],[$1])))])
case ax_intel_sc_analysis_default in #(
  "0" | "1" | "2" | "3")
     ;;#(
  *) AC_MSG_ERROR([Incorrect default value in AX_[]INTEL_SC_ANALYSIS])
     ;;
esac

dnl Process the incoming argument
AC_ARG_ENABLE(sc-analysis,AS_HELP_STRING([--enable-sc-analysis],[Source code analysis level in {0,1,2,3} @<:@default=ax_intel_sc_analysis_default@:>@]),[
    case $enableval in #(
      "no")
          ax_intel_sc_analysis=0
          ;;#(
      "yes")
          case ax_intel_sc_analysis_default in
            0) ax_intel_sc_analysis=2
               ;;#(
            *) ax_intel_sc_analysis=ax_intel_sc_analysis_default
               ;;
          esac
          ;;#(
      "0" | "1" | "2" | "3")
          ax_intel_sc_analysis=$enableval
          ;;#(
      *)  AC_MSG_WARN([Unknown --enable-sc-analysis value @<:@$enableval@:>@])
          ax_intel_sc_analysis=ax_intel_sc_analysis_default
          ;;
    esac
  ],[
    ax_intel_sc_analysis=ax_intel_sc_analysis_default
  ])

dnl Process the incoming flag
AC_MSG_CHECKING(whether to enable source code analysis)
case $ax_intel_sc_analysis in #(
  "1" | "2" | "3")
     ax_intel_sc_analysis_flag="-diag-enable sv$ax_intel_sc_analysis"
     AC_MSG_RESULT([yes at level $ax_intel_sc_analysis])
     AC_LANG_PUSH([C])
     case $ax_cv_c_compiler_vendor in #(
      intel) AX_CHECK_COMPILER_FLAGS([$ax_intel_sc_analysis_flag], [
               # Assume working for C implies it works for C++
               CFLAGS="$CFLAGS $ax_intel_sc_analysis_flag"
               CXXFLAGS="$CXXFLAGS $ax_intel_sc_analysis_flag"
               LDFLAGS="$LDFLAGS $ax_intel_sc_analysis_flag"
             ])
             ;;#(
       *)    AC_MSG_WARN([Cannot enable source code analysis for C compiler vendor $ax_cv_c_compiler_vendor])
             ;;
     esac
     AC_LANG_POP()
     ;;#(
  *) AC_MSG_RESULT([no])
     ;;
esac
])
