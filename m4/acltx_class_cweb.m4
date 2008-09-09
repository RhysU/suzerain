# ===========================================================================
#              http://autoconf-archive.cryp.to/acltx_class_cweb.html
# ===========================================================================
#
# SYNOPSIS
#
#   ACLTX_CLASS_CWEB(VARIABLETOSET[,ACTION-IF-FOUND[,ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macros test if the class cweb exists and works.  It sets
#   VARIABLETOSET to yes or no If ACTION-IF-FOUND (and ACTION-IF-NOT-FOUND)
#   are set, do the correct action.
#
#   The cweb package is used to provide LaTeX support atop Knuth's CWEB
#   literate programming environment.  The test LaTeX document requires
#   some additional logic beyond ACLTX_CLASS because the cweb package expects
#   cweave to insert some boilerplate.  Both macros for the original
#   CWEB and the compatible offshoot CWEBx are tested here.  See
#   http://www.ctan.org/tex-archive/help/Catalogue/entries/cweb-latex.html.
#
# LAST MODIFICATION
#
#   2008-08-03
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

AC_DEFUN([ACLTX_CLASS_CWEB],[

ac_cv_latex_class_cweb=no
ACLTX_PACKAGE_LOCATION(cweb.cls,ac_cv_latex_class_cweb_location)

if test -z "$ac_cv_latex_class_cweb_location" ; then
    AC_MSG_WARN([Unable to locate the cweb.cls file])
else
    ACLTX_PACKAGE_LOCATION(cwebmac.tex,ac_cv_latex_class_cweb_cwebmac)
    if test -n "$ac_cv_latex_class_cweb_cwebmac"; then
        AC_CACHE_CHECK(
            [for usability of class cweb with Levy/Knuth CWEB],
            [ac_cv_latex_class_cweb_CWEB],[
            _ACLTX_TEST([
            \input cwebmac
            \documentclass{cweb}
            \begin{document}
            \M{1}
            \end{document}
            \fi
            \fin
            \end dnl Additional ending commands required for
            \fi  dnl pdfeTeX, Version 3.141592-1.21a-2.2 (Web2C 7.5.4).
            \fin dnl Runs cleanly without on versions like 1.4x
            ], [ac_cv_latex_class_cweb_CWEB])
        ])
        if test "$ac_cv_latex_class_cweb_CWEB" = yes; then
            ac_cv_latex_class_cweb=yes
        fi
    fi
    ACLTX_PACKAGE_LOCATION(cwebcmac.tex,ac_cv_latex_class_cweb_cwebcmac)
    if test -n "$ac_cv_latex_class_cweb_cwebcmac"; then
        AC_CACHE_CHECK(
            [for usability of class cweb with van Leeuwen CWEBx],
            [ac_cv_latex_class_cweb_CWEBx],[
            _ACLTX_TEST([
            \input cwebcmac
            \documentclass{cweb}
            \begin{document}
            \M{1}
            \end{document}
            \fi
            \fin
            \end dnl Additional ending commands required for
            \fi  dnl pdfeTeX, Version 3.141592-1.21a-2.2 (Web2C 7.5.4).
            \fin dnl Runs cleanly without on versions like 1.4x
            ],[ac_cv_latex_class_cweb_CWEBx])
        ])
        if test "$ac_cv_latex_class_cweb_CWEBx" = yes; then
            ac_cv_latex_class_cweb=yes
        fi
    fi
fi

$1=$[ac_cv_latex_class_cweb] ; export $1
AC_SUBST($1)

if test "$[ac_cv_latex_class_cweb]" = yes ; then
    dnl NOP required
    :
    ifelse([$2], , ,[$2])
else
    dnl NOP required
    :
    ifelse([$3], , ,[$3])
fi
])
