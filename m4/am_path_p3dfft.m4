# SYNOPSIS
#
#   Test for P3DFFT (http://www.sdsc.edu/us/resources/p3dfft/index.php)
#
#   AM_PATH_P3DFFT([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-p3dfft=DIR option. Searches --with-p3dfft,
#   $P3DFFT_HOME, and the usual places for P3DFFT headers and libraries.
#
#   Supports separately specifying --with-p3dfft-include or
#   --with-p3dfft-libdir to override default locations underneath
#   either --with-p3dfft or $P3DFFT_HOME.
#
#   On success, sets P3DFFT_CFLAGS, P3DFFT_LIBS, and
#   #defines HAVE_P3DFFT.  When ACTION-IF-NOT-FOUND is not specified, 
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   2009-02-11
#
# COPYLEFT
#
#   Copyright (c) 2008 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_P3DFFT],
[
AC_ARG_VAR(P3DFFT_HOME,[root directory of P3DFFT installation])

AC_ARG_WITH(p3dfft, [AS_HELP_STRING([--with-p3dfft[=DIR]],[root directory of P3DFFT installation])],[
with_p3dfft=$withval
if test "${with_p3dfft}" != yes; then
    P3DFFT_HOME=$withval
    p3dfft_include="$withval/include"
    p3dfft_libdir="$withval/lib"
fi
],[
with_p3dfft=$withval
if test "x${P3DFFT_HOME}" != "x"; then
    p3dfft_include="${P3DFFT_HOME}/include"
    p3dfft_libdir="${P3DFFT_HOME}/lib"
fi
])

AC_ARG_WITH(p3dfft-include,
[AS_HELP_STRING([--with-p3dfft-include=DIR],[specify exact directory for P3DFFT headers])],[
if test -d "$withval"; then
    p3dfft_include="$withval"
else
    AC_MSG_ERROR([--with-p3dfft-include expected directory name])
fi
])

AC_ARG_WITH(p3dfft-libdir, [AS_HELP_STRING([--with-p3dfft-libdir=DIR],[specify exact directory for P3DFFT libraries])],[
if test -d "$withval"; then
    p3dfft_libdir="$withval"
else
    AC_MSG_ERROR([--with-p3dfft-libdir expected directory name])
fi
])

if test "${with_p3dfft}" != no ; then
    P3DFFT_LIBS="-lp3dfft"
    if test -d "${p3dfft_libdir}" ; then
        P3DFFT_LIBS="-L${p3dfft_libdir} ${P3DFFT_LIBS}"
    fi
    P3DFFT_CFLAGS=""
    if test -d "${p3dfft_include}" ; then
        P3DFFT_CFLAGS="-I${p3dfft_include} ${P3DFFT_CFLAGS}"
    fi

    ac_save_CFLAGS="$CFLAGS"
    ac_save_LIBS="$LIBS"
    CFLAGS="${P3DFFT_CFLAGS} ${CFLAGS}"
    LDFLAGS="${P3DFFT_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([p3dfft.h],[found_header=yes],[found_header=no])
    AC_HAVE_LIBRARY([p3dfft],[found_library=yes],[found_library=no])
    AC_LANG_POP([C])
    LIBS="$ac_save_LIBS"
    CFLAGS="$ac_save_CFLAGS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$found_library" = yes; then
            succeeded=yes
        fi
    fi

    if test "$succeeded" = no; then
        ifelse([$2],,AC_MSG_ERROR([P3DFFT not found.  Try either --with-p3dfft or setting P3DFFT_HOME.]),
            [$2])
    else
        AC_DEFINE(HAVE_P3DFFT,1,[Define if P3DFFT is available])
        AC_SUBST(P3DFFT_CFLAGS)
        AC_SUBST(P3DFFT_LIBS)
        ifelse([$1],,,[$1])
    fi

fi

])
