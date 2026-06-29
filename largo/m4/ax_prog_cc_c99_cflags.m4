# ===========================================================================
#         http://autoconf-archive.cryp.to/ax_prog_cc_c99_cflags.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PROG_CC_C99_CFLAGS([ACTION-IF-AVAILABLE], [ACTION-IF-UNAVAILABLE])
#
# DESCRIPTION
#
#   Ensure the C compiler operates in C99 (or later) mode.
#
#   Historically this macro wrapped the now-obsolete AC_PROG_CC_C99 in order
#   to append any required standard-selection flag to CFLAGS rather than to
#   CC.  As of Autoconf 2.70, AC_PROG_CC itself selects the newest supported
#   C standard and folds any required option directly into CC, so no extra
#   handling is needed here.  This macro now simply requires AC_PROG_CC and
#   runs ACTION-IF-AVAILABLE.  ACTION-IF-UNAVAILABLE is retained for
#   backwards compatibility but is unreachable: AC_PROG_CC aborts configure
#   earlier if no working, sufficiently modern C compiler is found.
#
# LAST MODIFICATION
#
#   2026-06-28
#
# COPYLEFT
#
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PROG_CC_C99_CFLAGS],[
    AC_REQUIRE([AC_PROG_CC])
    m4_default([$1],[:])
])
