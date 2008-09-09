# ===========================================================================
#         http://autoconf-archive.cryp.to/ac_lib_hdf5.html
# ===========================================================================
#
# SYNOPSIS
#
#   Test for the HDF5 (http://hdf.ncsa.uiuc.edu/HDF5/) library.
#
#   AC_LIB_HDF5([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   On success, #defines HAVE_LIBHDF5 due to AC_CHECK_LIB invocation.
#
# LAST MODIFICATION
#
#   2008-08-06
#
# COPYLEFT
#
#   Copyright (c) 2008 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AC_LIB_HDF5],[
    ac_lib_hdf5=yes
    AC_CHECK_HEADER([hdf5.h],,     [ac_lib_hdf5=no])
    AC_CHECK_LIB(   [hdf5],[main],,[ac_lib_hdf5=no])
    if test "$ac_lib_hdf5" = yes; then
        dnl NOP required
        : 
		ifelse([$1], , , [$1])
    else
        dnl NOP required
        :
		ifelse([$2], , , [$2])
    fi
])
