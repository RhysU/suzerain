#
#   AX_FFTW3([FFTW3_MODULES])
#   where FFTW3_MODULES is a PKG_CHECK_MODULES string like [fftw3 >= 3.3].
#   This macro behaves poorly on pre-FFTW 3.3 final releases.
#
AC_DEFUN([AX_FFTW3], [
AC_PREREQ(2.60)
AC_REQUIRE([AC_PROG_SED])
AC_LANG_PUSH([C])

# Search for FFTW itself using pkg-config
PKG_CHECK_MODULES(FFTW3,m4_default([$1],[fftw3 >= 3.3]))
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)

dnl Note that FFTW3_THREADS_{CFLAGS,LIBS} and FFTW3_MPI_THREADS_{CFLAGS,LIBS}
dnl are always AC_SUBSTituted (the latter whenever MPI is detected).  This is
dnl so applications wanting threading if it is available can simply use these
dnl with an automatic fallback to non-threaded compilation if no threading is
dnl possible.

# Prepare for FFTW threading search
ax_fftw3_libs_pre=`echo  $FFTW3_LIBS | $SED -e 's/-lfftw3.*$//'`
ax_fftw3_libs_post=`echo $FFTW3_LIBS | $SED -e 's/^.*-lfftw3/-lfftw3/'`
# Default values overridden on success
FFTW3_THREADS_CFLAGS=$FFTW3_CFLAGS
FFTW3_THREADS_LIBS=$FFTW3_LIBS
ax_fftw3_found_threads=no
ax_fftw3_found_openmp=no
ax_fftw3_found_pthreads=no

dnl Notice that POSIX Threads trumps OpenMP and will used as the default if both
dnl are detected.  That choice was made because autoconf/libtool are more robust
dnl in their handling of pthread-related flags versus OpenMP-related flags.

# Check if we can link in OpenMP-based threading without any additional help...
AC_CHECK_LIB([fftw3_omp],[fftw_init_threads],[
    ax_fftw3_found_openmp=yes
    FFTW3_OPENMP_CFLAGS="$OPENMP_CFLAGS $FFTW3_CFLAGS"
    FFTW3_OPENMP_LIBS="$ax_fftw3_libs_pre -lfftw3_omp $ax_fftw3_libs_post"
    AC_SUBST(FFTW3_OPENMP_CFLAGS)
    AC_SUBST(FFTW3_OPENMP_LIBS)
    AC_MSG_NOTICE([Detected FFTW OpenMP threading not requiring any external dependencies])
    FFTW3_THREADS_CFLAGS="$FFTW3_OPENMP_CFLAGS"
    FFTW3_THREADS_LIBS="$FFTW3_OPENMP_LIBS"
])

# ...or if explicitly requesting OpenMP is required.
if test "$ax_fftw3_found_openmp" != yes; then
    ax_fftw3_save_LDFLAGS=$LDFLAGS
    LDFLAGS="$LDFLAGS $ax_fftw3_libs_pre"
    AX_OPENMP([
        AC_DEFINE(HAVE_OPENMP,[1],[Define if OpenMP is enabled])
        # Note the new function name to thwart caching from earlier AC_CHECK_LIB
        AC_CHECK_LIB([fftw3_omp],[fftw_cleanup_threads], [
                ax_fftw3_found_openmp=yes
                FFTW3_OPENMP_CFLAGS="$OPENMP_CFLAGS $FFTW3_CFLAGS"
                FFTW3_OPENMP_LIBS="$ax_fftw3_libs_pre -lfftw3_omp $ax_fftw3_libs_post $OPENMP_CFLAGS"
                AC_SUBST(FFTW3_OPENMP_CFLAGS)
                AC_SUBST(FFTW3_OPENMP_LIBS)
                AC_MSG_NOTICE([Detected OpenMP-enabled FFTW threading])
                FFTW3_THREADS_CFLAGS="$FFTW3_OPENMP_CFLAGS"
                FFTW3_THREADS_LIBS="$FFTW3_OPENMP_LIBS"
            ], [
                AC_MSG_NOTICE([Detected OpenMP but not OpenMP-enabled FFTW threading])
            ], [$ax_fftw3_libs_post $OPENMP_CFLAGS]
        )
    ])
    LDFLAGS=$ax_fftw3_save_LDFLAGS
fi

# Check if POSIX Threads is another possible way to go.
ax_fftw3_save_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS $ax_fftw3_libs_pre"
ACX_PTHREAD([
    AC_DEFINE(HAVE_PTHREAD,[1],
            [Define if you have POSIX threads libraries and header files])
    # Note the new function name to thwart caching from earlier AC_CHECK_LIB
    AC_CHECK_LIB([fftw3_threads],[fftw_plan_with_nthreads], [
        ax_fftw3_found_pthreads=yes
        FFTW3_PTHREADS_CFLAGS="$PTHREAD_CFLAGS $FFTW3_CFLAGS"
        FFTW3_PTHREADS_LIBS="$ax_fftw3_libs_pre -lfftw3_threads $ax_fftw3_libs_post $PTHREAD_LIBS"
        AC_SUBST(FFTW3_PTHREADS_CFLAGS)
        AC_SUBST(FFTW3_PTHREADS_LIBS)
        AC_MSG_NOTICE([Detected pthread-enabled FFTW threading])
        if test "$ax_fftw3_found_openmp" != yes; then
            AC_MSG_NOTICE([Both pthreads and OpenMP FFTW threading found; pthreads will be the default])
        fi
        FFTW3_THREADS_CFLAGS="$FFTW3_PTHREADS_CFLAGS"
        FFTW3_THREADS_LIBS="$FFTW3_PTHREADS_LIBS"
    ], [
        AC_MSG_NOTICE([Detected pthreads but not pthread-enabled FFTW threading])
    ], [$ax_fftw3_libs_post $PTHREAD_LIBS])
])
LDFLAGS=$ax_fftw3_save_LDFLAGS

# Finish up FFTW threading details
if test "$ax_fftw3_found_pthreads" = yes -o "$ax_fftw3_found_openmp" = yes; then
    AC_DEFINE([HAVE_FFTW3_THREADS],[1],[Defined if FFTW3 threads available])
    ax_fftw3_found_threads=yes
fi
AC_SUBST(FFTW3_THREADS_CFLAGS)
AC_SUBST(FFTW3_THREADS_LIBS)
AM_CONDITIONAL([HAVE_FFTW3_THREADS],  [test "$ax_fftw3_found_threads"  = yes])
AM_CONDITIONAL([HAVE_FFTW3_PTHREADS], [test "$ax_fftw3_found_pthreads" = yes])
AM_CONDITIONAL([HAVE_FFTW3_OPENMP],   [test "$ax_fftw3_found_openmp"   = yes])

# Search for FFTW MPI independent of FFTW threading result
ax_fftw3_found_mpi=no
AC_REQUIRE([ACX_MPI])
ax_fftw3_save_CC=$CC
CC=$MPICC
ax_fftw3_save_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS $ax_fftw3_libs_pre"
AC_CHECK_LIB([fftw3_mpi],[fftw_mpi_init], [
    ax_fftw3_found_mpi=yes
    AC_DEFINE([HAVE_FFTW3_MPI],[1],[Defined if FFTW3 MPI available])

    # Prepare libraries for FFTW MPI only
    FFTW3_MPI_LIBS="$ax_fftw3_libs_pre -lfftw3_mpi $ax_fftw3_libs_post"
    AC_SUBST(FFTW3_MPI_LIBS)

    # Prepare libraries for FFTW MPI with possibly both threading options
    if test "$ax_fftw3_found_openmp" = yes; then
        FFTW3_MPI_OPENMP_CFLAGS="$FFTW3_OPENMP_CFLAGS"
        ax_fftw3_openmp_libs_pre=`echo  $FFTW3_OPENMP_LIBS | $SED -e 's/-lfftw3_omp.*$//'`
        ax_fftw3_openmp_libs_post=`echo $FFTW3_OPENMP_LIBS | $SED -e 's/^.*-lfftw3_omp/-lfftw3_omp/'`
        FFTW3_MPI_OPENMP_LIBS="$ax_fftw3_openmp_libs_pre -lfftw3_mpi $ax_fftw3_openmp_libs_post"
        AC_SUBST(FFTW3_MPI_OPENMP_CFLAGS)
        AC_SUBST(FFTW3_MPI_OPENMP_LIBS)
    fi
    if test "$ax_fftw3_found_pthreads" = yes; then
        FFTW3_MPI_PTHREADS_CFLAGS="$FFTW3_PTHREADS_CFLAGS"
        ax_fftw3_pthreads_libs_pre=`echo  $FFTW3_PTHREADS_LIBS | $SED -e 's/-lfftw3_threads.*$//'`
        ax_fftw3_pthreads_libs_post=`echo $FFTW3_PTHREADS_LIBS | $SED -e 's/^.*-lfftw3_threads/-lfftw3_threads/'`
        FFTW3_MPI_PTHREADS_LIBS="$ax_fftw3_pthreads_libs_pre -lfftw3_mpi $ax_fftw3_pthreads_libs_post"
        AC_SUBST(FFTW3_MPI_OPENMP_CFLAGS)
        AC_SUBST(FFTW3_MPI_PTHREADS_LIBS)
    fi
    # Again, POSIX threads trumps OpenMP as the default if both are detected
    if test "$ax_fftw3_found_pthreads" = yes; then
        FFTW3_MPI_THREADS_CFLAGS="$FFTW3_MPI_PTHREADS_CFLAGS"
        FFTW3_MPI_THREADS_LIBS="$FFTW3_MPI_PTHREADS_LIBS"
    elif test "$ax_fftw3_found_openmp" = yes; then
        FFTW3_MPI_THREADS_CFLAGS="$FFTW3_MPI_OPENMP_CFLAGS"
        FFTW3_MPI_THREADS_LIBS="$FFTW3_MPI_OPENMP_LIBS"
    else
        FFTW3_MPI_THREADS_CFLAGS=$FFTW3_CFLAGS
        FFTW3_MPI_THREADS_LIBS="$ax_fftw3_libs_pre -lfftw3_mpi $ax_fftw3_libs_post"
    fi
    AC_SUBST(FFTW3_MPI_THREADS_CFLAGS)
    AC_SUBST(FFTW3_MPI_THREADS_LIBS)
], [
    AC_MSG_ERROR([Could not find FFTW MPI functionality])
], [$ax_fftw3_libs_post])
AM_CONDITIONAL([HAVE_FFTW3_MPI],  [test "$ax_fftw3_found_mpi" = yes])
LDFLAGS=$ax_fftw3_save_LDFLAGS
CC=$ax_fftw3_save_CC

AC_LANG_POP([C])
])dnl AX_FFTW3
