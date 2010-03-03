#
#   AX_FFTW3([FFTW3_MODULES])
#   where FFTW3_MODULES is a PKG_CHECK_MODULES string like [fftw3 > 3.2].
#
AC_DEFUN([AX_FFTW3], [
AC_PREREQ(2.60)
AC_REQUIRE([AC_PROG_SED])
AC_LANG_PUSH([C])

PKG_CHECK_MODULES(FFTW3,m4_default([$1],[fftw3 > 3.2]))
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)

ax_fftw3_link_pre=`echo  $FFTW3_LIBS | $SED -e 's/-lfftw3.*//'`
ax_fftw3_link_post=`echo $FFTW3_LIBS | $SED -e 's/^.*-lfftw3/-lfftw3/'`
ax_fftw3_found_threads=no

# Check if we can link in threading without any additional help
AC_CHECK_LIB([fftw3_threads],[fftw_init_threads],[
    ax_fftw3_found_threads=yes 
    FFTW3_THREADS_CFLAGS="$OPENMP_CFLAGS $FFTW3_CFLAGS"
    FFTW3_THREADS_LIBS="$ax_fftw3_link_pre -lfftw3_threads $ax_fftw3_link_post"
    AC_MSG_NOTICE([Detected FFTW threading not requiring any external dependencies])
])

# Check if OpenMP is the way to go...
if test "${ax_fftw3_found_threads}" != yes; then
    ax_fftw3_save_LDFLAGS=$LDFLAGS
    LDFLAGS="$LDFLAGS $ax_fftw3_link_pre"
    AX_OPENMP([
        AC_DEFINE(HAVE_OPENMP,[1],[Define if OpenMP is enabled])
        # Note the new function name to thwart caching from earlier AC_CHECK_LIB
        AC_CHECK_LIB([fftw3_threads],[fftw_cleanup_threads], [
                ax_fftw3_found_threads=yes
                FFTW3_THREADS_CFLAGS="$OPENMP_CFLAGS $FFTW3_CFLAGS"
                FFTW3_THREADS_LIBS="$ax_fftw3_link_pre -lfftw3_threads $ax_fftw3_link_post $OPENMP_CFLAGS"
                AC_MSG_NOTICE([Detected OpenMP-enabled FFTW threading])
            ], [
                AC_MSG_NOTICE([Detected OpenMP but not OpenMP-enabled FFTW threading])
            ], [$ax_fftw3_link_post $OPENMP_CFLAGS]
        )
    ])
    LDFLAGS=$ax_fftw3_save_LDFLAGS
fi

# Check if POSIX Threads is the way to go...
if test "${ax_fftw3_found_threads}" != yes; then
    ax_fftw3_save_LDFLAGS=$LDFLAGS
    LDFLAGS="$LDFLAGS $ax_fftw3_link_pre"
    ACX_PTHREAD([
        AC_DEFINE(HAVE_PTHREAD,[1],
                [Define if you have POSIX threads libraries and header files])
        # Note the new function name to thwart caching from earlier AC_CHECK_LIB
        AC_CHECK_LIB([fftw3_threads],[fftw_plan_with_nthreads], [
            ax_fftw3_found_threads=yes
            FFTW3_THREADS_CFLAGS="$PTHREAD_CFLAGS $FFTW3_CFLAGS"
            FFTW3_THREADS_LIBS="$ax_fftw3_link_pre -lfftw3_threads $ax_fftw3_link_post $PTHREAD_LIBS"
            AC_MSG_NOTICE([Detected pthread-enabled FFTW threading])
        ], [
            AC_MSG_NOTICE([Detected pthread but not pthread-enabled FFTW threading])
        ], [$ax_fftw3_link_post $PTHREAD_LIBS])
    ])
    LDFLAGS=$ax_fftw3_save_LDFLAGS
fi

if test "${ax_fftw3_found_threads}" = yes; then
    AC_DEFINE([HAVE_FFTW3_THREADS],[1],[Defined if FFTW3 threads available])
else
    # Give up and provide serial-only versions
    FFTW3_THREADS_CFLAGS=$FFTW3_CFLAGS
    FFTW3_THREADS_LIBS=$FFTW3_LIBS
fi
AC_SUBST(FFTW3_THREADS_CFLAGS)
AC_SUBST(FFTW3_THREADS_LIBS)

# TODO Complete fftw3_found_mpi checks
ax_fftw3_found_mpi=no
AC_LANG_POP([C])
])dnl AX_FFTW3
