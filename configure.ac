dl ------------------------------------------------
dnl Package name and version number
dnl ------------------------------------------------
AC_INIT(suzerain, [0.0.1], [rhys@ices.utexas.edu])
AC_REVISION([$Id$])

dnl ------------------------------------------------
dnl Initialization macros
dnl ------------------------------------------------
AC_CONFIG_SRCDIR([suzerain/pencil_grid.hpp])
AM_CONFIG_HEADER([suzerain/config.h])
AX_PREFIX_CONFIG_H([suzerain/suzerain-config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.9 nostdinc no-define tar-ustar -Wall -Werror foreign)
AC_PREREQ([2.61])
AX_AM_MACROS
AX_ADD_AM_CTAGSFLAGS

dnl ------------------------------------------------
dnl LaTeX capabilities
dnl ------------------------------------------------
AX_ADD_AM_LATEX
ACLTX_CLASSES([amsart],AMSART_CLASS,,AC_MSG_WARN([Unable to find LaTeX class 'amsart']))
AC_CACHE_SAVE

dnl ------------------------------------------------------
dnl Check for compiler toolchain and optimization level
dnl Ordering very specific: AX_CC_MAXOPT usage is fragile
dnl ------------------------------------------------------
AX_CHECK_ENABLE_DEBUG
AC_LANG([C])
AC_PROG_CC
AX_COMPILER_VENDOR
if test "x$ax_enable_debug" != "xyes"; then
    dnl Determine good optimization flags for C and store in CFLAGS_MAXOPT
    CFLAGS_pre_AX_CC_MAXOPT=${CFLAGS}
    CFLAGS=
    ac_test_CFLAGS_pre_AX_CC_MAXOPT=${ac_test_CFLAGS}
    ac_test_CFLAGS=
    AX_CC_MAXOPT
    ac_test_CFLAGS=${ac_test_CFLAGS_pre_AX_CC_MAXOPT}
    CFLAGS_MAXOPT=${CFLAGS}
    CFLAGS="${CFLAGS_pre_AX_CC_MAXOPT} ${CFLAGS_MAXOPT}"
    dnl Presume optimization flags useful for other languages
    CXXFLAGS="${CXXFLAGS} ${CFLAGS_MAXOPT}"
    FFLAGS="${FFLAGS} ${CFLAGS_MAXOPT}"
    FCFLAGS="${FCFLAGS} ${CFLAGS_MAXOPT}"
fi
AC_LANG([C++])
AC_PROG_CXX
AX_COMPILER_VENDOR
AC_PROG_RANLIB
AC_PROG_SED
AX_AM_JOBSERVER([no])
AC_CACHE_SAVE

dnl ----------------------------------------------
dnl Enable rudimentary precompiled header support
dnl Currently very, very Intel-specific
dnl ----------------------------------------------
dnl AC_LANG([C])
dnl AX_CHECK_COMPILER_FLAGS([-pch],[CFLAGS="${CFLAGS} -pch"])
dnl AC_LANG([C++])
dnl AX_CHECK_COMPILER_FLAGS([-pch],[CXXFLAGS="${CXXFLAGS} -pch"])

dnl ------------------------------------------------
dnl Enable GNU libtool
dnl Current version checks also in Makefile.am
dnl ------------------------------------------------
LT_INIT
LT_PREREQ([2.2.6])
AC_SUBST(LIBTOOL_DEPS)
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C library requirements
dnl ------------------------------------------------
AC_LANG_PUSH([C])
AX_PROG_CC_C99_CFLAGS(,AC_MSG_ERROR([Compiler cannot support ISO C99]))
ACX_MKL()
export GSL_CBLAS_LIB="${MKL_LIBS} ${FLIBS}"
AX_PATH_GSL(1.13,,AC_MSG_ERROR([Could not find required GSL version.]))
dnl --------BEGIN FFTW3 checks---------------------
PKG_CHECK_MODULES(FFTW3, fftw3 > 3.2)
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)
found_fftw3_threads=no
if test "${found_fftw3_threads}" != yes; then
AX_OPENMP([
    AC_DEFINE(HAVE_OPENMP,[1],[Define if OpenMP is enabled])
    AC_CHECK_LIB([fftw3_threads],[fftw_init_threads],
        [
          AC_DEFINE([HAVE_FFTW3_THREADS],[1],[Defined if FFTW3 threads available])
          FFTW3_THREADS_CFLAGS="$OPENMP_CFLAGS $FFTW3_CFLAGS"
          FFTW3_THREADS_LIBS="-lfftw3_threads $FFTW3_LIBS $OPENMP_CFLAGS"
          AC_MSG_NOTICE([Detected OpenMP-enabled FFTW threading])
          found_fftw3_threads=yes
        ],
        [AC_MSG_NOTICE([Detected OpenMP but not OpenMP-enabled FFTW threading])],
        [$FFTW3_LIBS $OPENMP_CFLAGS])
    ])
fi
if test "${found_fftw3_threads}" != yes; then
ACX_PTHREAD([
    AC_DEFINE(HAVE_PTHREAD,[1],
              [Define if you have POSIX threads libraries and header files])
    AC_CHECK_LIB([fftw3_threads],[fftw_init_threads],
        [
          AC_DEFINE([HAVE_FFTW3_THREADS],[1],[Defined if FFTW3 threads available])
          FFTW3_THREADS_CFLAGS="$PTHREAD_CFLAGS $FFTW3_CFLAGS"
          FFTW3_THREADS_LIBS="-lfftw3_threads $FFTW3_LIBS $PTHREAD_LIBS"
          AC_MSG_NOTICE([Detected pthread-enabled FFTW threading])
          found_fftw3_threads=yes
        ],
        [AC_MSG_NOTICE([Detected pthread but not pthread-enabled FFTW threading])],
        [$FFTW3_LIBS $PTHREAD_LIBS])
    ])
fi
if test "${found_fftw3_threads}" != yes; then
  FFTW3_THREADS_CFLAGS=$FFTW3_CFLAGS
  FFTW3_THREADS_LIBS=$FFTW3_LIBS
  AC_MSG_NOTICE([FFTW threading not detected])
fi
AC_SUBST(FFTW3_THREADS_CFLAGS)
AC_SUBST(FFTW3_THREADS_LIBS)
dnl --------END FFTW3 checks-----------------------
dnl END FFTW3 checks
PKG_CHECK_MODULES(P3DFFT, suzerain-p3dfft)
AC_SUBST(P3DFFT_CFLAGS)
AC_SUBST(P3DFFT_LIBS)
PKG_CHECK_MODULES(HPCT,
                  [hpct >= 0.26],
                  AC_DEFINE([HAVE_HPCT],[1],[Defined if HPC Toolkit found]),
                  AC_MSG_WARN([Could not find HPC Toolkit; Continuing...]))
AC_SUBST(HPCT_CFLAGS)
AC_SUBST(HPCT_LIBS)
AC_LANG_POP([C])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C library requirements
dnl ------------------------------------------------


dnl ------------------------------------------------
dnl C++ library requirements
dnl ------------------------------------------------
AC_LANG_PUSH([C++])
PKG_CHECK_MODULES(LOG4CXX,
    [liblog4cxx >= 0.10 apr-1 >= 1.3.3 apr-util-1 >= 1.3.4],
    AC_DEFINE([HAVE_LOG4CXX],1,[Define if log4cxx is available])
    AC_SUBST(LOG4CXX_CFLAGS)
    AC_SUBST(LOG4CXX_LIBS))
BOOST_REQUIRE([1.39])
BOOST_FOREACH
BOOST_FORMAT
BOOST_IOSTREAMS
BOOST_MULTIARRAY
BOOST_NUMERIC_CONVERSION
BOOST_PREPROCESSOR
BOOST_PROGRAM_OPTIONS
BOOST_STATICASSERT
BOOST_TEST
BOOST_TYPETRAITS
PKG_CHECK_MODULES(EIGEN,
    [eigen2 >= 2.0.5],
    AC_DEFINE([HAVE_EIGEN],1,[Define if Eigen is available])
    dnl EIGEN_LIBS should be empty and is so no AC_SUBST is used
    AC_SUBST(EIGEN_CFLAGS))
AC_LANG_POP([C++])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C++ compilation strictness and compatibility
dnl ------------------------------------------------
dnl Strictness must come after libraries to avoid breaking AC_CHECK_LIB
AC_CXX_HAVE_STL
AC_CXX_HEADER_STDCXX_98
AC_CXX_HEADER_STDCXX_TR1
AC_CXX_COMPILE_STDCXX_0X
AC_CXX_CPPFLAGS_STD_LANG(AM_CPPFLAGS)
AC_CXX_CXXFLAGS_STD_LANG(AM_CXXFLAGS)
AC_CXX_LDFLAGS_STD_LANG(AM_LDFLAGS)
dnl Need to determine how to enable some but not of Intel's ICPC babble
dnl AX_CFLAGS_WARN_ALL_ANSI(CFLAGS)
dnl AX_CFLAGS_WARN_ALL_ANSI(CXXFLAGS)
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl MPI toolchain checks
dnl ------------------------------------------------
AC_LANG_PUSH([C])
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for C, e.g. mpicc.))
AX_MPIP()
AC_LANG_POP([C])
AC_LANG_PUSH([C++])
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for C++, e.g. mpicxx.))
AC_LANG_POP([C++])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Documentation generation
dnl ------------------------------------------------
DX_DOXYGEN_FEATURE([ON])
DX_HTML_FEATURE([ON])
DX_CHM_FEATURE([OFF])
DX_CHI_FEATURE([OFF])
DX_MAN_FEATURE([OFF])
DX_RTF_FEATURE([OFF])
DX_XML_FEATURE([OFF])
DX_PDF_FEATURE([OFF])
DX_PS_FEATURE([OFF])
DX_INIT_DOXYGEN([suzerain],[doxygen.cfg],[docs])

dnl ------------------------------------------------
dnl Generate Makefiles
dnl ------------------------------------------------
AC_CONFIG_FILES([
    Makefile
    apps/Makefile
    writeup/Makefile
    lib/Makefile
    suzerain/Makefile
    tests/Makefile
])

AC_OUTPUT()
