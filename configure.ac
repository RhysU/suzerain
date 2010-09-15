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
AX_CHECK_ENABLE_DEBUG([info])
AC_LANG([C])
AC_PROG_CC
AX_PROG_CC_C99_CFLAGS(,AC_MSG_ERROR([Compiler cannot support ISO C99]))
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
dnl dnl Enable IPO if the (Intel) compiler supports it, assumes that
dnl dnl C/C++/Fortran dnl all support IPO if the compiler supports it.
dnl dnl Assumption required because we must replace linker ar with xiar.
dnl AC_LANG_PUSH([C])
dnl AX_CHECK_COMPILER_FLAGS([-ipo], [
dnl     AR="xiar"
dnl     CFLAGS="$CFLAGS -ipo"
dnl     CXXFLAGS="$CXXFLAGS -ipo"
dnl     FFLAGS="$FFLAGS -ipo"
dnl ])
dnl AC_LANG_POP()
fi
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AC_LANG([C++])
AC_PROG_CXX
AX_COMPILER_VENDOR
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AC_LANG([C++])
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
dnl C library requirements
dnl ------------------------------------------------
AC_LANG_PUSH([C])
ACX_MKL()
export GSL_CBLAS_LIB="${MKL_LIBS} ${FLIBS}"
AX_PATH_GSL(1.13,,AC_MSG_ERROR([Could not find required GSL version.]))
AX_FFTW3()
PKG_CHECK_MODULES(P3DFFT, suzerain-p3dfft)
AC_SUBST(P3DFFT_CFLAGS)
AC_SUBST(P3DFFT_LIBS)
dnl TODO Update for post 0.27 release with corrected pkg-config
dnl This check is currently a mess to avoid a broken name change in HPCT
PKG_CHECK_MODULES(HPCT_DONOTUSE, [hpct >= 0.26], [
    AC_CHECK_LIB([hpct],[hpct_timer_init],[
        AC_CHECK_HEADER([hpct.h],[
            AC_DEFINE([HAVE_HPCT],[1],[Defined if HPC Toolkit found])
            AC_SUBST(HPCT_CFLAGS,[${HPCT_DONOTUSE_CFLAGS}])
            AC_SUBST(HPCT_LIBS,[${HPCT_DONOTUSE_LIBS}])
        ],[
            AC_MSG_WARN([Found HPCT's pkg-config file but not hpct.h; Continuing...])
        ],[${HPCT_DONOTUSE_CFLAGS}])
    ],[
        AC_MSG_WARN([Found HPCT's pkg-config file but not -lhpct; Continuing...])
    ],[${HPCT_DONOTUSE_LIBS}])
], [
    AC_MSG_WARN([Could not find HPCT pkg-config file; Continuing...])
])
AC_LANG_POP([C])
AC_CACHE_SAVE

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
BOOST_BIND
BOOST_FOREACH
BOOST_FORMAT
BOOST_FUNCTION
BOOST_IOSTREAMS
BOOST_MATH
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
dnl C/C++ compilation strictness and compatibility
dnl ------------------------------------------------
dnl Strictness must come after libraries to avoid breaking AC_CHECK_LIB
AX_CFLAGS_WARN_ALL(CFLAGS)
AX_CFLAGS_WARN_ALL(CXXFLAGS)
AX_WARNINGS_SANITIZE
AX_INTEL_SC_ANALYSIS
AC_CXX_HAVE_STL
AC_CXX_HEADER_STDCXX_98
AC_CXX_HEADER_STDCXX_TR1
AC_CXX_COMPILE_STDCXX_0X
AC_CXX_CPPFLAGS_STD_LANG(AM_CPPFLAGS)
AC_CXX_CXXFLAGS_STD_LANG(AM_CXXFLAGS)
AC_CXX_LDFLAGS_STD_LANG(AM_LDFLAGS)
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
