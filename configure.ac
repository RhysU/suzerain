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
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.9 nostdinc no-define tar-ustar -Wall -Werror foreign)
AC_PREREQ([2.61])
AX_AM_MACROS
AX_ADD_AM_CTAGSFLAGS

dnl ------------------------------------------------
dnl LaTeX capabilities (not required but suggested)
dnl ------------------------------------------------
build_writeups=true
ACLTX_PROG_TEX(      [build_writeups=false])
ACLTX_PROG_PDFTEX(   [build_writeups=false])
ACLTX_PROG_LATEX(    [build_writeups=false])
ACLTX_PROG_PDFLATEX( [build_writeups=false])
ACLTX_PROG_BIBTEX(   [build_writeups=false])
ACLTX_PROG_MAKEINDEX([build_writeups=false])
if [ test x$build_writeups = xtrue ]; then
    AX_ADD_AM_LATEX
    ACLTX_CLASSES([amsart],AMSART_CLASS,,[build_writeups=false])
else
    AC_MSG_WARN([LaTeX prerequisites not satisfied so writeups will not be built.])
fi
AM_CONDITIONAL([BUILD_WRITEUPS], [test x$build_writeups = xtrue])
AC_CACHE_SAVE

dnl ------------------------------------------------------
dnl Check for compiler toolchain and optimization level
dnl Ordering very specific: AX_CC_MAXOPT usage is fragile
dnl ------------------------------------------------------
AX_CHECK_ENABLE_DEBUG([info])
AC_LANG_PUSH([C])
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
dnl AX_CHECK_COMPILE_FLAG([-ipo], [
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
AC_LANG_POP([C])
AC_LANG_PUSH([C++])
AC_PROG_CXX
AX_COMPILER_VENDOR
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AC_PROG_SED
AX_AM_JOBSERVER([no])
AC_LANG_POP([C++])
AC_CACHE_SAVE

dnl ----------------------------------------------
dnl Optionally enable code coverage tool
dnl ----------------------------------------------
AX_CODE_COVERAGE

dnl ----------------------------------------------
dnl Allow user to disable the tests/ subdirectory
dnl ----------------------------------------------
AC_ARG_ENABLE([unittests],
    [AS_HELP_STRING([--disable-unittests],
                    [Do not build tests/ subdirectory])],
    [enable_unittests=no], [enable_unittests=yes])
if test x$enable_unittests = xno ; then
    AC_MSG_WARN([Unit tests under tests/ subdirectory will not be built])
fi
AM_CONDITIONAL(BUILD_TESTS, test x$enable_unittests = xyes)

dnl ----------------------------------------------
dnl Enable rudimentary precompiled header support
dnl Currently very, very Intel-specific
dnl ----------------------------------------------
dnl AC_LANG_PUSH([C])
dnl AX_APPEND_COMPILE_FLAGS([-pch],[CFLAGS])
dnl AC_LANG_POP([C])
dnl AC_LANG_PUSH([C++])
dnl AX_APPEND_COMPILE_FLAGS([-pch],[CXXFLAGS])
dnl AC_LANG_POP([C++])

dnl ------------------------------------------------
dnl Enable GNU libtool for all of configure/build
dnl Current version checks also in Makefile.am
dnl ------------------------------------------------
LT_INIT
LT_OUTPUT
LT_PREREQ([2.2.6])
AC_SUBST(LIBTOOL_DEPS)
m4_append([AC_LANG(C)],[
    ac_compile="./libtool --mode=compile --tag=CC $ac_compile"
    ac_link="./libtool --mode=link --tag=CC $ac_link"
])[]dnl
m4_append([AC_LANG(C++)],[
    ac_compile="./libtool --mode=compile --tag=CXX $ac_compile"
    ac_link="./libtool --mode=link --tag=CXX $ac_link"
])[]dnl
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
export GSL_CBLAS_LIB="${MKL_LDFLAGS} ${MKL_LIBS} ${FLIBS}"
AX_PATH_GSL(1.14,,AC_MSG_ERROR([Could not find required GSL version.]))
AX_FFTW3()
PKG_CHECK_MODULES([P3DFFT],[suzerain-p3dfft])
PKG_CHECK_MODULES([GRVY], [grvy >= 0.27.1],[
        AC_DEFINE([HAVE_GRVY],[1],[Defined if GRVY Toolkit found])
    ], [
        AC_MSG_WARN([Could not find GRVY pkg-config file; Continuing...])
    ])
PKG_CHECK_MODULES([ESIO],[esio-hdf5 >= 0.1.4])
AC_LANG_POP([C])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C++ library requirements
dnl ------------------------------------------------
AC_LANG_PUSH([C++])
PKG_CHECK_MODULES([LOG4CXX],
    [liblog4cxx >= 0.10 apr-1 >= 1.3.3 apr-util-1 >= 1.3.4],
    AC_DEFINE([HAVE_LOG4CXX],1,[Define if log4cxx is available]))
AC_CACHE_CHECK([if we can link APR using LOG4CXX_LIBS], [my_cv_apr_linkable],[
    LIBS=$LOG4CXX_LIBS
    AC_TRY_LINK_FUNC([apr_initialize],
                     [my_cv_apr_linkable=yes], [my_cv_apr_linkable=no])
    LIBS=$saved_LIBS
])
if [ test x$my_cv_apr_linkable = xno ]; then
    AC_MSG_ERROR([Unable to link APR using LOG4CXX_LIBS=$LOG4CXX_LIBS])
fi
PKG_CHECK_MODULES([EIGEN],
    [eigen3 >= 3.0.0],
    AC_DEFINE([HAVE_EIGEN],1,[Define if Eigen is available]))
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
BOOST_STRING_ALGO
BOOST_TEST
BOOST_TYPETRAITS
AC_LANG_POP([C++])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C/C++ compilation strictness and compatibility
dnl ------------------------------------------------
dnl Strictness must come after libraries to avoid breaking AC_CHECK_LIB
AC_LANG_PUSH([C])
AX_CFLAGS_WARN_ALL(CFLAGS)
AX_CFLAGS_WARN_ALL(CXXFLAGS)
AX_WARNINGS_SANITIZE
AX_INTEL_SC_ANALYSIS
AC_LANG_POP([C])
AC_LANG_PUSH([C++])
AC_CXX_HAVE_STL
AC_CXX_HEADER_STDCXX_98
AC_CXX_HEADER_STDCXX_TR1
AC_CXX_COMPILE_STDCXX_0X
AC_CXX_CPPFLAGS_STD_LANG(AM_CPPFLAGS)
AC_CXX_CXXFLAGS_STD_LANG(AM_CXXFLAGS)
AC_CXX_LDFLAGS_STD_LANG(AM_LDFLAGS)
AX_CXX_EXTERN_TEMPLATE
AC_LANG_POP([C++])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Display optimization details
dnl ------------------------------------------------
if test "x$ax_enable_debug" != "xno"; then
    AX_VECTORIZATION_REPORT()
fi

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
    fields/Makefile
    lib/Makefile
    suzerain/Makefile
    tests/Makefile
    writeups/Makefile
])

AC_OUTPUT()
