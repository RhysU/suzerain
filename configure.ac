dnl ------------------------------------------------
dnl Package name and version number
dnl ------------------------------------------------
AC_INIT(suzerain, [0.0.1], [rhys@ices.utexas.edu])
AC_REVISION([$Id$])

dnl ------------------------------------------------
dnl Initialization macros
dnl ------------------------------------------------
AC_CONFIG_SRCDIR([src/helloworld.cc])
AM_CONFIG_HEADER(src/config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(1.9 no-define tar-ustar -Wall -Werror foreign)
AC_PREREQ([2.61])
AX_AM_MACROS
AX_ADD_AM_LATEX

dnl ------------------------------------------------
dnl Check for compiler toolchain
dnl ------------------------------------------------
AX_CHECK_ENABLE_DEBUG
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AC_PROG_SED
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Enable GNU libtool
dnl Current version checks also in Makefile.am
dnl ------------------------------------------------
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

dnl TODO Automatically add optimization if debugging disabled per AX_CC_MAXOPT
dnl Issues: 1) AX_CC_MAXOPT doesn't play nicely with AX_CHECK_ENABLE_DEBUG
dnl         2) AX_CC_MAXOPT only sets CFLAGS, not CXXFLAGS

dnl ------------------------------------------------
dnl C library requirements
dnl ------------------------------------------------
AC_LANG([C])
ACX_BLAS([
    dnl Workaround for bogus FLIBS
    FLIBS=`echo ${FLIBS} | ${SED} 's/-lgfortranbegin//'`
    export GSL_CBLAS_LIB="${BLAS_LIBS} ${FLIBS}"
    AC_SUBST(GSL_CBLAS_LIB)
],AC_MSG_WARN([Will use gslcblas if it is present.  Try --with-blas=<lib>.]))
AX_PATH_GSL(1.12,,AC_MSG_ERROR([Could not find required GSL version.]))
PKG_CHECK_MODULES(FFTW3, fftw3 >= 3.2)
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)
AX_PATH_P3DFFT(,AC_MSG_ERROR([Could not find the P3DFFT library.]))
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C++ library requirements
dnl ------------------------------------------------
AC_LANG([C++])
PKG_CHECK_MODULES(LOG4CXX, liblog4cxx >= 0.10 apr-1 >= 1.3.3 apr-util-1 >= 1.3.4)
AC_SUBST(LOG4CXX_CFLAGS)
AC_SUBST(LOG4CXX_LIBS)
BOOST_REQUIRE([1.37])
BOOST_FOREACH
BOOST_FORMAT
BOOST_PROGRAM_OPTIONS
BOOST_TEST
dnl TODO Check for boost::numeric::ublas
AC_DEFINE(BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR,1,[Enable boost::numeric::ublas::shallow_array_adaptor<T>])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C++ compilation strictness and compatibility
dnl ------------------------------------------------
dnl Strictness must come after libraries to avoid breaking AC_CHECK_LIB
AC_CXX_BOOL
AC_CXX_HAVE_STL
AC_CXX_HEADER_STDCXX_98
AC_CXX_HEADER_STDCXX_TR1
AC_CXX_COMPILE_STDCXX_0X
AC_CXX_CPPFLAGS_STD_LANG(AM_CPPFLAGS)
AC_CXX_CXXFLAGS_STD_LANG(AM_CXXFLAGS)
AC_CXX_LDFLAGS_STD_LANG(AM_LDFLAGS)
dnl TODO Find a nice way to enable/disable -Wall at configure time
dnl AX_CFLAGS_WARN_ALL_ANSI(CFLAGS)
dnl AX_CFLAGS_WARN_ALL_ANSI(CXXFLAGS)
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl MPI toolchain checks
dnl ------------------------------------------------
AC_LANG([C])
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for C, e.g. mpicc.))
AC_LANG([C++])
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for C++, e.g. mpicxx.))
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Generate Makefiles
dnl ------------------------------------------------
AC_CONFIG_FILES([
    Makefile
    doc/Makefile
    lib/Makefile
    src/Makefile
])

AC_OUTPUT()
