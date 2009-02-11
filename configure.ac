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
AX_PATH_P3DFFT(,AC_MSG_ERROR([Could not find the P3DFFT library.]))
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C++ library requirements
dnl ------------------------------------------------
AC_LANG([C++])
PKG_CHECK_MODULES(LOG4CXX, liblog4cxx >= 0.10)
AC_SUBST(LOG4CXX_CFLAGS)
AC_SUBST(LOG4CXX_LIBS)

AX_BOOST_BASE(1.34)
AX_BOOST_UNIT_TEST_FRAMEWORK
AX_BOOST_PROGRAM_OPTIONS
## Potential future expansion
# AX_TRILINOS_BASE
# AX_TRILINOS_TEUCHOS
# AX_TRILINOS_EPETRA
# AX_TRILINOS_EPETRAEXT
# AX_TRILINOS_THYRA
# AX_TRILINOS_THYRA_EPETRA
# AX_TRILINOS_THYRA_EPETRAEXT
# AX_TRILINOS_RTOP
# AX_TRILINOS_RYTHMOS
## Potential future expansion
# AX_TRILINOS_AMESOS
# AX_TRILINOS_EPETRAEXT_HDF5
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Checks for the literate programming toolchain
dnl ------------------------------------------------
## Potential future expansion
# AX_ADD_AM_CWEB_LATEX
# ACLTX_CLASS_ARTICLE(
#     AC_MSG_ERROR([LaTeX article class not found.]))
# ACLTX_PACKAGE(rcs,article,rcs_latex,,
#     AC_MSG_ERROR([LaTeX package rcs not functional.]))

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
AX_CFLAGS_WARN_ALL(CFLAGS)
AX_CFLAGS_WARN_ALL(CXXFLAGS)
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Prepare log4cxx properties in the build tree
dnl ------------------------------------------------
AC_CONFIG_LINKS([src/log4cxx.properties:src/log4cxx.properties.in])

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
