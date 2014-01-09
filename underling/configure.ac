dnl ----------------------------------------------------
dnl Package name, release number, description, and URL
dnl ----------------------------------------------------
AC_INIT([underling], [0.3.1], [rhys@ices.utexas.edu])
AC_REVISION([$Id$])
PACKAGE_DESCRIPTION="an FFTW MPI-based library for 3D pencil decompositions"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="http://red.ices.utexas.edu/projects/underling"
AC_SUBST([PACKAGE_URL])

dnl ------------------------------------------------
dnl Initialization macros
dnl ------------------------------------------------
AC_CONFIG_COMMANDS_PRE([m4_ifdef([AH_HEADER],
    [AC_SUBST([CONFIG_INCLUDE], m4_defn([AH_HEADER]))])])
AC_CONFIG_SRCDIR([underling/underling.h])
AM_CONFIG_HEADER(underling/config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AX_PREFIX_CONFIG_H([underling/underling-config.h])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.11 no-define tar-ustar color-tests -Wall foreign check-news)
AC_PREREQ([2.64])
AX_AM_MACROS
AX_ADD_AM_CTAGSFLAGS
AX_AM_JOBSERVER([no])
AX_SPLIT_VERSION
AC_SUBST([AX_MAJOR_VERSION])
AC_SUBST([AX_MINOR_VERSION])
AC_SUBST([AX_POINT_VERSION])
AC_PATH_PROG([PERL],[perl],[no])
AM_CONDITIONAL([HAVE_PERL],[test x$PERL != xno])
AC_CACHE_SAVE

dnl ------------------------------------------------------
dnl Check for MPI-ready C/C++ compiler toolchain
dnl ------------------------------------------------------
AX_CHECK_ENABLE_DEBUG([info])
AC_LANG_PUSH([C])
AC_PROG_CC
gl_EARLY
gl_INIT
AX_PROG_CC_C99_CFLAGS(,AC_MSG_ERROR([Compiler cannot support ISO C99]))
AX_CFLAGS_WARN_ALL
AX_COMPILER_VENDOR()
AX_VISIBILITY([hidden])
AX_TLS()
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for C, e.g. mpicc.))
MPICC_LOCATION=$(eval command -v $MPICC)
AC_SUBST([MPICC_LOCATION])
AC_LANG_POP([C])
AC_CACHE_SAVE
AC_LANG_PUSH([C++])
AC_PROG_CXX
AX_COMPILER_VENDOR()
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for C++, e.g. mpicxx.))
MPICXX_LOCATION=$(eval command -v $MPICXX)
AC_SUBST([MPICXX_LOCATION])
AC_LANG_POP([C++])
AX_WARNINGS_SANITIZE()
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Enable GNU libtool
dnl Current version checks also in Makefile.am
dnl ------------------------------------------------
LT_INIT([])
LT_PREREQ([2.2.6])
AC_SUBST(LIBTOOL_DEPS)
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Check for non-MPI dependencies for C/C++
dnl ------------------------------------------------
AC_LANG_PUSH([C])
AX_PATH_GRVY_NEW([0.27],[no])
AX_MPIP()
AX_FFTW3()
if test "x$ax_fftw3_found_mpi" != "xyes"; then
    AC_MSG_ERROR([FFTW3 installation with MPI support not detected.])
fi
saved_CC=$CC
CC=$MPICC
AC_CHECK_LIB([fftw3_mpi],[fftw_mpi_execute_dft],[
        AC_DEFINE([HAVE_FFTW33_BETA1_OR_LATER],
                  [1], [Are we building against 3.3-beta1 or later?])
    ],[
       AC_MSG_WARN([FFTW 3.3-beta1 or later was not detected and is strongly recommended!])
    ],[$FFTW3_MPI_THREADS_LIBS])
CC=$saved_CC
AC_LANG_POP([C])
AC_CACHE_SAVE
dnl Boost is only required to build test cases
AC_LANG_PUSH([C++])
have_boost="yes"
BOOST_REQUIRE([1.39],[have_boost="no"])
if test "$have_boost" = "yes"; then
    BOOST_TEST
else
    AC_MSG_WARN([Boost.Test not usable so unit tests will not be built])
fi
AM_CONDITIONAL([BUILD_TESTS], [test "$have_boost" = yes])
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
DX_INIT_DOXYGEN([underling],[doxygen/doxygen.cfg],[docs])

#--------------------------
# Checks for code coverage
#-------------------------
AX_CODE_COVERAGE

dnl ------------------------------------------------
dnl Generate Makefiles
dnl ------------------------------------------------
AC_SUBST([CONFIG_STATUS_DEPENDENCIES], ['$(top_srcdir)/LICENSE.in'])
AC_CONFIG_FILES([
    LICENSE
    Makefile
    underling.pc
    underling-fftw.pc
    underling-fftw-pthreads.pc
    underling-fftw-openmp.pc
    apps/Makefile
    doxygen/Makefile
    underling/Makefile
    underling/version.h
    examples/Makefile
    gnulib/Makefile
    tests/Makefile
])

AC_OUTPUT()
