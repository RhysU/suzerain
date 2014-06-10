dnl ----------------------------------------------------
dnl Package name, release number, description, and URL
dnl ----------------------------------------------------
AC_INIT([largo], [0.0.2], [topalian@ices.utexas.edu])
AC_REVISION([$Id$])
PACKAGE_DESCRIPTION="largo - slow growth terms for turbulence simulations"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="http://pecos.ices.utexas.edu/"
AC_SUBST([PACKAGE_URL])

dnl ------------------------------------------------
dnl Initialization macros
dnl ------------------------------------------------
AC_CONFIG_SRCDIR([largo/largo.h])
AM_CONFIG_HEADER(largo/config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.11 no-define tar-ustar -Wall foreign color-tests silent-rules parallel-tests check-news)
AC_PREREQ([2.68])
AX_AM_MACROS
AX_SPLIT_VERSION
AC_SUBST([AX_MAJOR_VERSION])
AC_SUBST([AX_MINOR_VERSION])
AC_SUBST([AX_POINT_VERSION])
AC_PATH_PROG([PERL],[perl],[no])
AM_CONDITIONAL([HAVE_PERL],[test x$PERL != xno])
AC_CACHE_SAVE

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

dnl ---------------------
dnl Manipulate C settings
dnl ---------------------
AX_CHECK_ENABLE_DEBUG([info])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_CC
AC_LANG_PUSH([C])
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AX_CFLAGS_WARN_ALL
AX_PROG_CC_C99_CFLAGS(,AC_MSG_ERROR([Compiler cannot support ISO C99]))
AX_COMPILER_VENDOR()
AX_WARNINGS_SANITIZE()
AC_LANG_POP([C])
AC_CACHE_SAVE

dnl ------------------
dnl Manipulate Fortran
dnl ------------------
AC_PROG_FC
AC_LANG_PUSH([Fortran])
dnl Poor man's PP_DEFINE allows Autoconf 2.68+ instead of requiring 2.69+
m4_ifdef([AC_FC_PP_DEFINE], [AC_FC_PP_DEFINE],[AC_SUBST([FC_DEFINE],[-D])])
AC_FC_SRCEXT([f90]) dnl PP_DEFINE must be before SRCEXT for some reason
AC_FC_FREEFORM()
AC_FC_LINE_LENGTH([unlimited])
dnl BEGIN: Enable the preprocessor for all Fortran .f90 sources
dnl Stupidly complicated to get correct.  Modify FCFLAGS_f90 with result.
dnl
ax_save_FCFLAGS=$FCFLAGS
AC_LANG_WERROR
AC_MSG_CHECKING([for flags required to enable Fortran preprocessing])
for fppflag in '' -fpp -cpp '-x f95-cpp-input' none
do
    if test "x$fppflag" == xnone; then
        AC_MSG_RESULT([unknown])
        AC_MSG_ERROR([Could not enable Fortran preprocessing])
        break
    fi
    ax_fpp_save_FCFLAGS=$FCFLAGS
    FCFLAGS="$FCFLAGS $fppflag"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[#define TESTDEFINE])],[
                         FCFLAGS=$ax_fpp_save_FCFLAGS
                         break
                      ],[
                         FCFLAGS=$ax_fpp_save_FCFLAGS
                      ])
done
AC_MSG_RESULT([$fppflag])
FCFLAGS_f90="$FCFLAGS_f90 $fppflag"
dnl END: Enable the preprocessor for all Fortran sources
AC_FC_WRAPPERS
AX_FCFLAGS_WARN_ALL
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AC_LANG_POP([Fortran])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Enable GNU libtool
dnl Current version checks also in Makefile.am
dnl ------------------------------------------------
LT_INIT()
LT_PREREQ([2.2.6])
AC_SUBST(LIBTOOL_DEPS)
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Documentation generation
dnl ------------------------------------------------
DX_DOXYGEN_FEATURE([ON])
DX_HTML_FEATURE([ON])
DX_CHM_FEATURE([OFF])
DX_CHI_FEATURE([OFF])
DX_MAN_FEATURE([ON])
DX_RTF_FEATURE([OFF])
DX_XML_FEATURE([OFF])
DX_PDF_FEATURE([OFF])
DX_PS_FEATURE([OFF])
DX_INIT_DOXYGEN([largo],[doxygen/doxygen.cfg],[docs])

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
    largo.pc
    doxygen/Makefile
    largo/Makefile
    largo/version.h
    tests/Makefile
    writeups/Makefile
])

AC_OUTPUT()
