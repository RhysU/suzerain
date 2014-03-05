dnl ------------------------------------------------
dnl Package overview
dnl ------------------------------------------------
AC_INIT([suzerain],[0.1.2],[rhys@ices.utexas.edu],[suzerain],[https://red.ices.utexas.edu/projects/suzerain])
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
AC_ARG_PROGRAM
AX_ENABLE_BUILDDIR
AM_INIT_AUTOMAKE(1.11 nostdinc no-define tar-ustar color-tests silent-rules parallel-tests -Wall foreign check-news subdir-objects)
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])
AC_PREREQ([2.68])
AX_CHECK_ENABLE_DEBUG([info])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AX_AM_MACROS
AX_AM_JOBSERVER([no])
AX_ADD_AM_CTAGSFLAGS
AC_PROG_AWK
AC_PROG_GREP
AC_PROG_SED
AC_CACHE_SAVE

dnl ----------------------------------------------
dnl Allow user to disable the writeups/ subdirectory
dnl ----------------------------------------------
AC_ARG_ENABLE([writeups],
    [AS_HELP_STRING([--disable-writeups],
                    [Do not build writeups/ subdirectory])],
    [build_writeups=false], [build_writeups=true])

dnl ------------------------------------------------
dnl LaTeX capabilities (not required but suggested)
dnl ------------------------------------------------
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
    AC_MSG_WARN([LaTeX prerequisites not satisfied so writeups cannot be built.])
fi
AM_CONDITIONAL([BUILD_WRITEUPS], [test x$build_writeups = xtrue])
AC_CACHE_SAVE

dnl ------------------------------------------------------
dnl Check for compiler toolchain and optimization level
dnl Ordering very specific: AX_CC_MAXOPT usage is fragile
dnl ------------------------------------------------------
AC_LANG_PUSH([C])
AC_PROG_CC
AM_PROG_CC_C_O
AM_PROG_AR
AC_OPENMP
AX_PROG_CC_C99_CFLAGS(,AC_MSG_ERROR([Compiler cannot support ISO C99]))
AX_COMPILER_VENDOR
dnl TODO get AX_APPEND_LINK_FLAGS([-Wl,--as-needed],[LDFLAGS]) working
case $ax_cv_c_compiler_vendor in dnl Avoid overlinking when possible
  *)   :                         dnl Using --as-needed with MKL breaks
       ;;
esac
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
dnl dnl Enable IPO if the (Intel) compiler supports it, assumes that
dnl dnl C/C++/Fortran all support IPO if the compiler supports it.
dnl dnl Assumption required because we must replace linker ar with xiar.
dnl AC_LANG_PUSH([C])
dnl AX_CHECK_COMPILE_FLAG([-ipo], [
dnl     AR="xiar"
dnl     CFLAGS="$CFLAGS -ipo"
dnl     CXXFLAGS="$CXXFLAGS -ipo"
dnl     FFLAGS="$FFLAGS -ipo"
dnl     FCFLAGS="$FCFLAGS -ipo"
dnl ])
fi
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AC_LANG_POP([C])
AC_LANG_PUSH([C++])
AC_PROG_CXX
AC_OPENMP
AX_COMPILER_VENDOR
AX_APPEND_COMPILE_FLAGS([${CFLAGS_MAXOPT}]) dnl Check if MAXOPT good on C++
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AC_LANG_POP([C++])
AC_LANG_PUSH([Fortran])
AC_PROG_FC
AC_OPENMP
AX_COMPILER_VENDOR
dnl ------------------------------------------------------------------------
dnl Ensure any optimization flags we maintain produce no warnings
dnl Required because Intel ifort barfs on -malign-double (see Ticket #2829)
ac_[]_AC_LANG_ABBREV[]_werror_flag=yes;      dnl Push AC_LANG_WERROR
AX_APPEND_COMPILE_FLAGS([${CFLAGS_MAXOPT}])  dnl Run tests
ac_[]_AC_LANG_ABBREV[]_werror_flag=;         dnl Pop AC_LANG_WERROR
dnl ------------------------------------------------------------------------
AC_FC_LIBRARY_LDFLAGS
AC_FC_FREEFORM
AC_FC_WRAPPERS
dnl Autoconf 2.69 got markedly better at Fortran handling
m4_version_prereq([2.69],[
    AC_FC_MODULE_EXTENSION
    AC_FC_PP_SRCEXT([F90])
],[
    AC_MSG_WARN([Configure generated older autoconf version. Please use 2.69 or later.])

    m4_pattern_allow([AC_FC_MODULE_EXTENSION])
    AC_MSG_NOTICE([BEGIN: Stubbed version of AC_FC_MODULE_EXTENSION])
    AC_SUBST([FC_MODEXT],[mod])
    AC_MSG_NOTICE([END:   Stubbed version of AC_FC_MODULE_EXTENSION])

    m4_pattern_allow([AC_FC_PP_SRCEXT])
    AC_MSG_NOTICE([BEGIN: Stubbed version of AC_FC_PP_SRCEXT])
    AC_FC_SRCEXT([F90])
    ax_save_FCFLAGS=$FCFLAGS
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
    FCFLAGS_F90="$FCFLAGS_F90 $fppflag"
    AC_MSG_NOTICE([END:   Stubbed version of AC_FC_PP_SRCEXT])
])
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AC_LANG_POP([Fortran])
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
m4_append([AC_LANG(Fortran)],[
    ac_compile="./libtool --mode=compile --tag=FC $ac_compile"
    ac_link="./libtool --mode=link --tag=FC $ac_link"
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
AC_LANG_PUSH([Fortran])
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for Fortran, e.g. mpf90.))
AC_LANG_POP([Fortran])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl ESIO requirements and associated HDF5 utilities
dnl ------------------------------------------------
AC_LANG_PUSH([C])
esio_version="esio-hdf5 >= 0.1.6"
PKG_CHECK_MODULES([ESIO], [$esio_version])
AC_MSG_CHECKING([h5pcc path used to build ESIO])
PKG_CHECK_EXISTS(
    [$esio_version],
    [pkg_cv_H5PCC=`$PKG_CONFIG --variable=h5pcc "$esio_version" 2>/dev/null`],
    [pkg_cv_H5PCC=h5pcc])
AC_MSG_RESULT($pkg_cv_H5PCC)
AC_CACHE_VAL([pkg_cv_H5PCC_path],[pkg_cv_H5PCC_path=`AS_DIRNAME(["$pkg_cv_H5PCC"])`])
AC_PATH_PROG(H5LS,   [h5ls],   [h5ls],   [$PATH$PATH_SEPARATOR$pkg_cv_H5PCC_path])
AC_PATH_PROG(H5DUMP, [h5dump], [h5dump], [$PATH$PATH_SEPARATOR$pkg_cv_H5PCC_path])
AC_CACHE_CHECK([for h5diff supporting --exclude-path], [ac_cv_path_H5DIFF], [
    AC_PATH_PROGS_FEATURE_CHECK([H5DIFF], [h5diff],
        [
         echo "Looking for --exclude-path support in $ac_path_H5DIFF" >&AS_MESSAGE_LOG_FD ; \
         $ac_path_H5DIFF --help | $GREP -q -- --exclude-path && ac_cv_path_H5DIFF=$ac_path_H5DIFF ac_cv_path_H5DIFF_found=:
        ],
        [AC_MSG_WARN([could not h5diff supporting --exclude-path])
         ac_cv_path_H5DIFF=h5diff],
        [$PATH$PATH_SEPARATOR$pkg_cv_H5PCC_path])
])
AC_SUBST([H5DIFF], [$ac_cv_path_H5DIFF])
AC_LANG_POP([C])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C library requirements
dnl ------------------------------------------------
AC_LANG_PUSH([C])
ACX_MKL()
export GSL_CBLAS_LIB="${MKL_LDFLAGS} ${MKL_LIBS}"
AX_PATH_GSL(1.15,,AC_MSG_ERROR([Could not find required GSL version.]))
AX_FFTW3()
PKG_CHECK_MODULES([GRVY], [grvy >= 0.27.1],[
        AC_DEFINE([HAVE_GRVY],[1],[Defined if GRVY Toolkit found])
    ], [
        AC_MSG_WARN([Could not find GRVY pkg-config file; Continuing...])
    ])
PKG_CHECK_MODULES([UNDERLING],[underling-fftw >= 0.3.0],[
        AC_DEFINE([HAVE_UNDERLING],[1],
                  [Defined if underling library found])
        AM_CONDITIONAL([HAVE_UNDERLING], [true])
    ], [
        AM_CONDITIONAL([HAVE_UNDERLING], [false])
        AC_MSG_NOTICE([Did not find the underling library])
    ])
AC_LANG_POP([C])
AC_CACHE_SAVE

dnl -------------------------------------------------------------------------
dnl Our on-off version of P3DFFT has been internalized as of revision r38519.
dnl For this reason, we are always certain to have a pencil library at build.
dnl Therefore the following two checks automatically succeed.  They are left
dnl here in anticipation of the day when P3DFFT can be stripped out entirely.
dnl When that happens, the underling checks above become a firm requirement.
dnl -------------------------------------------------------------------------
AM_CONDITIONAL([HAVE_P3DFFT], [true])
AC_DEFINE([HAVE_P3DFFT],[1],
          [Defined if Suzerain-specific version of P3DFFT found])

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
BOOST_REQUIRE([1.47])
BOOST_ARRAY
BOOST_BIND
BOOST_FUNCTION
BOOST_IOSTREAMS
BOOST_MATH
BOOST_MULTIARRAY
BOOST_NUMERIC_CONVERSION
BOOST_POINTER_CONTAINER
BOOST_PREPROCESSOR
BOOST_PROGRAM_OPTIONS
BOOST_REF
BOOST_SMART_PTR
BOOST_STATICASSERT
BOOST_STRING_ALGO
BOOST_TEST
BOOST_TYPETRAITS
BOOST_UTILITY
AC_LANG_POP([C++])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Fortran library requirements
dnl ------------------------------------------------
dnl None

dnl ----------------------------------------------
dnl Chemistry library requirements
dnl ----------------------------------------------
build_reacting=false
PKG_CHECK_MODULES([ANTIOCH],
    [antioch >= 0.0.0], [
        AC_DEFINE([HAVE_ANTIOCH],1,[Define if Antioch is available])
        build_reacting=true
    ],[
       : dnl Do nothing as not finding Antioch is not an error
    ])
AM_CONDITIONAL([BUILD_REACTING], [test x$build_reacting != xfalse])
if test x$build_reacting = xfalse; then
    AC_MSG_NOTICE([Because no chemistry library was detected, some functionality may not be built.])
fi
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Compilation strictness and compatibility
dnl ------------------------------------------------
dnl Strictness must come after libraries to avoid breaking AC_CHECK_LIB
AC_LANG_PUSH([C])
AX_CFLAGS_WARN_ALL
case $ax_cv_c_compiler_vendor in #(
    gnu) AX_APPEND_COMPILE_FLAGS([-Wsuggest-attribute=pure    ],[CFLAGS])
         AX_APPEND_COMPILE_FLAGS([-Wsuggest-attribute=const   ],[CFLAGS])
         AX_APPEND_COMPILE_FLAGS([-Wsuggest-attribute=noreturn],[CFLAGS])
         ;;
esac
AC_ARG_VAR([CFLAGS_PRECISE],[C compiler flags specifying no unsafe math optimizations])
case $ax_cv_c_compiler_vendor in #(
    intel)    AX_CHECK_COMPILE_FLAG([-fp-model precise],
                                    [AX_APPEND_FLAG([-fp-model precise],
                                                    [CFLAGS_PRECISE])])                  ;;
    portland) AX_APPEND_COMPILE_FLAGS([-Kieee -Mnofpapprox],           [CFLAGS_PRECISE]) ;;
    gnu|*)    AX_APPEND_COMPILE_FLAGS([-fno-unsafe-math-optimizations],[CFLAGS_PRECISE]) ;;
esac
if test "x$ax_enable_debug" = "xyes"; then
    AC_MSG_NOTICE([Debug Build Detected: Adding Precision Flags to CFLAGS.])
    CFLAGS="$CFLAGS $CFLAGS_PRECISE"
fi
AC_LANG_POP([C])
AC_LANG_PUSH([C++])
AX_CXXFLAGS_WARN_ALL
AC_CXX_CPPFLAGS_STD_LANG(AM_CPPFLAGS)
AC_CXX_CXXFLAGS_STD_LANG(AM_CXXFLAGS)
AC_CXX_LDFLAGS_STD_LANG(AM_LDFLAGS)
AX_CXX_EXTERN_TEMPLATE
case $ax_cv_cxx_compiler_vendor in #(
    gnu) AX_APPEND_COMPILE_FLAGS([-Wsuggest-attribute=pure    ],[CXXFLAGS])
         AX_APPEND_COMPILE_FLAGS([-Wsuggest-attribute=const   ],[CXXFLAGS])
         AX_APPEND_COMPILE_FLAGS([-Wsuggest-attribute=noreturn],[CXXFLAGS])
         ;;
esac
AC_ARG_VAR([CXXFLAGS_PRECISE],[C++ compiler flags specifying no unsafe math optimizations])
case $ax_cv_cxx_compiler_vendor in #(
    intel)    AX_CHECK_COMPILE_FLAG([-fp-model precise],
                                    [AX_APPEND_FLAG([-fp-model precise],
                                                    [CXXFLAGS_PRECISE])])                  ;;
    portland) AX_APPEND_COMPILE_FLAGS([-Kieee -Mnofpapprox],           [CXXFLAGS_PRECISE]) ;;
    gnu|*)    AX_APPEND_COMPILE_FLAGS([-fno-unsafe-math-optimizations],[CXXFLAGS_PRECISE]) ;;
esac
if test "x$ax_enable_debug" = "xyes"; then
    AC_MSG_NOTICE([Debug Build Detected: Adding Precision Flags to CXXFLAGS.])
    CXXFLAGS="$CXXFLAGS $CXXFLAGS_PRECISE"
fi
AC_LANG_POP([C++])
AC_LANG_PUSH([Fortran])
AX_FCFLAGS_WARN_ALL
AC_LANG_POP([Fortran])
AX_WARNINGS_SANITIZE
AX_INTEL_SC_ANALYSIS
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
DX_CHI_FEATURE([OFF])
DX_CHM_FEATURE([OFF])
DX_DOT_FEATURE([ON])
DX_HTML_FEATURE([ON])
DX_MAN_FEATURE([OFF])
DX_PDF_FEATURE([OFF])
DX_PS_FEATURE([OFF])
DX_RTF_FEATURE([OFF])
DX_XML_FEATURE([OFF])
DX_INIT_DOXYGEN([suzerain],[doxygen.cfg],[docs])

dnl ------------------------------------------------
dnl Generate Makefiles
dnl ------------------------------------------------
AC_CONFIG_FILES([
    Makefile
    apps/Makefile
    apps/perfect/Makefile
    apps/reacting/Makefile
    build-aux/Makefile
    fields/Makefile
    jobscripts/Makefile
    lib/Makefile
    lib/suzerain-p3dfft/Makefile
    postproc/Makefile
    preproc/Makefile
    suzerain/Makefile
    tests/Makefile
    writeups/Makefile
])
AC_MSG_NOTICE([Configuring Largo slow growth library as a nested autoconf package])
AC_CONFIG_SUBDIRS([largo])
AC_OUTPUT()
