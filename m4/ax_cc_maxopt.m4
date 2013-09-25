dnl @synopsis AX_CC_MAXOPT
dnl @summary turn on optimization flags for the C compiler
dnl @category C
dnl
dnl Try to turn on "good" C optimization flags for various compilers
dnl and architectures, for some definition of "good".  (In our case,
dnl good for FFTW and hopefully for other scientific codes.  Modify
dnl as needed.)
dnl
dnl The user can override the flags by setting the CFLAGS environment
dnl variable.
dnl
dnl Note also that the flags assume that ANSI C aliasing rules are
dnl followed by the code (e.g. for gcc's -fstrict-aliasing), and that
dnl floating-point computations can be re-ordered as needed.
dnl
dnl Modified by Rhys Ulerich for Suzerain's requirements.  Specifically,
dnl Intel workarounds and no requirements to preserve the ABI.
dnl
dnl Requires macros: AX_CHECK_COMPILE_FLAG, AX_COMPILER_VENDOR,
dnl
dnl @version 2012-01-05
dnl @license GPLWithACException
dnl @author Steven G. Johnson <stevenj@alum.mit.edu> and Matteo Frigo.
AC_DEFUN([AX_CC_MAXOPT],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AX_COMPILER_VENDOR])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine "good" native compiler flags if none specified via CFLAGS
if test "$ac_test_CFLAGS" != "set"; then
  CFLAGS=""
  case $ax_cv_c_compiler_vendor in
    dec) CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
         ;;

    sun) CFLAGS="-native -fast -xO5 -dalign"
         ;;

    hp)  CFLAGS="+Oall +Optrs_ansi +DSnative"
         ;;

    ibm) xlc_opt="-qtune=auto"
         AX_CHECK_COMPILE_FLAG($xlc_opt,
                CFLAGS="-O3 -qansialias -w $xlc_opt",
               [CFLAGS="-O3 -qansialias -w"
                echo "******************************************************"
                echo "*  You seem to have the IBM  C compiler.  It is      *"
                echo "*  recommended for best performance that you use:    *"
                echo "*                                                    *"
                echo "*    CFLAGS=-O3 -qarch=xxx -qtune=xxx -qansialias -w *"
                echo "*                      ^^^        ^^^                *"
                echo "*  where xxx is pwr2, pwr3, 604, or whatever kind of *"
                echo "*  CPU you have.  (Set the CFLAGS environment var.   *"
                echo "*  and re-run configure.)  For more info, man cc.    *"
                echo "******************************************************"])
         ;;

    intel) CFLAGS="-O3"

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
        [[#if __INTEL_COMPILER >= 1110 && __INTEL_COMPILER < 1200
         ah ha: icc is exactly version 11.1
         introducing a version 11.1 workaround from Intel support #602718
        #endif
        ]])], [], [CFLAGS="$CFLAGS -mP2OPT_cndxform_max_new=-1"])

        # Intel seems to have changed the spelling of this flag recently
        # ANSI aliasing is problematic after Intel 11.1 and enabling it
        # seems to make trivial performance differences according to Karl's
        # testing within Redmine #2729.  Explicitly turning on -no-ansi-alias
        # to document this decision and prevent new Intel -O3 hidden behavior
        # from turning it on again.
        icc_no_ansi_alias="unknown"
        for flag in -no-ansi-alias -no-ansi_alias; do
          AX_CHECK_COMPILE_FLAG($flag, [icc_no_ansi_alias=$flag; break])
        done
        if test "x$icc_no_ansi_alias" != xunknown; then
            CFLAGS="$CFLAGS $icc_no_ansi_alias"
        fi
        AX_CHECK_COMPILE_FLAG(-malign-double, CFLAGS="$CFLAGS -malign-double")

        # Check for architecture flags here, e.g. -xHost etc.  Beware of poor
        # interaction between -mavx and -xHost per FFTW ax_cc_maxopt and issues
        # with the resulting ABI changing.  Note that just checking for -xHost
        # is a no-go as pre-11.1 versions ignore it in a not-so-silent fashion.
        icc_archflag=unknown
        icc_flags=""
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
        [[#if __INTEL_COMPILER >= 1110
          ah ha: icc is at least version 11.1
        #endif
        ]])], [], [icc_flags=-xHost])
        case $host_cpu in
          i686*|x86_64*)
            # icc accepts gcc assembly syntax, so these should work:
            AX_GCC_X86_CPUID(0)
            AX_GCC_X86_CPUID(1)
            case $ax_cv_gcc_x86_cpuid_0 in # TODO see AX_GCC_ARCHFLAG
              *:756e6547:*:*) # Intel
                case $ax_cv_gcc_x86_cpuid_1 in
                  *6a?:*[[234]]:*:*|*6[[789b]]?:*:*:*) icc_flags="-xK";;
                  *f3[[347]]:*:*:*|*f4[1347]:*:*:*) icc_flags="-xP -xN -xW -xK";;
                  *f??:*:*:*) icc_flags="-xN -xW -xK";;
                esac ;;
            esac ;;
        esac
        if test "x$icc_flags" != x; then
          for flag in $icc_flags; do
            AX_CHECK_COMPILE_FLAG($flag, [icc_archflag=$flag; break])
          done
        fi
        AC_MSG_CHECKING([for icc architecture flag])
        AC_MSG_RESULT($icc_archflag)
        if test "x$icc_archflag" != xunknown; then
          CFLAGS="$CFLAGS $icc_archflag"
        fi
        ;;

    gnu)
        # Default optimization flags for gcc on all systems.
        CFLAGS="-O3"

        # -O3 does not imply -fomit-frame-pointer on ia32
        AX_APPEND_COMPILE_FLAGS([-fomit-frame-pointer])

        # tune for the host by default
        AX_APPEND_COMPILE_FLAGS([-mtune=native])

        # -malign-double for x86 systems
        AX_APPEND_COMPILE_FLAGS([-malign-double])

        # -fstrict-aliasing for gcc-2.95+
        AX_APPEND_COMPILE_FLAGS([-fstrict-aliasing])

        # We enable "unsafe" fp optimization with other compilers, too
        # We do rely on some IEEE promises hence no -ffast-math
        AX_APPEND_COMPILE_FLAGS([-funsafe-math-optimizations])
        ;;
  esac

  if test -z "$CFLAGS"; then
        echo ""
        echo "********************************************************"
        echo "* WARNING: Don't know the best CFLAGS for this system  *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
        echo "* (otherwise, a default of CFLAGS=-O3 will be used)    *"
        echo "********************************************************"
        echo ""
        CFLAGS="-O3"
  fi

  AX_CHECK_COMPILE_FLAG($CFLAGS, [], [
        echo ""
        echo "********************************************************"
        echo "* WARNING: The guessed CFLAGS don't seem to work with  *"
        echo "* your compiler.                                       *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
        echo "********************************************************"
        echo ""
        CFLAGS=""
  ])

fi
])
