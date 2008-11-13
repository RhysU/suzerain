# Configure path for GSLwrap, a C++ wrapper for the GNU Scientific Library
# Christopher R. Gabriel <cgabriel@linux.it>, April 2000
# Modified by Rhys Ulerich <rhys.ulerich@gmail.com>, November 2008

AC_DEFUN([AX_PATH_GSLWRAP],
[
AC_ARG_WITH(gslwrap-prefix,[  --with-gslwrap-prefix=PFX   Prefix where GSLwrap is installed (optional)],
            gslwrap_prefix="$withval", gslwrap_prefix="")
AC_ARG_WITH(gslwrap-exec-prefix,[  --with-gslwrap-exec-prefix=PFX Exec prefix where GSLwrap is installed (optional)],
            gslwrap_exec_prefix="$withval", gslwrap_exec_prefix="")
AC_ARG_ENABLE(gslwraptest, [  --disable-gslwraptest       Do not try to compile and run a test GSLwrap program],
		    , enable_gslwraptest=yes)

  if test "x${GSLWRAP_CONFIG+set}" != xset ; then
     if test "x$gslwrap_prefix" != x ; then
         GSLWRAP_CONFIG="$gslwrap_prefix/bin/gslwrap-config"
     fi
     if test "x$gslwrap_exec_prefix" != x ; then
        GSL_CONFIG="$gslwrap_exec_prefix/bin/gslwrap-config"
     fi
  fi

  AC_PATH_PROG(GSLWRAP_CONFIG, gslwrap-config, no)
  min_gslwrap_version=ifelse([$1], ,0.2.0,$1)
  AC_MSG_CHECKING(for GSLwrap - version >= $min_gslwrap_version)
  no_gslwrap=""
  if test "$GSLWRAP_CONFIG" = "no" ; then
    no_gslwrap=yes
  else
    GSLWRAP_CFLAGS=`$GSLWRAP_CONFIG --cflags`
    GSLWRAP_LIBS=`$GSLWRAP_CONFIG --libs`

    gslwrap_major_version=`$GSLWRAP_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${gslwrap_major_version}" = "x" ; then
       gslwrap_major_version=0
    fi

    gslwrap_minor_version=`$GSLWRAP_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${gslwrap_minor_version}" = "x" ; then
       gslwrap_minor_version=0
    fi

    gslwrap_micro_version=`$GSLWRAP_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${gslwrap_micro_version}" = "x" ; then
       gslwrap_micro_version=0
    fi

    if test "x$enable_gslwraptest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GSLWRAP_CFLAGS"
      LIBS="$LIBS $GSLWRAP_LIBS"

      rm -f conf.gslwraptest
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* my_strdup (const char *str);

char*
my_strdup (const char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = (char *)malloc ((strlen (str) + 1) * sizeof(char));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main (void)
{
  int major = 0, minor = 0, micro = 0;
  int n;
  char *tmp_version;

  system ("touch conf.gslwraptest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_gslwrap_version");

  n = sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) ;

  if (n != 2 && n != 3) {
     printf("%s, bad version string\n", "$min_gslwrap_version");
     exit(1);
   }

   if (($gslwrap_major_version > major) ||
      (($gslwrap_major_version == major) && ($gslwrap_minor_version > minor)) ||
      (($gslwrap_major_version == major) && ($gslwrap_minor_version == minor) && ($gslwrap_micro_version >= micro)))
    {
      exit(0);
    }
  else
    {
      printf("\n*** 'gslwrap-config --version' returned %d.%d.%d, but the minimum version\n", $gslwrap_major_version, $gslwrap_minor_version, $gslwrap_micro_version);
      printf("*** of GSLwrap required is %d.%d.%d. If gslwrap-config is correct, then it is\n", major, minor, micro);
      printf("*** best to upgrade to the required version.\n");
      printf("*** If gslwrap-config was wrong, set the environment variable GSLWRAP_CONFIG\n");
      printf("*** to point to the correct copy of gslwrap-config, and remove the file\n");
      printf("*** config.cache before re-running configure\n");
      exit(1);
    }
}

],, no_gslwrap=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gslwrap" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GSLWRAP_CONFIG" = "no" ; then
       echo "*** The gslwrap-config script installed by GSLwrap could not be found"
       echo "*** If GSLwrap was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GSLWRAP_CONFIG environment variable to the"
       echo "*** full path to gslwrap-config."
     else
       if test -f conf.gslwraptest ; then
        :
       else
          echo "*** Could not run GSLwrap test program, checking why..."
          CFLAGS="$CFLAGS $GSLWRAP_CFLAGS"
          LIBS="$LIBS $GSLWRAP_LIBS"
          AC_TRY_LINK([
#include <stdio.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GSLwrap or finding the wrong"
          echo "*** version of GSLwrap. If it is not finding GSLwrap, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GSLwrap was incorrectly installed"
          echo "*** or that you have moved GSLwrap since it was installed. In the latter case, you"
          echo "*** may want to edit the gslwrap-config script: $GSLWRAP_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
#     GSLWRAP_CFLAGS=""
#     GSLWRAP_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GSLWRAP_CFLAGS)
  AC_SUBST(GSLWRAP_LIBS)
  rm -f conf.gslwraptest
])

AU_ALIAS([AM_PATH_GSLWRAP], [AX_PATH_GSLWRAP])
