/*  SvnRev
 *
 *  This utility retrieves the highest number that follows the "$Id: $" keyword
 *  or a combination of the $Rev: $ and $Date: $ keywords. The Subversion
 *  version control system expands these keywords and keeps them up to date.
 *  For an example of the tag, see the end of this comment.
 *
 *  Details on the usage and the operation of this utility is available on-line
 *  at http://www.compuphase.com/svnrev.htm.
 *
 *
 *  Acknowledgements
 *
 *  The support for .java files is contributed by Tom McCann (tommc@spoken.com).
 *  The support for C# is contributed by Gunther Zander (gzander@gmx.li).
 *  The option for prefixing and/or suffixing the build number (in the string
 *  constant SVN_REVSTR) was suggested by Robert Nitzel.
 *
 *
 *  License
 *
 *  Copyright (c) 2005-2011, ITB CompuPhase (www.compuphase.com).
 *
 *  This software is provided "as-is", without any express or implied warranty.
 *  In no event will the authors be held liable for any damages arising from
 *  the use of this software.
 *
 *  Permission is granted to anyone to use this software for any purpose,
 *  including commercial applications, and to alter it and redistribute it
 *  freely, subject to the following restrictions:
 *
 *  1.  The origin of this software must not be misrepresented; you must not
 *      claim that you wrote the original software. If you use this software in
 *      a product, an acknowledgment in the product documentation would be
 *      appreciated but is not required.
 *  2.  Altered source versions must be plainly marked as such, and must not be
 *      misrepresented as being the original software.
 *  3.  This notice may not be removed or altered from any source distribution.
 *
 * Version: $Id: svnrev.c 4499 2011-05-03 13:18:01Z thiadmer $
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "svnrev.h"


#if defined __WIN32__ || defined _Win32 || defined _WIN32 || defined WIN32
  #define DIRSEP '\\'
  #define strcasecmp  stricmp
#elif defined macintosh
  #define DIRSEP ':'
#else
  /* assume Linux/Unix */
  #define DIRSEP '/'
#endif

#define MAX_LINELENGTH      512
#define MAX_SYMBOLLENGTH    32
enum {
  LANG_C_CPP,   /* C/C++ format, a ".h" file */
  LANG_JAVA,
  LANG_CS,      /* C# format */
  LANG_PROP,
  LANG_ORACLE   /* Oracle package specification .sql */
};

static void about(void)
{
  printf("svnrev 1.12." SVN_REVSTR "\n\n");
  printf("Copyright 2005-2011 ITB CompuPhase. See http://www.compuphase.com/ for source\ncode and documentation\n\n");
  printf("Usage: svnrev [options] <input> [input [...]]\n\n"
         "Options:\n"
         "-jname\t\tJava: this option writes a java package file instead of a C/C++\n"
         "\t\theader file. The name of the Java package must follow the\n"
         "\t\toption (this is not the filename).\n\n"
         "-cname\t\tC#: this option writes a c-sharp class file instead of a C/C++\n"
         "\t\theader file. The name of the namespace and class must follow\n"
         "\t\tthe option. Namespace and class names are separated by a\n"
         "\t\tperiod (namespace.class).\n\n"
         "-Oname\t\tOracle: this option creates a Oracle package file instead of\n"
         "\t\ta C/C++ header file. The name of the package must follow the\n"
         "\t\toption.\n\n"
         "-P\t\tProperty file: This option creates a \"property file\", that can\n"
         "\t\tthen be used with ant.\n\n"
         "-pfilename\tC#: Assembly Property Filename (AssemblyInfo.cs)\n"
         "\t\tthe AssemblyVersion and AssemblyFileVersion will be updated.\n\n"
         "-ofilename\tOutput filename for the file with the build number. When no\n"
         "\t\tfilename follows \"-o\", the result is written to stdout. The\n"
         "\t\tdefault filename is \"svnrev.h\" for C/C++, \"SvnRevision.java\"\n"
         "\t\tfor Java and \"SvnRevision.cs\" for C#.\n\n"
         "-i\t\tIncremental: this option should be used when the list of input\n"
         "\t\tfiles is a subset of all files in the project. When -i is\n"
         "\t\tpresent, svnrev also scans the output file that was generated\n"
         "\t\ton a previous run.\n\n"
         "-fpattern\tFormat: adds text before or after the build number in the\n"
         "\t\tconstant SVN_REVSTR. The pattern has the form \"text#text\"\n"
         "\t\t(without the quotes) where \"text\" is arbitrary text and \"#\"\n"
         "\t\twill be replaced by the build number.\n\n"
         "-mname\t\tMacro prefix: by default all generated constants start with\n"
         "\t\tSVN_REV, but this can be changed with this option.\n\n"
         "-n\t\tIgnore line endings (CR versus LF versus CR+LF) when comparing\n"
         "\t\tthe working copy to the base.\n\n"
         "-v\t\tVerbose: prints the names of files that are modified since the\n"
         "\t\tlast commit (into version control) to stderr.\n");
  exit(1);
}

static void processfile(const char *name, int failsilent, int ignore_eol,
                        int *max_build, int *accum_build,
                        int *max_year, int *max_month, int *max_day,
                        int *ismodified)

{
  char str[MAX_LINELENGTH], str_base[MAX_LINELENGTH];
  char name_base[MAX_LINELENGTH];
  char *p1;
  FILE *fp, *fp_base;
  int build, maj_build;
  int year, month, day;
  int cnt;
  char modchar;

  /* since we also want to verify whether the file is modified in version
   * control, get the path to the working copy name
   * for every source file "<path>\<filename>, the "working copy" base can
   * be found in "<path>\.svn\text-base\<filename>.svn-base"
   */
  if ((p1 = strrchr(name, DIRSEP)) != NULL) {
    ++p1; /* skip directory separator character ('\' in Windows, '/' in Linux) */
    strncpy(name_base, name, (int)(p1 - name));
    name_base[(int)(p1 - name)] = '\0';
  } else {
    name_base[0] = '\0';
    p1 = (char*)name;
  } /* if */
  sprintf(name_base + strlen(name_base), ".svn%ctext-base%c%s.svn-base",
          DIRSEP, DIRSEP, p1);

  /* first extract the revision keywords */
  fp = fopen(name, "r");
  if (fp == NULL) {
    if (!failsilent)
      fprintf(stderr, "Failed to open input file '%s'\n", name);
    return;
  } /* if */
  fp_base = fopen(name_base, "r");  /* fail silently */
  build = 0;
  maj_build = 0;      /* RCS / CVS */
  year = month = day = 0;

  while (fgets(str, sizeof str, fp) != NULL) {
    if (fp_base == NULL || fgets(str_base, sizeof str_base, fp_base) == NULL)
      str_base[0] = '\0';
    if ((p1 = strstr(str, "$Id:")) != NULL && strchr(p1+1, '$') != NULL) {
      if (sscanf(p1, "$Id: %*s %d %d-%d-%d", &build, &year, &month, &day) < 4
          && sscanf(p1, "$Id: %*s %d %d/%d/%d", &build, &year, &month, &day) < 4)
        if (sscanf(p1, "$Id: %*s %d.%d %d-%d-%d", &maj_build, &build, &year, &month, &day) < 5)
          sscanf(p1, "$Id: %*s %d.%d %d/%d/%d", &maj_build, &build, &year, &month, &day);
    } else if ((p1 = strstr(str, "$Rev:")) != NULL && strchr(p1+1, '$') != NULL) {
      if (sscanf(p1, "$Rev: %d.%d", &maj_build, &build) < 2) {
        sscanf(p1, "$Rev: %d", &build);
        maj_build = 0;
      } /* if */
    } else if ((p1 = strstr(str, "$Revision:")) != NULL && strchr(p1+1, '$') != NULL) {
      if (sscanf(p1, "$Revision: %d.%d", &maj_build, &build) < 2) {
        /* SvnRev also writes this keyword in its own generated file; read it
         * back for partial updates
         */
        cnt = sscanf(p1, "$Revision: %d%c", &build, &modchar);
        if (cnt == 2 && modchar == 'M' && ismodified != NULL)
          *ismodified = 1;
        maj_build = 0;
      } /* if */
    } else if ((p1 = strstr(str, "$Date:")) != NULL && strchr(p1+1, '$') != NULL) {
      if (sscanf(p1, "$Date: %d-%d-%d", &year, &month, &day) < 3)
        sscanf(p1, "$Date: %d/%d/%d", &year, &month, &day);
    } else if (ismodified != NULL && *ismodified == 0 && fp_base != NULL) {
      /* no keyword present, compare the lines for equivalence,
       * but (optionally) first cut off the strings at the first CR/LF character
       * (for cross-platform compatibility)
       */
      if (ignore_eol) {
        for (p1 = strchr(str, '\0'); p1 != NULL && p1 > str && (*(p1 - 1) == '\r' || *(p1 - 1) == '\n'); p1--)
          *--p1 = '\0';
        for (p1 = strchr(str_base, '\0'); p1 != NULL && p1 > str_base && (*(p1 - 1) == '\r' || *(p1 - 1) == '\n'); p1--)
          *--p1 = '\0';
      } /* if */
      *ismodified = strcmp(str, str_base) != 0;
    } /* if */

    if (maj_build)
      *accum_build += build;            /* RCS / CVS */
    else if (build > *max_build)
      *max_build = build;               /* Subversion */
    if (year > *max_year
        || (year == *max_year && month > *max_month)
        || (year == *max_year && month == *max_month && day > *max_day))
    {
        *max_year = year;
        *max_month = month;
        *max_day = day;
    } /* if */
    if (build > 0 && year > 0 && (fp_base == NULL || ismodified == NULL || *ismodified != 0))
      break;      /* both build # and date found, not comparing or modification
                   * already found => no need to search further */

  } /* while */
  fclose(fp);
  if (fp_base != NULL)
    fclose(fp_base);
}

static int scanmodifications(const char *filename, const char *basename)
{
  FILE *fp;
  char str[MAX_LINELENGTH];
  char constant[128];
  char *ptr;
  size_t constlen;
  int modifcount = 0;

  if ((fp = fopen(filename, "r")) == NULL)
    return 0;

  sprintf(constant, "%sMODIFIED", basename);
  constlen = strlen(constant);
  while (fgets(str, sizeof str, fp) != NULL) {
    if ((ptr = strstr(str, constant)) != NULL && ptr > str && *(ptr - 1) <= ' ' && *(ptr + constlen) <= ' ')
      modifcount = strtol(ptr + constlen, NULL, 10);
  } /* while */

  fclose(fp);
  return modifcount;
}

static void get_buildno_ptr(char **start, char **stop)
{
  int rem_points = 3;
  char *stop_old = *stop;

  while(rem_points > 0 && *start < stop_old)
  {
    while(**start != '.' && *start < stop_old)
    {
      (*start)++;
    }
    (*start)++;
    rem_points--;
  }
  *stop = *start;
  while(**stop != '"' && *stop < stop_old) (*stop)++;
}

static void update_prop_file(char *filename, int rev)
{
  char line[512];
  char fileout[128];
  char *out;
  char *outComment;
  char *outVersion;
  char *outFileVersion;
  FILE *fsource;
  FILE *fdest;

  if(filename == NULL || strlen(filename) == 0)
  {
    return;
  }

  if((fsource = fopen(filename, "r+")) == NULL)
  {
    fprintf(stderr, "Failed to open output file '%s'\n", filename);
    return;
  }

  sprintf(fileout, "%s.tmp", filename);
  if((fdest = fopen(fileout, "w")) == NULL)
  {
    fprintf(stderr, "Failed to create output file '%s'\n", fileout);
    fclose(fsource);
    return;
  }

  /*
  look for these two lines:
  [assembly: AssemblyVersion("1.2.3.4")]
  [assembly: AssemblyFileVersion("1.2.3.4")]
  */

  while((out = fgets(line, 512, fsource)) != NULL)
  {
    outComment = strstr(line, "//");
    outVersion = strstr(line, "AssemblyVersion");
    outFileVersion = strstr(line, "AssemblyFileVersion");
    if(outVersion == NULL && outFileVersion == NULL)
    {
      fputs(line, fdest);
      continue;
    }

    if(outVersion != NULL)
    {
      if(outComment == NULL || outComment > outVersion)
      {
        char *start, *stop;
        char buffer[512];
        start = outVersion + strlen("AssemblyVersion") + 2;
        stop = start + strlen(start) - 1;
        get_buildno_ptr(&start, &stop);
        strcpy(buffer, stop);
        sprintf(start, "%i%s", rev, buffer);
        fputs(line, fdest);
      }
      else
      {
        fputs(line, fdest);
        continue;
      }
    }

    if(outFileVersion != NULL)
    {
      if(outComment == NULL || outComment > outFileVersion)
      {
        char *start, *stop;
        char buffer[512];
        start = outFileVersion + strlen("AssemblyFileVersion") + 2;
        stop = start + strlen(start) - 1;
        get_buildno_ptr(&start, &stop);
        strcpy(buffer, stop);
        sprintf(start, "%i%s", rev, buffer);
        fputs(line, fdest);
      }
      else
      {
        fputs(line, fdest);
        continue;
      }
    }
  }

  fclose(fsource);
  fclose(fdest);

  remove(filename);
  rename(fileout, filename);

  (void)out;  /* to avoid a compiler warning */
}

int main(int argc, char *argv[])
{
  char *outname = NULL;
  FILE *fp;
  int index;
  int process_self = 0;
  int verbose = 0;
  int ignore_eol = 0;
  int max_build, accum_build;
  int max_year, max_month, max_day;
  int ismodified, filemodified, modificationcount;
  char prefix[MAX_SYMBOLLENGTH], suffix[MAX_SYMBOLLENGTH];
  char *const_basename = NULL;
  char modified_suffix[2];
  int output_language = LANG_C_CPP;   /* flag for C/C++, Java or C# output */
  char *namespace_name = NULL;        /* java package or C# namespace to put revision info in */
  char *propname = NULL;
  char *startComment = NULL;
  char *endComment = NULL;
  char *continueComment = NULL;

  if (argc <= 1)
    about();

  /* collect the options */
  prefix[0] = '\0';
  suffix[0] = '\0';

  for (index = 1; index < argc; index++) {
    /* check for options */
    if (argv[index][0] == '-'
#if defined __WIN32__ || defined _Win32 || defined _WIN32
     || argv[index][0] == '/'
#endif
    )
    {
      switch (argv[index][1]) {
      case 'f': {
        size_t len;
        char *ptr = strchr(&argv[index][2], '#');
        len = (ptr != NULL) ? (int)(ptr - &argv[index][2]) : (int)strlen(&argv[index][2]);
        if (len >= MAX_SYMBOLLENGTH)
          len = MAX_SYMBOLLENGTH - 1;
        strncpy(prefix, &argv[index][2], len);
        prefix[len] = '\0';
        ptr = (ptr != NULL) ? ptr + 1 : strchr(argv[index], '\0');
        len = strlen(ptr);
        if (len >= MAX_SYMBOLLENGTH)
          len = MAX_SYMBOLLENGTH - 1;
        strncpy(suffix, ptr, len);
        suffix[len] = '\0';
        break;
      } /* case */
      case 'i':
        process_self = 1;
        break;
      case 'j':
        output_language = LANG_JAVA;
        namespace_name = &argv[index][2];
        break;
      case 'c':
        output_language = LANG_CS;
        namespace_name = &argv[index][2];
        break;
      case 'P':
        output_language = LANG_PROP;
        break;
      case 'O':
        output_language = LANG_ORACLE;
        namespace_name = &argv[index][2];
        break;
      case 'p':
        propname = &argv[index][2];
        break;
      case 'm':
        const_basename = &argv[index][2];
        break;
      case 'o':
        outname = &argv[index][2];
        break;
      case 'v':
        verbose = 1;
        break;
      case 'n':
        ignore_eol = 1;
        break;
      default:
        fprintf(stderr, "Invalid option '%s'\n", argv[index]);
        about();
      } /* switch */
    } /* if */
  } /* for */

  if (outname == NULL) {
    switch (output_language) {
    case LANG_C_CPP:
      outname = "svnrev.h";
      break;
    case LANG_JAVA:
      outname = "SvnRevision.java";
      break;
    case LANG_CS:
      outname = "SvnRevision.cs";
      break;
    case LANG_PROP:
      outname = "svnrev.property";
      break;
    case LANG_ORACLE:
      outname = "svnrev.sql";
      break;
    } /* switch */
  } /* if */

  /* scan through the existing svnrev file to find the "build" count */
  if (const_basename == NULL || *const_basename == '\0')
    const_basename = "SVN_REV";
  assert(outname != NULL);
  modificationcount = scanmodifications(outname, const_basename);
  if (!process_self && *outname != '\0')
    remove(outname);

  /* phase 1: scan through all files and get the highest build number */
  max_build = 0;
  accum_build = 0;      /* for RCS / CVS */
  max_year = max_month = max_day = 0;
  ismodified = 0;
  for (index = 1; index < argc; index++) {
    /* skip the options (already handled) */
    if (argv[index][0] == '-'
#if defined __WIN32__ || defined _Win32 || defined _WIN32
     || argv[index][0] == '/'
#endif
    )
      continue;

    filemodified = 0;
    if (strcasecmp(argv[index], outname)!=0)
      processfile(argv[index], 0, ignore_eol, &max_build, &accum_build, &max_year, &max_month, &max_day, &filemodified);
    if (filemodified && verbose)
      fprintf(stderr, "\tNotice: modified file '%s'\n", argv[index]);
    ismodified = ismodified || filemodified;
  } /* for */

  /* also run over the existing header file, if any */
  if (process_self && *outname != '\0')
    processfile(outname, 1, ignore_eol, &max_build, &accum_build, &max_year, &max_month, &max_day, NULL/*&ismodified*/);

  if (accum_build > max_build)
    max_build = accum_build;
  modified_suffix[0] = (char)(ismodified ? 'M' : '\0');
  modified_suffix[1] = '\0';

  /* phase 2: write a file with this highest build number */
  if (*outname == '\0') {
    fp = stdout;
  } else if ((fp = fopen(outname, "w")) == NULL) {
    fprintf(stderr, "Failed to create output file '%s'\n", outname);
    return 2;
  } /* if */
  if (*outname != '\0') {
    switch (output_language) {
    case LANG_C_CPP:
    case LANG_JAVA:
    case LANG_CS:
      startComment = "/* ";
      endComment = " */";
      continueComment = " * ";
      break;
    case LANG_ORACLE:
      startComment = "-- ";
      endComment = "-- ";
      continueComment = "-- ";
      break;
    case LANG_PROP:
      startComment = "# ";
      endComment = "# ";
      continueComment = "# ";
      break;
    } /* switch */
    /* don't print the comments to stdout */
    fprintf(fp, "%sThis file was generated by the \"svnrev\" utility\n"
                "%s(http://www.compuphase.com/svnrev.htm).\n"
                "%sYou should not modify it manually, as it may be re-generated.\n"
                "%s\n"
                "%s$Revision: %d%s$\n"
                "%s$Date: %04d-%02d-%02d$\n"
                "%s\n\n", startComment, continueComment, continueComment, continueComment,
                  continueComment,  max_build,
              modified_suffix, continueComment, max_year, max_month, max_day,endComment);
  } /* if */

  switch (output_language) {
  case LANG_C_CPP:
    fprintf(fp, "#ifndef %s_H\n", const_basename);
    fprintf(fp, "#define %s_H\n\n", const_basename);
    fprintf(fp, "#define %s\t\t%d\n", const_basename, max_build);
    fprintf(fp, "#define %sSTR\t\"%s%d%s%s\"\n", const_basename, prefix, max_build, modified_suffix, suffix);
    fprintf(fp, "#define %sDATE\t\"%04d-%02d-%02d\"\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "#define %sSTAMP\t%04d%02d%02dL\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "#define %sMODIFIED\t%d\n", const_basename, ismodified ? modificationcount + 1 : 0);
    fprintf(fp, "\n#endif /* %s_H */\n", const_basename);
    break;
  case LANG_JAVA:
    if (namespace_name == NULL || *namespace_name == '\0')
      namespace_name = "com.compuphase";
    fprintf(fp, "package %s;\n\n", namespace_name);
    fprintf(fp, "public interface SvnRevision\n");
    fprintf(fp, "{\n");
    fprintf(fp, "    public final static int %s = %d;\n", const_basename, max_build);
    fprintf(fp, "    public final static String %sSTR = \"%s%d%s%s\";\n", const_basename, prefix, max_build, modified_suffix, suffix);
    fprintf(fp, "    public final static String %sDATE = \"%04d-%02d-%02d\";\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "    public final static long %sSTAMP = %04d%02d%02dL;\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "    public final static int %sMODIFIED = %d;\n", const_basename, ismodified ? modificationcount + 1 : 0);
    fprintf(fp, "}\n\n");
    break;
  case LANG_CS: {
    char *cls;
    if (namespace_name == NULL || (cls = strrchr(namespace_name, '.')) == NULL) {
      namespace_name = "compuphase";
      cls = "svnrev";
    } else {
      *cls++ = '\0';
    } /* if */
    fprintf(fp, "using System;\n\n");
    fprintf(fp, "namespace %s\n{\n", namespace_name);
    fprintf(fp, "\tpublic static partial class %s\n", cls);
    fprintf(fp, "\t{\n");
    fprintf(fp, "\t\tpublic const Int32  %s = %d;\n", const_basename, max_build);
    fprintf(fp, "\t\tpublic const String %sSTR = \"%s%d%s%s\";\n", const_basename, prefix, max_build, modified_suffix, suffix);
    fprintf(fp, "\t\tpublic const String %sDATE = \"%04d-%02d-%02d\";\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "\t\tpublic const Int64  %sSTAMP = %04d%02d%02d;\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "\t\tpublic const Int32  %sMODIFIED = %d;\n", const_basename, ismodified ? modificationcount + 1 : 0);
    fprintf(fp, "\t}\n}\n");
    update_prop_file(propname, max_build);
    break;
  } /* case */
  case LANG_PROP:
    fprintf(fp, "%s = %d\n", const_basename, max_build);
    fprintf(fp, "%sSTR = %s%d%s%s\n", const_basename, prefix, max_build, modified_suffix, suffix);
    fprintf(fp, "%sDATE = %04d-%02d-%02d\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "%sSTAMP = %04d%02d%02dL\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "%sMODIFIED = %d\n", const_basename, ismodified ? modificationcount + 1 : 0);
    break;
  case LANG_ORACLE:
    if (namespace_name == NULL || *namespace_name == '\0')
      namespace_name = "compuphase";
    fprintf( fp, "create or replace package %s\n", namespace_name);
    fprintf( fp, "as\n" );
    fprintf(fp, "\t%s CONSTANT INT := %d;\n", const_basename, max_build);
    fprintf(fp, "\t%sSTR CONSTANT VARCHAR2(200) := \"%s%d%s%s\";\n", const_basename, prefix, max_build, modified_suffix, suffix);
    fprintf(fp, "\t%sDATE CONSTANT DATE := TO_DATE('%04d-%02d-%02d','YYYY-MM-DD');\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "\t%sSTAMP CONSTANT INT := %04d%02d%02d;\n", const_basename, max_year, max_month, max_day);
    fprintf(fp, "\t%sMODIFIED CONSTANT INT := %d;\n", const_basename, ismodified ? modificationcount + 1 : 0);
    fprintf(fp, "END %s;\n/\n\nshow errors\n\n", namespace_name);
    break;
  } /* switch */

  if (*outname != '\0')
    fclose(fp);

  return 0;
}
