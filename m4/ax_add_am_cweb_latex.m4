# ===========================================================================
#            http://autoconf-archive.cryp.to/ax_add_am_cweb_latex.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_ADD_AM_CWEB_LATEX
#
# DESCRIPTION
#
#   Using the AX_AM_MACROS framework, add CWEB compilation and 
#   CWEB/TeX and CWEB/LaTeX to DVI/PDF capabilities to every 
#   Automake file containing @INC_AMINCLUDE@.  Requires that the
#   LaTeX cweb package be usable.
#   
#   By default, rules include options to CWEB programs to reduce verbosity,
#   this can be tweaked by modifying CWEAVE_FLAGS and CTANGLE_FLAGS.
#   LATEX_FLAGS, TEX_FLAGS, PDFLATEX_FLAGS, and PDFTEX_FLAGS
#   may be set to modify document processing.
#
# LAST MODIFICATION
#
#   2008-08-07
#
# COPYLEFT
#
#   Copyright (c) 2008 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_ADD_AM_CWEB_LATEX],[

AC_REQUIRE([AC_CONFIG_AUX_DIR_DEFAULT])
AC_REQUIRE_AUX_FILE([tex-it])
AC_SUBST([tex_it_dir],[$ac_aux_dir])

AC_REQUIRE([AX_AM_MACROS])
AC_REQUIRE([ACLTX_PROG_TEX])
AC_REQUIRE([ACLTX_PROG_PDFTEX])
AC_REQUIRE([ACLTX_PROG_LATEX])
AC_REQUIRE([ACLTX_PROG_PDFLATEX])
AC_REQUIRE([AX_PROG_CWEAVE])
AC_REQUIRE([AX_PROG_CTANGLE])
ACLTX_CLASS_CWEB(cweb_latex,,AC_MSG_ERROR([LaTeX package cweb not working.]))

AC_ARG_VAR(CWEAVE_FLAGS,[Options for cweave used in Makefile rules])
if test "x$CWEAVE_FLAGS" = x; then
    CWEAVE_FLAGS=-bhp          # Default to quiet mode
fi
AC_ARG_VAR(CTANGLE_FLAGS,[Options for ctangle used in Makefile rules])
if test "x$CTANGLE_FLAGS" = x; then
    CTANGLE_FLAGS=-bhp         # Default to quiet mode
fi
AC_ARG_VAR(LATEX_FLAGS,[Options for latex used in Makefile rules])
if test "x$LATEX_FLAGS" = x; then
    LATEX_FLAGS=-src-specials  # Default to DVI specials insertion
fi
AC_ARG_VAR(TEX_FLAGS,[Options for tex used in Makefile rules])
if test "x$TEX_FLAGS" = x; then
    TEX_FLAGS=
fi
AC_ARG_VAR(PDFLATEX_FLAGS,[Options for pdflatex used in Makefile rules])
if test "x$PDFLATEX_FLAGS" = x; then
    PDFLATEX_FLAGS=
fi
AC_ARG_VAR(PDFTEX_FLAGS,[Options for pdftex used in Makefile rules])
if test "x$PDFTEX_FLAGS" = x; then
    PDFTEX_FLAGS=
fi

dnl TODO Does not handle .bib files when building in VPATH
AX_ADD_AM_MACRO([
%%.c %%.h: %%.w
	\$(CTANGLE) \$(CTANGLE_FLAGS) \$< - \@S|@@

%%.cc %%.h: %%.w
	\$(CTANGLE) \$(CTANGLE_FLAGS) \$< - \@S|@@

%%.tex : %%.w
	\$(CWEAVE) \$(CWEAVE_FLAGS) \$< - \@S|@@
	echo >> \@S|@@
	echo %% Workaround bug with CWEB and older pdfeTeX 1.2x releases. >> \@S|@@
	echo \\\\\\\\fi\\\\\\\\fin >> \@S|@@

%%.dvi : %%.tex
	if file \$< | grep LaTeX &>/dev/null ; then \\
		\$(top_builddir)/\$(tex_it_dir)/tex-it \$(latex) \$(LATEX_FLAGS) \$< ;\\
	else \\
		\$(top_builddir)/\$(tex_it_dir)/tex-it \$(TEX) \$(TEX_FLAGS) \$< ;\\
	fi

%%.pdf : %%.tex
	if file \$< | grep LaTeX &>/dev/null ; then \\
		\$(top_builddir)/\$(tex_it_dir)/tex-it \$(pdflatex) \$(PDFLATEX_FLAGS) \$< ;\\
	else \\
		\$(top_builddir)/\$(tex_it_dir)/tex-it \$(pdftex) \$(PDFTEX_FLAGS) \$< ;\\
	fi
])
])
