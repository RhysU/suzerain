# ===========================================================================
#            http://autoconf-archive.cryp.to/ax_add_am_latex.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_ADD_AM_LATEX
#
# DESCRIPTION
#
#  Using the AX_AM_MACROS framework, add TeX- and LaTeX-based
#  DVI/PDF capabilities to every Automake file containing @INC_AMINCLUDE@.
#
#   LATEX_FLAGS, TEX_FLAGS, PDFLATEX_FLAGS, and PDFTEX_FLAGS
#   may be set to modify document processing.
#
# LAST MODIFICATION
#
#   2009-09-14
#
# COPYLEFT
#
#   Copyright (c) 2008 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_ADD_AM_LATEX],[

AC_REQUIRE([AC_CONFIG_AUX_DIR_DEFAULT])
AC_REQUIRE_AUX_FILE([tex-it])
AC_SUBST([tex_it_dir],[`basename "$ac_aux_dir"`])

AC_REQUIRE([AX_AM_MACROS])
AC_REQUIRE([ACLTX_PROG_TEX])
AC_REQUIRE([ACLTX_PROG_PDFTEX])
AC_REQUIRE([ACLTX_PROG_LATEX])
AC_REQUIRE([ACLTX_PROG_PDFLATEX])
AC_REQUIRE([ACLTX_PROG_BIBTEX])
AC_REQUIRE([ACLTX_PROG_MAKEINDEX])

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

AX_ADD_AM_MACRO([
# Added by AX_ADD_AM_LATEX m4 macro
%%.dvi : %%.tex
	if file \$< | grep LaTeX &>/dev/null ; then \\
		\$(abs_top_srcdir)/\$(tex_it_dir)/tex-it \$(latex) \$(LATEX_FLAGS) \$< ;\\
	else \\
		\$(abs_top_srcdir)/\$(tex_it_dir)/tex-it \$(TEX) \$(TEX_FLAGS) \$< ;\\
	fi

%%.pdf : %%.tex
	if file \$< | grep LaTeX &>/dev/null ; then \\
		\$(abs_top_srcdir)/\$(tex_it_dir)/tex-it \$(pdflatex) \$(PDFLATEX_FLAGS) \$< ;\\
	else \\
		\$(abs_top_srcdir)/\$(tex_it_dir)/tex-it \$(pdftex) \$(PDFTEX_FLAGS) \$< ;\\
	fi
])
])
