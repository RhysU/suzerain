## Process this file with automake to produce Makefile.in

# The source directory layout is as follows:
#
#     - suzerain/*          Logic that is tested by tests/ and which could be
#                           conceivably picked up for use outside of Suzerain
#                           because it has minimal third-party dependencies.
#
#     - suzerain/support/*  Logic that builds atop suzerain/* and which is used
#                           and tested by primarily by apps/.  It would be
#                           difficult to use outside of Suzerain because of its
#                           third-party dependencies.
#
# All logic is installed as libsuzerain.la using non-recursive Automake.
# This final target is built from objects and multiple convenience libraries.

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_LDFLAGS     =

BUILT_SOURCES             = # Append below
CLEANFILES                = # Append below
lib_LTLIBRARIES           = # Append below
nobase_pkginclude_HEADERS = # Append below
noinst_LTLIBRARIES        = # Append below

# MPI-based compilation and linking
CXX = @MPICXX@
CC  = @MPICC@

BUILT_SOURCES += suzerain-config.h
BUILT_SOURCES += suzerain-revision.h

suzerain-revision.h: $(srcdir)/*.h                            \
                     $(srcdir)/*.c                            \
                     $(srcdir)/*.hpp  $(srcdir)/support/*.hpp \
                     $(srcdir)/*.cpp  $(srcdir)/support/*.cpp \
                     $(srcdir)/*.def                          \
                     $(srcdir)/Makefile*
	$(top_builddir)/build-aux/svnrev -n -o$@ -mSUZERAIN_REVISION_ "-f$(PACKAGE_VERSION).#" -i $?

# Public API headers
nobase_pkginclude_HEADERS += allocator.hpp
nobase_pkginclude_HEADERS += blas_et_al.h
nobase_pkginclude_HEADERS += blas_et_al.hpp
nobase_pkginclude_HEADERS += bl.h
nobase_pkginclude_HEADERS += blasius.h
nobase_pkginclude_HEADERS += bsmbsm.h
nobase_pkginclude_HEADERS += bsmbsm_solver.hpp
nobase_pkginclude_HEADERS += bspline.h
nobase_pkginclude_HEADERS += bspline.hpp
nobase_pkginclude_HEADERS += bsplineop.h
nobase_pkginclude_HEADERS += coalescing_pool.hpp
nobase_pkginclude_HEADERS += common.h
nobase_pkginclude_HEADERS += common.hpp
nobase_pkginclude_HEADERS += complex.h
nobase_pkginclude_HEADERS += complex.hpp
nobase_pkginclude_HEADERS += constraint.hpp
nobase_pkginclude_HEADERS += constraint_treatment.hpp
nobase_pkginclude_HEADERS += contiguous_memory.hpp
nobase_pkginclude_HEADERS += countof.h
nobase_pkginclude_HEADERS += diffwave.h
nobase_pkginclude_HEADERS += diffwave.hpp
nobase_pkginclude_HEADERS += error.h
nobase_pkginclude_HEADERS += exprgrammar.hpp
nobase_pkginclude_HEADERS += exprparse.hpp
nobase_pkginclude_HEADERS += fftw.hpp
nobase_pkginclude_HEADERS += filterop.h
nobase_pkginclude_HEADERS += format.hpp
nobase_pkginclude_HEADERS += functional.hpp
nobase_pkginclude_HEADERS += function.h
nobase_pkginclude_HEADERS += gbdddddmv.h
nobase_pkginclude_HEADERS += gbddddmv.h
nobase_pkginclude_HEADERS += gbdddmv.h
nobase_pkginclude_HEADERS += gbddmv.h
nobase_pkginclude_HEADERS += gbdmv.h
nobase_pkginclude_HEADERS += gbmatrix.h
nobase_pkginclude_HEADERS += gbmv.h
nobase_pkginclude_HEADERS += grid_specification.hpp
nobase_pkginclude_HEADERS += htstretch.h
nobase_pkginclude_HEADERS += hybrid_residual_operator.hpp
nobase_pkginclude_HEADERS += inorder.h
nobase_pkginclude_HEADERS += inorder.hpp
nobase_pkginclude_HEADERS += isothermal_mass_operator.hpp
nobase_pkginclude_HEADERS += isothermal_specification.hpp
nobase_pkginclude_HEADERS += iterator.hpp
nobase_pkginclude_HEADERS += kahan.h
nobase_pkginclude_HEADERS += l2.hpp
nobase_pkginclude_HEADERS += largo_formulation.hpp
nobase_pkginclude_HEADERS += lowstorage.hpp
nobase_pkginclude_HEADERS += mass_operator.hpp
nobase_pkginclude_HEADERS += math.hpp
nobase_pkginclude_HEADERS += mpi_datatype.hpp
nobase_pkginclude_HEADERS += mpi.h
nobase_pkginclude_HEADERS += mpi.hpp
nobase_pkginclude_HEADERS += mpl.hpp
nobase_pkginclude_HEADERS += multi_array.hpp
nobase_pkginclude_HEADERS += ndx.hpp
nobase_pkginclude_HEADERS += noise_specification.hpp
nobase_pkginclude_HEADERS += operator_base.hpp
nobase_pkginclude_HEADERS += operator_tools.hpp
nobase_pkginclude_HEADERS += os.h
nobase_pkginclude_HEADERS += parcel.h
nobase_pkginclude_HEADERS += pencilfft.hpp
nobase_pkginclude_HEADERS += pencil_grid.hpp
nobase_pkginclude_HEADERS += pencil.hpp
nobase_pkginclude_HEADERS += physical_view.hpp
nobase_pkginclude_HEADERS += pre_gsl.h
nobase_pkginclude_HEADERS += reacting_imexop.h
nobase_pkginclude_HEADERS += rholt.hpp
nobase_pkginclude_HEADERS += rholut.hpp
nobase_pkginclude_HEADERS += rholut_imexop.h
nobase_pkginclude_HEADERS += radial_nozzle.h
nobase_pkginclude_HEADERS += richardson.h
nobase_pkginclude_HEADERS += rngstream.hpp
nobase_pkginclude_HEADERS += running_statistics.hpp
nobase_pkginclude_HEADERS += safe_bool.hpp
nobase_pkginclude_HEADERS += sbmv.h
nobase_pkginclude_HEADERS += shared_range.hpp
nobase_pkginclude_HEADERS += state_fwd.hpp
nobase_pkginclude_HEADERS += state.hpp
nobase_pkginclude_HEADERS += storage.hpp
nobase_pkginclude_HEADERS += support/application_base.hpp
nobase_pkginclude_HEADERS += support/definition_base.hpp
nobase_pkginclude_HEADERS += support/driver_base.hpp
nobase_pkginclude_HEADERS += support/esio_fwd.hpp
nobase_pkginclude_HEADERS += support/fftw_definition.hpp
nobase_pkginclude_HEADERS += support/field.hpp
nobase_pkginclude_HEADERS += support/grid_definition.hpp
nobase_pkginclude_HEADERS += support/isothermal_definition.hpp
nobase_pkginclude_HEADERS += support/largo_definition.hpp
nobase_pkginclude_HEADERS += support/loadable.hpp
nobase_pkginclude_HEADERS += support/logging.hpp
nobase_pkginclude_HEADERS += support/noise_definition.hpp
nobase_pkginclude_HEADERS += support/overridable.hpp
nobase_pkginclude_HEADERS += support/populatable.hpp
nobase_pkginclude_HEADERS += support/program_options.hpp
nobase_pkginclude_HEADERS += support/restart_definition.hpp
nobase_pkginclude_HEADERS += support/savable.hpp
nobase_pkginclude_HEADERS += support/signal_definition.hpp
nobase_pkginclude_HEADERS += support/statistics_definition.hpp
nobase_pkginclude_HEADERS += support/support.hpp
nobase_pkginclude_HEADERS += support/time_definition.hpp
nobase_pkginclude_HEADERS += suzerain-config.h
nobase_pkginclude_HEADERS += suzerain-revision.h
nobase_pkginclude_HEADERS += svehla.h
nobase_pkginclude_HEADERS += timecontroller.hpp
nobase_pkginclude_HEADERS += timers.h
nobase_pkginclude_HEADERS += traits.hpp
nobase_pkginclude_HEADERS += utility.hpp
nobase_pkginclude_HEADERS += validation.hpp
nobase_pkginclude_HEADERS += version.hpp
nobase_pkginclude_HEADERS += zgbsv_specification.hpp

# Public library is assembled from multiple convenience libraries
# combined with carefully-ordered external dependencies.
#
# The ordering of the external libraries on which libsuzerain.la depends places
# the most-uncommon prerequisities first in the hope that introduction of
# any system-wide /usr/lib, /usr/lib64, etc. is delayed as long as possible.
lib_LTLIBRARIES          += libsuzerain.la
libsuzerain_la_LDFLAGS    = -rpath $(libdir)
libsuzerain_la_LDFLAGS   += -release $(PACKAGE_VERSION)
libsuzerain_la_SOURCES    = version.cpp
libsuzerain_la_CPPFLAGS   = $(AM_CPPFLAGS)
libsuzerain_la_LIBADD     = libblasetal.la
libsuzerain_la_LIBADD    += libbspline.la
libsuzerain_la_LIBADD    += libbase.la
libsuzerain_la_LIBADD    += libexprparse.la
libsuzerain_la_LIBADD    += libprecise.la
libsuzerain_la_LIBADD    += libsupport.la
libsuzerain_la_LIBADD    += $(ESIO_LIBS)
libsuzerain_la_LIBADD    += $(UNDERLING_LIBS)
libsuzerain_la_LIBADD    += $(top_builddir)/lib/suzerain-p3dfft/libsuzerain-p3dfft.la
libsuzerain_la_LIBADD    += $(top_builddir)/largo/largo/liblargo.la
libsuzerain_la_LIBADD    += $(FCLIBS)   # Should be after libsuzerain-p3dfft and largo
libsuzerain_la_LIBADD    += $(FFTW3_MPI_THREADS_LIBS)
libsuzerain_la_LIBADD    += $(GRVY_LIBS)
libsuzerain_la_LIBADD    += $(BOOST_IOSTREAMS_LDFLAGS)
libsuzerain_la_LIBADD    += $(BOOST_IOSTREAMS_LIBS)
libsuzerain_la_LIBADD    += $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
libsuzerain_la_LIBADD    += $(BOOST_PROGRAM_OPTIONS_LIBS)
libsuzerain_la_LIBADD    += $(GSL_LIBS)  # Includes MKL via BLAS hooks
libsuzerain_la_LIBADD    += $(LOG4CXX_LIBS)

# BLAS and LAPACK wrappers, including custom BLAS- and LAPACK-like kernels.
# These routines ("BLAS et. al.") should be directly usable elsewhere
# caveat a Boost.Preprocessor dependency.
#
# When in doubt, check libblasetal_la_CPPFLAGS to see if something belongs.  If
# it won't compile without libbase_la_CPPFLAGS, it belongs in libbase.la.
noinst_LTLIBRARIES           += libblasetal.la
libblasetal_la_LDFLAGS        = # Convenience-only
libblasetal_la_CPPFLAGS       = $(libsuzerain_la_CPPFLAGS) # Append!
libblasetal_la_CPPFLAGS      += $(BOOST_CPPFLAGS)
libblasetal_la_CPPFLAGS      += $(MKL_CFLAGS)
libblasetal_la_SOURCES        = blas_et_al.c
libblasetal_la_SOURCES       += gbdddddmv.c
libblasetal_la_SOURCES       += gbddddmv.c
libblasetal_la_SOURCES       += gbdddmv.c
libblasetal_la_SOURCES       += gbddmv.c
libblasetal_la_SOURCES       += gbdmv.c
libblasetal_la_SOURCES       += gbmv.c
libblasetal_la_SOURCES       += sbmv.c
EXTRA_libblasetal_la_SOURCES  = dsgbsvx.def
EXTRA_libblasetal_la_SOURCES += gbdddddmv.def
EXTRA_libblasetal_la_SOURCES += gbddddmv.def
EXTRA_libblasetal_la_SOURCES += gbdddmv.def
EXTRA_libblasetal_la_SOURCES += gbddmv.def
EXTRA_libblasetal_la_SOURCES += gbdmv.def
EXTRA_libblasetal_la_SOURCES += gbmv.def
EXTRA_libblasetal_la_SOURCES += sbmv.def

# Functionality related to B-splines built atop the GSL.
# Includes a variety of banded operator functionality and compute kernels.
#
# When in doubt, check libbspline_la_CPPFLAGS to see if something belongs.  If
# it won't compile without libbase_la_CPPFLAGS, it belongs in libbase.la.
noinst_LTLIBRARIES          += libbspline.la
libbspline_la_LDFLAGS        = # Convenience-only
libbspline_la_CPPFLAGS       = $(libblasetal_la_CPPFLAGS) # Append!
libbspline_la_CPPFLAGS      += $(GSL_CFLAGS)
libbspline_la_SOURCES        = bsmbsm.c
libbspline_la_SOURCES       += bspline.c
libbspline_la_SOURCES       += bsplineop.c
libbspline_la_SOURCES       += error.c
libbspline_la_SOURCES       += filterop.c
libbspline_la_SOURCES       += gbmatrix.c
libbspline_la_SOURCES       += htstretch.c
libbspline_la_SOURCES       += pre_gsl.c
libbspline_la_SOURCES       += reacting_imexop.c
libbspline_la_SOURCES       += rholut_imexop.c
libbspline_la_SOURCES       += rholut_imexop00.c
EXTRA_libbspline_la_SOURCES  = bsmbsm_aPxpby_complex.def
EXTRA_libbspline_la_SOURCES += bsmbsm_aPxpby_real.def
EXTRA_libbspline_la_SOURCES += bsmbsm_pack.def
EXTRA_libbspline_la_SOURCES += filterop.def
EXTRA_libbspline_la_SOURCES += reacting_imexop.def
EXTRA_libbspline_la_SOURCES += rholut_imexop.def
EXTRA_libbspline_la_SOURCES += rholut_imexop00.def

# Arithmetic expression parsing functionality written as an outgrowth of
# http://agentzlerich.blogspot.com/2011/06/using-boost-spirit-21-to-evaluate.html
noinst_LTLIBRARIES       += libexprparse.la
libexprparse_la_LDFLAGS   = # Convenience-only
libexprparse_la_CPPFLAGS  = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS)
libexprparse_la_SOURCES   = exprgrammar.hpp
libexprparse_la_SOURCES  += exprparse_dchar.cpp
libexprparse_la_SOURCES  += exprparse_dstr.cpp
libexprparse_la_SOURCES  += exprparse_fchar.cpp
libexprparse_la_SOURCES  += exprparse_fstr.cpp
libexprparse_la_SOURCES  += exprparse.hpp
libexprparse_la_SOURCES  += exprparse_impl.hpp

# Logic with more extensive external dependencies, in particular
# depending on parallel decomposition details.
#
# When in doubt, check libbase_la_CPPFLAGS to see if something belongs.  If it
# won't compile without libsupport_la_CPPFLAGS, it belongs in libsupport.la.
noinst_LTLIBRARIES       += libbase.la
libbase_la_LDFLAGS        = # Convenience-only
libbase_la_CPPFLAGS       = $(libbspline_la_CPPFLAGS) # Append!
libbase_la_CPPFLAGS      += $(EIGEN_CFLAGS)
libbase_la_CPPFLAGS      += $(UNDERLING_CFLAGS)
libbase_la_CPPFLAGS      += $(FFTW3_MPI_THREADS_CFLAGS)
libbase_la_CPPFLAGS      += $(GRVY_CFLAGS)
libbase_la_SOURCES        = bl.c
libbase_la_SOURCES       += blasius.c
libbase_la_SOURCES       += bsmbsm_solver.cpp
libbase_la_SOURCES       += constraint.cpp
libbase_la_SOURCES       += diffwave.c
libbase_la_SOURCES       += fftw.cpp
libbase_la_SOURCES       += grid_specification.cpp
libbase_la_SOURCES       += hybrid_residual_operator.cpp
libbase_la_SOURCES       += inorder.c
libbase_la_SOURCES       += isothermal_mass_operator.cpp
libbase_la_SOURCES       += isothermal_specification.cpp
libbase_la_SOURCES       += l2.cpp
libbase_la_SOURCES       += largo_formulation.cpp
libbase_la_SOURCES       += mass_operator.cpp
libbase_la_SOURCES       += mpi.cpp
libbase_la_SOURCES       += noise_specification.cpp
libbase_la_SOURCES       += ndx.cpp
libbase_la_SOURCES       += operator_base.cpp
libbase_la_SOURCES       += os.c
libbase_la_SOURCES       += pencil_grid.cpp
libbase_la_SOURCES       += radial_nozzle.c
libbase_la_SOURCES       += richardson.c
libbase_la_SOURCES       += rngstream.cpp
libbase_la_SOURCES       += svehla.c
libbase_la_SOURCES       += zgbsv_specification.cpp

# Sources requiring IEEE-faithful floating point behavior
# Built separately to utilize separate optimization settings
noinst_LTLIBRARIES          += libprecise.la
libprecise_la_LDFLAGS        = # Convenience-only
libprecise_la_SOURCES        = kahan.c
EXTRA_libprecise_la_SOURCES  = kahan.def
libprecise_la_CFLAGS         = $(AM_CFLAGS)   $(CFLAGS_PRECISE)
libprecise_la_CXXFLAGS       = $(AM_CXXFLAGS) $(CXXFLAGS_PRECISE)
libprecise_la_CPPFLAGS       = $(libbase_la_CPPFLAGS) # Append!

# Sources requiring additional external dependencies
# Built separately to encourage decoupling versus other targets
noinst_LTLIBRARIES          += libsupport.la
libsupport_la_LDFLAGS        = # Convenience-only
libsupport_la_CPPFLAGS       = $(libbase_la_CPPFLAGS) # Append!
libsupport_la_CPPFLAGS      += $(ESIO_CFLAGS)
libsupport_la_CPPFLAGS      += $(LOG4CXX_CFLAGS)
libsupport_la_CPPFLAGS      += -I$(top_builddir)/largo -I$(top_srcdir)/largo
libsupport_la_SOURCES        = support/application_base.cpp
libsupport_la_SOURCES       += support/definition_base.cpp
libsupport_la_SOURCES       += support/driver_base.cpp
libsupport_la_SOURCES       += support/fftw_definition.cpp
libsupport_la_SOURCES       += support/field.cpp
libsupport_la_SOURCES       += support/grid_definition.cpp
libsupport_la_SOURCES       += support/isothermal_definition.cpp
libsupport_la_SOURCES       += support/largo_definition.cpp
libsupport_la_SOURCES       += support/logging.cpp
libsupport_la_SOURCES       += support/noise_definition.cpp
libsupport_la_SOURCES       += support/overridable.cpp
libsupport_la_SOURCES       += support/populatable.cpp
libsupport_la_SOURCES       += support/program_options.cpp
libsupport_la_SOURCES       += support/restart_definition.cpp
libsupport_la_SOURCES       += support/signal_definition.cpp
libsupport_la_SOURCES       += support/statistics_definition.cpp
libsupport_la_SOURCES       += support/support.cpp
libsupport_la_SOURCES       += support/time_definition.cpp
EXTRA_libsupport_la_SOURCES  =

# Clean up from Intel Compiler's -prof-gen output
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi

# Clean up from buildaux/svnrev-generated headers
CLEANFILES += *_svnrev.h

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
