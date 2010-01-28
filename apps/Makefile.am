## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(MKL_CFLAGS) $(GSL_CFLAGS) $(LOG4CXX_CFLAGS)
AM_CPPFLAGS   += $(BOOST_CPPFLAGS) $(P3DFFT_CFLAGS) $(FFTW3_CFLAGS)
AM_CPPFLAGS   += $(EIGEN_CFLAGS) $(HPCT_CFLAGS)
AM_LDFLAGS     = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
AM_LDFLAGS    += $(BOOST_IOSTREAMS_LDFLAGS)
bin_PROGRAMS   = # Append below
check_PROGRAMS = # Append below
CLEANFILES     = # Append below

# FIXME: Ugly, separate MPI compilation into separate dir, maybe
CXX=@MPICXX@

# Clean up from Intel Compiler's -prof-gen output
CLEANFILES += pgopti.spi pgopti.spl *.dyn *.dpi
# Clean up from Intel Compiler's precompiled header output
CLEANFILES += *.pchi

# Progams
bin_PROGRAMS                   += driver_suzerain_p3dfft
driver_suzerain_p3dfft_SOURCES  = driver_suzerain_p3dfft.cpp
driver_suzerain_p3dfft_LDADD    = ../suzerain/libsuzerain.la
driver_suzerain_p3dfft_LDADD   += $(GSL_LIBS) $(LOG4CXX_LIBS)
driver_suzerain_p3dfft_LDADD   += $(P3DFFT_LIBS) $(FFTW3_LIBS)
driver_suzerain_p3dfft_LDADD   += $(BOOST_PROGRAM_OPTIONS_LIBS)

bin_PROGRAMS          += driver_pencil
driver_pencil_SOURCES  = driver_pencil.cpp
driver_pencil_LDADD    = ../suzerain/libsuzerain.la
driver_pencil_LDADD   += $(GSL_LIBS) $(LOG4CXX_LIBS)
driver_pencil_LDADD   += $(FFTW3_LIBS)
driver_pencil_LDADD   += $(BOOST_PROGRAM_OPTIONS_LIBS)

bin_PROGRAMS             += driver_underling
driver_underling_SOURCES  = driver_underling.cpp
driver_underling_LDADD    = ../suzerain/libsuzerain.la
driver_underling_LDADD   += ../suzerain/libunderling.la
driver_underling_LDADD   += -lfftw3_mpi
driver_underling_LDADD   += $(FFTW3_LIBS) $(HPCT_LIBS)
driver_underling_LDADD   += $(GSL_LIBS) $(LOG4CXX_LIBS)
driver_underling_LDADD   += $(BOOST_IOSTREAMS_LIBS)
driver_underling_LDADD   += $(BOOST_PROGRAM_OPTIONS_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
