## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_CPPFLAGS   += -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(MKL_CFLAGS) $(GSL_CFLAGS) $(LOG4CXX_CFLAGS)
AM_CPPFLAGS   += $(BOOST_CPPFLAGS) $(P3DFFT_CFLAGS) $(FFTW3_THREADS_CFLAGS)
AM_CPPFLAGS   += $(EIGEN_CFLAGS) $(GRVY_CFLAGS) $(ESIO_CFLAGS)
AM_LDFLAGS     = $(MPIP_LDFLAGS)
AM_LDFLAGS    += $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
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

# Test progams and/or performance drivers
bin_PROGRAMS                   += driver_suzerain_p3dfft
driver_suzerain_p3dfft_SOURCES  = driver_suzerain_p3dfft.cpp
driver_suzerain_p3dfft_LDADD    = ../suzerain/libsuzerain.la
driver_suzerain_p3dfft_LDADD   += $(MPIP_LIBS)
driver_suzerain_p3dfft_LDADD   += $(GSL_LIBS) $(LOG4CXX_LIBS)
driver_suzerain_p3dfft_LDADD   += $(P3DFFT_LIBS) $(FFTW3_THREADS_LIBS)
driver_suzerain_p3dfft_LDADD   += $(P3DFFT_LIBS) $(FFTW3_THREADS_LIBS)
driver_suzerain_p3dfft_LDADD   += $(BOOST_PROGRAM_OPTIONS_LIBS)

# Applications
bin_PROGRAMS       += channel_ex
channel_ex_SOURCES  = channel_ex.cpp
channel_ex_LDADD    = ../suzerain/libsuzerain.la
channel_ex_LDADD   += $(ESIO_LIBS) $(GRVY_LIBS)
channel_ex_LDADD   += $(GSL_LIBS) $(LOG4CXX_LIBS)
channel_ex_LDADD   += $(P3DFFT_LIBS) $(FFTW3_THREADS_LIBS)
channel_ex_LDADD   += $(BOOST_IOSTREAMS_LIBS)
channel_ex_LDADD   += $(BOOST_PROGRAM_OPTIONS_LIBS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
