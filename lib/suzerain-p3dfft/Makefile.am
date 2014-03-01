## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir)

# MPI-based compilation and linking
CC  = @MPICC@
CXX = @MPICXX@
FC  = @MPIFC@

# Building in parallel is too brittle on Intel and compilation is fast enough.
.NOTPARALLEL:

# Distribute the following files
EXTRA_DIST     = COPYING COPYRIGHT README
noinst_HEADERS = p3dfft_internal.h p3dfft_s.h p3dfft_d.h

# Cannot have per-target flags otherwise it breaks dependencies like
#   fft_spec.$(OBJEXT):  $(srcdir)/fft_spec.F90
# below!  See "Per-Object Flags" in the Automake manual
AM_CPPFLAGS     += $(FFTW3_CFLAGS)  # Make FFTW available for compilation
AM_CPPFLAGS     += -DFFTW           # Always use FFTW for the FFTs
AM_CPPFLAGS     += -DSTRIDE1        # Always use stride one packing option
AM_LIBTOOLFLAGS  = --tag=FC         # Link with F90 to get Fortran MPI support

# Build a convenience library containing our one-off version of P3DFFT 2.3.
# The Fortran runtime libraries (i.e. $(FCLIBS)) are explicitly not added.
noinst_LTLIBRARIES             = libsuzerain-p3dfft.la
libsuzerain_p3dfft_la_SOURCES  = fft_exec.F90
libsuzerain_p3dfft_la_SOURCES += fft_init.F90
libsuzerain_p3dfft_la_SOURCES += fft_spec.F90
libsuzerain_p3dfft_la_SOURCES += module.F90
libsuzerain_p3dfft_la_SOURCES += wrap.F90
libsuzerain_p3dfft_la_LIBADD   = $(FFTW3_LIBS)

# Get rid of module files on 'make clean'
CLEANFILES  = fft_spec.$(FC_MODEXT) p3dfft.$(FC_MODEXT)
CLEANFILES += *__genmod.f90 *__genmod.$(FC_MODEXT)

# Explicitly provide object dependencies to avoid computing them.
# These are catalogued in case someone later wants to build regular objects.
fft_spec.$(OBJEXT):  $(srcdir)/fft_spec.F90
module.$(OBJEXT):    $(srcdir)/module.F90   fft_spec.$(OBJEXT)
wrap.$(OBJEXT):      $(srcdir)/wrap.F90     module.$(OBJEXT)
fft_exec.$(OBJEXT):  $(srcdir)/fft_exec.F90 fft_spec.$(OBJEXT) module.$(OBJEXT)

# Explicitly provide libtool object dependencies to avoid computing them.
# These are identical to those above but are required for a LTLIBRARIES target.
fft_spec.lo:         $(srcdir)/fft_spec.F90
module.lo:           $(srcdir)/module.F90   fft_spec.lo
wrap.lo:             $(srcdir)/wrap.F90     module.lo
fft_exec.lo:         $(srcdir)/fft_exec.F90 fft_spec.lo module.lo

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
