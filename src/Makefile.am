## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(GSL_CFLAGS) $(BOOST_CPPFLAGS) $(MATIO_CFLAGS)
AM_LDFLAGS     = $(GSL_LIBS)   $(BOOST_LDFLAGS)  $(MATIO_LIBS) 
bin_PROGRAMS   =
check_PROGRAMS = 

bin_PROGRAMS        += helloworld
helloworld_SOURCES   = helloworld.cc

check_PROGRAMS      += test_example
test_example_SOURCES = test_example.cc
test_example_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
