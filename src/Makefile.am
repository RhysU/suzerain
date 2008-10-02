## Process this file with automake to produce Makefile.in

AM_CPPFLAGS    = -I$(top_builddir)/lib -I$(top_srcdir)/lib
AM_CPPFLAGS   += $(GSL_CFLAGS) $(BOOST_CPPFLAGS) $(MATIO_CFLAGS)
AM_LDFLAGS     = $(GSL_LIBS)   $(BOOST_LDFLAGS)  $(MATIO_LIBS) 
bin_PROGRAMS   =
check_PROGRAMS = 

# Progams
bin_PROGRAMS        += helloworld
helloworld_SOURCES   = helloworld.cc

bin_PROGRAMS     += burgers
burgers_SOURCES   = burgers.cc solver.cc residual.cc legendrePoly.cc tools.cc
burgers_LDADD     = $(BOOST_PROGRAM_OPTIONS_LIB)

# Unit test executables
check_PROGRAMS      += test_example
test_example_SOURCES = test_example.cc
test_example_LDADD   = $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

check_PROGRAMS            += legendrePoly_test
legendrePoly_test_SOURCES  = legendrePoly_test.cc legendrePoly.cc
legendrePoly_test_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

check_PROGRAMS        += residual_test
residual_test_SOURCES  = residual_test.cc residual.cc legendrePoly.cc
residual_test_LDADD    = $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
