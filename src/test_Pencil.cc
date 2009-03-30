#define BOOST_TEST_MODULE $Id$

#include "config.h"

#include <boost/test/included/unit_test.hpp>

#include "Pencil.h"

int add( int i, int j ) {
    return i+j;
}

BOOST_AUTO_TEST_CASE( declare_pointer )
{
  pecos::suzerain::Pencil<> *p = NULL;
}
