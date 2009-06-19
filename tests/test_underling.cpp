#define BOOST_TEST_MODULE $Id$

#include <config.h>

#include <boost/test/included/unit_test.hpp>
#include <suzerain/underling.h>

BOOST_AUTO_TEST_CASE( alloc )
{
    underling_workspace * w = underling_workspace_alloc(3);
    underling_workspace_free(w);
}
