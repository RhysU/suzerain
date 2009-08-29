#define BOOST_TEST_MODULE $Id$

#include <config.h>

#include <boost/test/included/unit_test.hpp>
#include <suzerain/underling.h>

BOOST_AUTO_TEST_CASE( alloc )
{
    /* In conjunction with valgrind, sanity checks allocation/deallocation */
    for (int ndim = 1; ndim <= 4; ++ndim) {
        for (int nstage = 1; nstage <= ndim+1; ++nstage) {
            underling_workspace * const w
                = underling_workspace_alloc(ndim, nstage);
            BOOST_CHECK_EQUAL(w->ndim, ndim);
            BOOST_CHECK_EQUAL(w->nstage, nstage);
            underling_workspace_free(w);
        }
    }
}
