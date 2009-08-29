#define BOOST_TEST_MODULE $Id$

#include <config.h>

#include <boost/test/included/unit_test.hpp>
#include <suzerain/underling.h>

BOOST_AUTO_TEST_CASE( alloc_and_prepare_workspace )
{
    const int ndim                = 3;
    underling_workspace * const w = underling_workspace_alloc(ndim);
    BOOST_CHECK_EQUAL(w->ndim, ndim);

    const int physical_size[ndim] = { 16, 24, 32 };
    underling_prepare_physical_size(w, physical_size);
    for (int i = 0; i < w->ndim; ++i) {
        BOOST_CHECK_EQUAL(w->dim_p[i].global_size, physical_size[i]);
    }

    underling_prepare_link(w, 0, 1);
    BOOST_CHECK_EQUAL(*(w->dim_p[0].transformed), w->dim_w[1]);
    BOOST_CHECK_EQUAL(*(w->dim_w[1].transformed), w->dim_p[0]);
    BOOST_CHECK_EQUAL(*((w->dim_p[0].transformed)->transformed), w->dim_p[0]);

    underling_prepare_link(w, 1, 2);
    BOOST_CHECK_EQUAL(*(w->dim_p[1].transformed), w->dim_w[2]);
    BOOST_CHECK_EQUAL(*(w->dim_w[2].transformed), w->dim_p[1]);
    BOOST_CHECK_EQUAL(*((w->dim_p[1].transformed)->transformed), w->dim_p[1]);

    underling_prepare_link(w, 2, 0);
    BOOST_CHECK_EQUAL(*(w->dim_p[2].transformed), w->dim_w[0]);
    BOOST_CHECK_EQUAL(*(w->dim_w[0].transformed), w->dim_p[2]);
    BOOST_CHECK_EQUAL(*((w->dim_p[2].transformed)->transformed), w->dim_p[2]);

    const underling_state state[ndim] = {
        underling_state_physical,
        underling_state_physical,
        underling_state_nottransformed
    };
    underling_prepare_state(w, state);
    for (int i = 0; i < w->ndim; ++i) {
        BOOST_CHECK_EQUAL(w->state[i], state[i]);
    }

    underling_workspace_free(w);
}
