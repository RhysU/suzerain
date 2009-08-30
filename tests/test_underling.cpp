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

    /* Check information across each stage with all directions transformed */
    {
        const int ndim   = 3;
        const int nstage = ndim+1;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);

        const underling_state states[nstage][ndim] = {
            { UNDERLING_STATE_PHYSICAL,
                UNDERLING_STATE_PHYSICAL,
                UNDERLING_STATE_PHYSICAL },
            { UNDERLING_STATE_PHYSICAL,
                UNDERLING_STATE_PHYSICAL,
                UNDERLING_STATE_WAVE },
            { UNDERLING_STATE_PHYSICAL,
                UNDERLING_STATE_WAVE,
                UNDERLING_STATE_WAVE },
            { UNDERLING_STATE_WAVE,
                UNDERLING_STATE_WAVE,
                UNDERLING_STATE_WAVE }
        };
        for (int i = 0; i < nstage; ++i) {
            for (int j = 0; j < ndim; ++j) {
                BOOST_CHECK_EQUAL(w->stage[i].dim[j].state, states[i][j]);
            }
        }
        underling_workspace_free(w);
    }

    /* Check information across each stage with one untransformed direction */
    {
        const int ndim   = 3;
        const int nstage = ndim;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);

        const underling_state states[nstage][ndim] = {
            { UNDERLING_STATE_PHYSICAL,
                UNDERLING_STATE_PHYSICAL,
                UNDERLING_STATE_NOTTRANSFORMED },
            { UNDERLING_STATE_PHYSICAL,
                UNDERLING_STATE_NOTTRANSFORMED,
                UNDERLING_STATE_WAVE },
            { UNDERLING_STATE_NOTTRANSFORMED,
                UNDERLING_STATE_WAVE,
                UNDERLING_STATE_WAVE }
        };
        for (int i = 0; i < nstage; ++i) {
            for (int j = 0; j < ndim; ++j) {
                BOOST_CHECK_EQUAL(w->stage[i].dim[j].state, states[i][j]);
            }
        }
        underling_workspace_free(w);
    }
}
