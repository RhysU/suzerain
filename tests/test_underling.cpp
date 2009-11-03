#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>
#include <suzerain/common.hpp>
#include <suzerain/underling.h>
#include <boost/test/included/unit_test.hpp>

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

    /* Check information across each stage for low dimension problem */
    {
        const int ndim   = 1;
        const int nstage = ndim+1;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);

        /* Check states are correct in each stage */
        const underling_state states[nstage][ndim] = {
            { UNDERLING_STATE_PHYSICAL },
            { UNDERLING_STATE_WAVE }
        };
        for (int i = 0; i < nstage; ++i) {
            for (int j = 0; j < ndim; ++j) {
                BOOST_CHECK_EQUAL(w->stage[i].dim[j].state, states[i][j]);
            }
        }

        /* Check next_r2c pointers are correct */
        BOOST_CHECK_EQUAL(w->stage[0].dim[0].next_r2c, &(w->stage[1].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[0].next_r2c, (void *) NULL);

        /* Check next_c2r pointers are correct */
        BOOST_CHECK_EQUAL(w->stage[1].dim[0].next_c2r, &(w->stage[0].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[0].dim[0].next_c2r, (void *) NULL);

        underling_workspace_free(w);
    }

    /* Check information across each stage for untransformed low dimension */
    {
        const int ndim   = 1;
        const int nstage = 1;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);

        /* Check states are correct in each stage */
        const underling_state states[nstage][ndim] = {
            { UNDERLING_STATE_NOTTRANSFORMED }
        };
        for (int i = 0; i < nstage; ++i) {
            for (int j = 0; j < ndim; ++j) {
                BOOST_CHECK_EQUAL(w->stage[i].dim[j].state, states[i][j]);
            }
        }

        /* Check next_r2c pointers are correct */
        BOOST_CHECK_EQUAL(w->stage[0].dim[0].next_r2c, (void *) NULL);

        /* Check next_c2r pointers are correct */
        BOOST_CHECK_EQUAL(w->stage[0].dim[0].next_c2r, (void *) NULL);

        underling_workspace_free(w);
    }

    /* Check information across each stage with all directions transformed */
    {
        const int ndim   = 3;
        const int nstage = ndim+1;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);

        /* Check states are correct in each stage */
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

        /* Check next_r2c pointers are correct */
        BOOST_CHECK_EQUAL(w->stage[0].dim[0].next_r2c, &(w->stage[1].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[0].dim[1].next_r2c, &(w->stage[1].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[0].dim[2].next_r2c, &(w->stage[1].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[0].next_r2c, &(w->stage[2].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[1].next_r2c, &(w->stage[2].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[2].next_r2c, &(w->stage[2].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[0].next_r2c, &(w->stage[3].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[1].next_r2c, &(w->stage[3].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[2].next_r2c, &(w->stage[3].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[3].dim[0].next_r2c, (void *) NULL);
        BOOST_CHECK_EQUAL(w->stage[3].dim[1].next_r2c, (void *) NULL);
        BOOST_CHECK_EQUAL(w->stage[3].dim[2].next_r2c, (void *) NULL);

        /* Check next_c2r pointers are correct */
        BOOST_CHECK_EQUAL(w->stage[3].dim[0].next_c2r, &(w->stage[2].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[3].dim[1].next_c2r, &(w->stage[2].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[3].dim[2].next_c2r, &(w->stage[2].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[0].next_c2r, &(w->stage[1].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[1].next_c2r, &(w->stage[1].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[2].next_c2r, &(w->stage[1].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[0].next_c2r, &(w->stage[0].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[1].next_c2r, &(w->stage[0].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[2].next_c2r, &(w->stage[0].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[0].dim[0].next_c2r, (void *) NULL);
        BOOST_CHECK_EQUAL(w->stage[0].dim[1].next_c2r, (void *) NULL);
        BOOST_CHECK_EQUAL(w->stage[0].dim[2].next_c2r, (void *) NULL);
        underling_workspace_free(w);
    }

    /* Check information across each stage with one untransformed direction */
    {
        const int ndim   = 3;
        const int nstage = ndim;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);

        /* Check states are correct in each stage */
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

        /* Check next_r2c pointers are correct */
        BOOST_CHECK_EQUAL(w->stage[0].dim[0].next_r2c, &(w->stage[1].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[0].dim[1].next_r2c, &(w->stage[1].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[0].dim[2].next_r2c, &(w->stage[1].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[0].next_r2c, &(w->stage[2].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[1].next_r2c, &(w->stage[2].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[2].next_r2c, &(w->stage[2].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[0].next_r2c, (void *) NULL);
        BOOST_CHECK_EQUAL(w->stage[2].dim[1].next_r2c, (void *) NULL);
        BOOST_CHECK_EQUAL(w->stage[2].dim[2].next_r2c, (void *) NULL);

        /* Check next_c2r pointers are correct */
        BOOST_CHECK_EQUAL(w->stage[2].dim[0].next_c2r, &(w->stage[1].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[1].next_c2r, &(w->stage[1].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[2].dim[2].next_c2r, &(w->stage[1].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[0].next_c2r, &(w->stage[0].dim[1]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[1].next_c2r, &(w->stage[0].dim[2]));
        BOOST_CHECK_EQUAL(w->stage[1].dim[2].next_c2r, &(w->stage[0].dim[0]));
        BOOST_CHECK_EQUAL(w->stage[0].dim[0].next_c2r, (void *) NULL);
        BOOST_CHECK_EQUAL(w->stage[0].dim[1].next_c2r, (void *) NULL);
        BOOST_CHECK_EQUAL(w->stage[0].dim[2].next_c2r, (void *) NULL);

        underling_workspace_free(w);
    }
}

BOOST_AUTO_TEST_CASE( name_dimension )
{
    /* Check dimension names across each stage for low dimension problem */
    {
        const int ndim   = 1;
        const int nstage = ndim+1;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);
        underling_name_dimension(w, 0, 0, "X");

        BOOST_CHECK_EQUAL(w->stage[0].dim[0].name, "X");
        BOOST_CHECK_EQUAL(w->stage[1].dim[0].name, "X");

        underling_workspace_free(w);
    }

    /* Check dimension names across each stage specified using last stage */
    {
        const int ndim   = 3;
        const int nstage = ndim;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);
        underling_name_dimension(w, 0, 2, "Y");
        underling_name_dimension(w, 1, 2, "X");
        underling_name_dimension(w, 2, 2, "Z");

        BOOST_CHECK_EQUAL(w->stage[2].dim[0].name, "Y");
        BOOST_CHECK_EQUAL(w->stage[2].dim[1].name, "X");
        BOOST_CHECK_EQUAL(w->stage[2].dim[2].name, "Z");

        BOOST_CHECK_EQUAL(w->stage[1].dim[0].name, "Z");
        BOOST_CHECK_EQUAL(w->stage[1].dim[1].name, "Y");
        BOOST_CHECK_EQUAL(w->stage[1].dim[2].name, "X");

        BOOST_CHECK_EQUAL(w->stage[0].dim[0].name, "X");
        BOOST_CHECK_EQUAL(w->stage[0].dim[1].name, "Z");
        BOOST_CHECK_EQUAL(w->stage[0].dim[2].name, "Y");

        underling_workspace_free(w);
    }

    /* Check dimension names across each stage specified using intermediate stage */
    {
        const int ndim   = 3;
        const int nstage = ndim;
        underling_workspace * const w
            = underling_workspace_alloc(ndim, nstage);
        underling_name_dimension(w, 0, 1, "Z");
        underling_name_dimension(w, 1, 1, "Y");
        underling_name_dimension(w, 2, 1, "X");

        BOOST_CHECK_EQUAL(w->stage[2].dim[0].name, "Y");
        BOOST_CHECK_EQUAL(w->stage[2].dim[1].name, "X");
        BOOST_CHECK_EQUAL(w->stage[2].dim[2].name, "Z");

        BOOST_CHECK_EQUAL(w->stage[1].dim[0].name, "Z");
        BOOST_CHECK_EQUAL(w->stage[1].dim[1].name, "Y");
        BOOST_CHECK_EQUAL(w->stage[1].dim[2].name, "X");

        BOOST_CHECK_EQUAL(w->stage[0].dim[0].name, "X");
        BOOST_CHECK_EQUAL(w->stage[0].dim[1].name, "Z");
        BOOST_CHECK_EQUAL(w->stage[0].dim[2].name, "Y");

        underling_workspace_free(w);
    }
}
