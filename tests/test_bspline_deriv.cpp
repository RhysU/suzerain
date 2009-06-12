#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <boost/foreach.hpp>
#include <log4cxx/logger.h>

#include <boost/test/included/unit_test.hpp>

using namespace log4cxx;

LoggerPtr logger = Logger::getRootLogger();

struct bspline_deriv_params {
    size_t i_ndx;
    size_t i_nderiv;
    size_t j_ndx;
    size_t j_nderiv;

    gsl_matrix *dB;
    gsl_bspline_workspace *w;
    gsl_bspline_deriv_workspace *dw;
};

double eval_bspline_product(double x, void *p) {
    struct bspline_deriv_params* const params = (struct bspline_deriv_params *) p;

    const size_t nderiv = std::max(params->i_nderiv, params->j_nderiv);

    assert(params->dB->size1 >= gsl_bspline_order(params->w));
    assert(params->dB->size2 >= nderiv + 1);

    size_t index_start, index_end;
    gsl_bspline_deriv_eval_nonzero(x,
                                   nderiv,
                                   params->dB,
                                   &index_start,
                                   &index_end,
                                   params->w,
                                   params->dw);

    double i_value = 0.0;
    if (params->i_ndx >= index_start && params->i_ndx <= index_end) {
        i_value = gsl_matrix_get(params->dB,
                                 params->i_ndx - index_start,
                                 params->i_nderiv);
    }

    double j_value = 0.0;
    if (params->j_ndx >= index_start && params->j_ndx <= index_end) {
        j_value = gsl_matrix_get(params->dB,
                                 params->j_ndx - index_start,
                                 params->j_nderiv);
    }

    return i_value * j_value;
}

BOOST_AUTO_TEST_CASE( main_test )
{

    const size_t nbreak = 8, k = 2;
    const double left_end = 0.0, right_end = 7.0;
    const double eval_points[] = { 1.0 / 8.0, 1.0 / 16.0, 1.0 / 32.0, 0 };

    gsl_bspline_workspace *bw = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_deriv_workspace *bdw = gsl_bspline_deriv_alloc(k);
    gsl_bspline_knots_uniform(left_end, right_end, bw);

    if (logger->isDebugEnabled()) {
        for (size_t i = 0; i < bw->knots->size; i++) {
            LOG4CXX_DEBUG(logger,
                          "knot[" << i << "] at " << gsl_vector_get(bw->knots, i)
                         )
        }
    }

    gsl_matrix *A = gsl_matrix_alloc(gsl_bspline_ncoeffs(bw),
            gsl_bspline_ncoeffs(bw));
    gsl_matrix *dB = gsl_matrix_alloc(k, 2);

    struct bspline_deriv_params f_params;
    f_params.i_nderiv = 0;
    f_params.j_nderiv = 0;
    f_params.dB       = dB;
    f_params.w        = bw;
    f_params.dw       = bdw;
    const gsl_function f = { eval_bspline_product, &f_params };

    for (f_params.i_ndx = 0; f_params.i_ndx < A->size1; ++f_params.i_ndx) {
        for (f_params.j_ndx = 0; f_params.j_ndx < A->size1; ++f_params.j_ndx) {

            double result = 0.0;
            double abserr = 0.0;
            size_t neval  = 0;

            if (abs(f_params.i_ndx - f_params.j_ndx) < k) {
                const size_t left_index  = std::max(f_params.i_ndx,
                                                    f_params.j_ndx);
                const size_t right_index = std::min(f_params.i_ndx,
                                                    f_params.j_ndx) + k;
                const double left  = gsl_vector_get(bw->knots, left_index);
                const double right = gsl_vector_get(bw->knots, right_index);
                gsl_integration_qng(&f, left, right, 0.0, 0.5, // FIXME suspect
                                    &result, &abserr, &neval);
            }
            gsl_matrix_set(A, f_params.i_ndx, f_params.j_ndx, result);

            printf(" %10g", result);
        }
        printf("\n");
    }


    gsl_vector *v = gsl_vector_alloc(gsl_bspline_ncoeffs(bw));

    BOOST_FOREACH(double eval_point, eval_points) {
        gsl_bspline_eval(eval_point, v, bw);

        if (logger->isDebugEnabled()) {
            for (size_t i = 0; i < v->size; i++) {
                LOG4CXX_DEBUG(logger,
                              "At x = " << eval_point
                              << " basis[" << i << "] = "
                              << gsl_vector_get(v, i)
                             )
            }
        }
    }

    gsl_matrix_free(dB);
    gsl_matrix_free(A);
    gsl_vector_free(v);
    gsl_bspline_deriv_free(bdw);
    gsl_bspline_free(bw);
}
