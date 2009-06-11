#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <assert.h>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <boost/foreach.hpp>
#include <log4cxx/logger.h>

#include <boost/test/included/unit_test.hpp>

using namespace log4cxx;

LoggerPtr logger = Logger::getRootLogger();

struct bspline_deriv_params {
    size_t i_index;
    size_t i_derivorder;
    size_t j_index;
    size_t j_derivorder;

    gsl_matrix *dB;
    gsl_bspline_workspace *w;
    gsl_bspline_deriv_workspace *dw;
};

double eval_bspline_product(double x, void *p) {
    struct bspline_deriv_params* const params = (struct bspline_deriv_params *) p;

    const size_t nderiv = std::max(params->i_derivorder, params->j_derivorder);

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
    if (params->i_index >= index_start && params->i_index <= index_end) {
        i_value = gsl_matrix_get(params->dB,
                                 params->i_index - index_start,
                                 params->i_derivorder);
    }

    double j_value = 0.0;
    if (params->j_index >= index_start && params->j_index <= index_end) {
        j_value = gsl_matrix_get(params->dB,
                                 params->j_index - index_start,
                                 params->j_derivorder);
    }

    return i_value * j_value;
}

BOOST_AUTO_TEST_CASE( main_test )
{

    const size_t nbreak = 8, order = 2;
    const double left_end = 0.0, right_end = 7.0;
    const double eval_points[] = { 1.0 / 8.0, 1.0 / 16.0, 1.0 / 32.0, 0 };

    gsl_bspline_workspace *bw = gsl_bspline_alloc(order, nbreak);
    gsl_bspline_deriv_workspace *bdw = gsl_bspline_deriv_alloc(order);
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
    gsl_matrix *dB = gsl_matrix_alloc(order, 2);

    struct bspline_deriv_params f_params;
    f_params.i_derivorder = 0;
    f_params.j_derivorder = 0;
    f_params.dB           = dB;
    f_params.w            = bw;
    f_params.dw           = bdw;
    const gsl_function f = { eval_bspline_product, &f_params };



    for (int i = 0; i < A->size1; ++i) {
        f_params.i_index = i;
        for (int j = 0; j < A->size2; ++j) {
            f_params.j_index = j;
            // NOP
        }
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

            // DEBUG
            f_params.i_index     = 0;
            f_params.j_index     = 1;
            LOG4CXX_DEBUG(logger,
                          "At x = " << eval_point
                          << " basis[" << f_params.i_index << "]*basis["
                          << f_params.j_index << "] = "
                          << GSL_FN_EVAL(&f, eval_point)
                         );

        }
    }

    gsl_matrix_free(dB);
    gsl_matrix_free(A);
    gsl_vector_free(v);
    gsl_bspline_deriv_free(bdw);
    gsl_bspline_free(bw);
}
