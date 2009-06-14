#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <log4cxx/logger.h>

#include <boost/test/included/unit_test.hpp>

using namespace log4cxx;

LoggerPtr logger = Logger::getRootLogger();

int
gsl_bspline_basis_product(const double x,
                          const size_t i_index,
                          const size_t i_nderiv,
                          const size_t j_index,
                          const size_t j_nderiv,
                          double * result,
                          gsl_matrix * dB,
                          gsl_bspline_workspace * w,
                          gsl_bspline_deriv_workspace * dw) {

    const size_t nderiv = (i_nderiv > j_nderiv) ? i_nderiv : j_nderiv;
    size_t nonzero_start, nonzero_end;
    double i_value = 0.0, j_value = 0.0;

    if (dB->size1 < gsl_bspline_order(w)) {
        GSL_ERROR("dB->size1 too small for given B-spline order",
                  GSL_EINVAL);
    }

    if (dB->size2 < nderiv + 1) {
        GSL_ERROR("dB->size2 too small for number of requested derivatives",
                  GSL_EINVAL);
    }

    gsl_bspline_deriv_eval_nonzero(x, nderiv, dB,
                                   &nonzero_start, &nonzero_end,
                                   w, dw);

    if (i_index >= nonzero_start && i_index <= nonzero_end) {
        i_value = gsl_matrix_get(dB, i_index - nonzero_start, i_nderiv);
    }

    if (j_index >= nonzero_start && j_index <= nonzero_end) {
        j_value = gsl_matrix_get(dB, j_index - nonzero_start, j_nderiv);
    }

    *result = i_value * j_value;

    return GSL_SUCCESS;
}

int
gsl_bspline_greville_abscissae(gsl_bspline_workspace *w, gsl_vector *abscissae)
{
    const size_t k      = gsl_bspline_order(w);
    const size_t n      = gsl_bspline_ncoeffs(w);
    const size_t stride = abscissae->stride;
    double * data       = w->knots->data;
    size_t i;

    if (abscissae->size < n) {
        GSL_ERROR("abscissae->size too small for number of basis functions",
                  GSL_EINVAL);
    }

    for (i = 0; i < n; ++i) {
        gsl_vector_set(abscissae, i, gsl_stats_mean(data, stride, k));
        data += stride;
    }

    return GSL_SUCCESS;
}

BOOST_AUTO_TEST_CASE( main_test )
{
/*
    // Greville abscissa example taken from R. W. Johnson paper
    // Applied Numerical Mathematics vol 52 (2005) pages 63--75
    const size_t k = 4;
    const double bpoint_data[]    = { 0.0, 0.2, 0.5, 0.75, 1.0 };
*/
    const size_t k = 1;
    const double bpoint_data[]    = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    const size_t nbreak           = sizeof(bpoint_data)/sizeof(bpoint_data[0]);
    gsl_vector_const_view bpoints = gsl_vector_const_view_array(bpoint_data, nbreak);

    gsl_bspline_workspace *bw = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_deriv_workspace *bdw = gsl_bspline_deriv_alloc(k);
    gsl_bspline_knots((const gsl_vector *) &bpoints, bw);

    if (logger->isDebugEnabled()) {
        for (size_t i = 0; i < bw->knots->size; ++i) {
            LOG4CXX_DEBUG(logger, "knot[" << i << "] at " 
                                  << gsl_vector_get(bw->knots, i));
        }
    }

    gsl_vector *abscissae = gsl_vector_alloc(gsl_bspline_ncoeffs(bw));
    gsl_bspline_greville_abscissae(bw, abscissae);
    if (logger->isDebugEnabled()) {
        for (size_t i = 0; i < abscissae->size; ++i) {
            LOG4CXX_DEBUG(logger, "abscissa[" << i << "] at " 
                                  << gsl_vector_get(abscissae, i));
        }
    }

    gsl_matrix *A = gsl_matrix_alloc(gsl_bspline_ncoeffs(bw), gsl_bspline_ncoeffs(bw));
    gsl_matrix *dB = gsl_matrix_alloc(k, 2);

    for (size_t i = 0; i < A->size1; ++i) {
        double x = gsl_vector_get(abscissae, i);
        for (size_t j = 0; j < A->size2; ++j) {
            double result = 0.0;
            if (abs(i - j) < k) {
                gsl_bspline_basis_product(x, i, 0, j, 0, &result, dB, bw, bdw);
            }
            gsl_matrix_set(A, i, j, result);

            printf(" %10g", result);
        }
        printf("\n");
    }

    gsl_vector_free(abscissae);
    gsl_matrix_free(dB);
    gsl_matrix_free(A);
    gsl_bspline_deriv_free(bdw);
    gsl_bspline_free(bw);
}

