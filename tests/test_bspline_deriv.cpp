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

BOOST_AUTO_TEST_CASE( main_test )
{
/*
    const size_t k = 3;
    const double bpoint_data[]    = { 0.0, 0.2, 0.5, 0.75, 1.0 };
*/
/*
    const size_t k = 3;
    const double bpoint_data[]    = { 0.0, 2.0, 4.0 };
*/
    const size_t k = 4;
    const double bpoint_data[]    = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };

    const size_t nbreak = sizeof(bpoint_data)/sizeof(bpoint_data[0]);
    gsl_vector_const_view bpoints = gsl_vector_const_view_array(bpoint_data, nbreak);
    gsl_bspline_workspace *bw = gsl_bspline_alloc(k, nbreak);
    gsl_bspline_deriv_workspace *bdw = gsl_bspline_deriv_alloc(k);
    gsl_bspline_knots((const gsl_vector *) &bpoints, bw);
    const size_t ncoeff = gsl_bspline_ncoeffs(bw);

    if (logger->isDebugEnabled()) {
        for (size_t i = 0; i < bw->knots->size; ++i) {
            LOG4CXX_DEBUG(logger, "knot[" << i << "] = "
                                  << gsl_vector_get(bw->knots, i));
        }
    }

    if (logger->isDebugEnabled()) {
        for (size_t i = 0; i < ncoeff; ++i) {
            LOG4CXX_DEBUG(logger, "abscissa[" << i << "] = "
                                  << gsl_bspline_greville_abscissa(i, bw));
        }
    }

    const size_t nderiv = 2;
    gsl_matrix *M  = gsl_matrix_alloc(ncoeff, ncoeff);
    gsl_matrix *D1 = gsl_matrix_alloc(ncoeff, ncoeff);
    gsl_matrix *D2 = gsl_matrix_alloc(ncoeff, ncoeff);
    gsl_matrix *dB = gsl_matrix_alloc(ncoeff, nderiv+1);

    for (size_t j = 0; j < M->size1; ++j) {
        double x = gsl_bspline_greville_abscissa(j, bw);
        gsl_bspline_deriv_eval(x, nderiv, dB, bw, bdw);
        gsl_matrix_set_col(M,  j, (const gsl_vector *) &gsl_matrix_column(dB, 0));
        gsl_matrix_set_col(D1, j, (const gsl_vector *) &gsl_matrix_column(dB, 1));
        gsl_matrix_set_col(D2, j, (const gsl_vector *) &gsl_matrix_column(dB, 2));
    }

    for (size_t i = 0; i < M->size1; ++i) {
        for (size_t j = 0; j < M->size2; ++j) {
            printf(" %10g", gsl_matrix_get(M, i, j));
        }
        printf("\n");
    }

    gsl_matrix_free(dB);
    gsl_matrix_free(M);
    gsl_matrix_free(D1);
    gsl_matrix_free(D2);
    gsl_bspline_deriv_free(bdw);
    gsl_bspline_free(bw);
}
