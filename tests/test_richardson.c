#include <suzerain/config.h>

#include <stdlib.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_test.h>
#include <suzerain/richardson.h>

void
test_richardson_extrapolation_step()
{
    const int    ki           = 1;
    const size_t n            = 1;
    gsl_vector *Ah   = gsl_vector_alloc(n);
    gsl_vector *Aht  = gsl_vector_alloc(n);

    {
        const double t = 2.0;
        gsl_vector_set(Ah,  0, 1.0);
        gsl_vector_set(Aht, 0, 2.0);
        gsl_test(suzerain_richardson_extrapolation_step(Ah, Aht, t, ki),
                "Unexpected error reported in %s for t=%f", __func__, t);
        gsl_test_abs(gsl_vector_get(Ah, 0), 3.0, GSL_DBL_EPSILON,
                "Extrapolation in %s for ki=%d, t=%f", __func__, ki, t);
    }

    {
        const double t = 3.0;
        gsl_vector_set(Ah,  0, 1.0);
        gsl_vector_set(Aht, 0, 2.0);
        gsl_test(suzerain_richardson_extrapolation_step(Ah, Aht, t, ki),
                "Unexpected error reported in %s for t=%f", __func__, t);
        gsl_test_abs(gsl_vector_get(Ah, 0), 5.0/2.0, GSL_DBL_EPSILON,
                "Extrapolation in %s for ki=%d, t=%f", __func__, ki, t);
    }

    gsl_vector_free(Aht);
    gsl_vector_free(Ah);
}

void
test_richardson_extrapolation_defaults()
{
    gsl_matrix * data1 = gsl_matrix_alloc(1, 2);
    gsl_matrix * data2 = gsl_matrix_alloc(data1->size1, data1->size2);
    {
        gsl_matrix_set(data1, 0, 0, 1.0);
        gsl_matrix_set(data1, 0, 1, 2.0);

        gsl_test(suzerain_richardson_extrapolation(data1, 2, NULL, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data1, 0, 0), 3.0, GSL_DBL_EPSILON,
                "%s scalar correct result", __func__);
    }

    {
        gsl_matrix_set(data1, 0, 0, 1.0);
        gsl_matrix_set(data1, 0, 1, 2.0);

        gsl_test(suzerain_richardson_extrapolation(data1, 3, NULL, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data1, 0, 0), 5.0/2.0, GSL_DBL_EPSILON,
                "%s scalar correct result", __func__);
    }

    {
        gsl_matrix_set(data1, 0, 0, 1.0);
        gsl_matrix_set(data1, 0, 1, 2.0);
        gsl_matrix_memcpy(data2, data1);

        gsl_test(suzerain_richardson_extrapolation(data1, 2, NULL, NULL, NULL),
                "Unexpected error reported in %s");
        const double result1 = gsl_matrix_get(data1, 0, 0);

        gsl_test_abs(result1, 3.0, GSL_DBL_EPSILON,
                "%s scalar correct result", __func__);

        gsl_vector_int * k = gsl_vector_int_alloc(1);
        gsl_vector_int_set(k, 0, 1);
        gsl_test(suzerain_richardson_extrapolation(data2, 2, k, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_vector_int_free(k);
        const double result2 = gsl_matrix_get(data2, 0, 0);

        gsl_test_abs(result1, result2, GSL_DBL_EPSILON,
                "%s with k == NULL", __func__);
    }

    {
        gsl_matrix_set(data1, 0, 0, 1.0);
        gsl_matrix_set(data1, 0, 1, 2.0);
        gsl_matrix_memcpy(data2, data1);

        gsl_matrix * normtable1 = gsl_matrix_alloc(data1->size2, data1->size2);
        gsl_matrix * normtable2 = gsl_matrix_alloc(data1->size2, data1->size2);

        gsl_test(suzerain_richardson_extrapolation(
                    data1, 2, NULL, normtable1, NULL),
                "Unexpected error reported in %s");
        const double result1 = gsl_matrix_get(data1, 0, 0);

        gsl_test_abs(result1, 3.0, GSL_DBL_EPSILON,
                "%s correct extrapolation answer");

        for (int i = 0; i < normtable1->size1 - 1; ++i) {
            for (int j = i+1; j < normtable1->size2; ++j) {
                gsl_test(!gsl_isnan(gsl_matrix_get(normtable1, i, j)),
                        "%s expected NAN in normtable1 at (%d, %d)",
                        __func__, i, j);
            }
        }

        gsl_test_abs(gsl_matrix_get(normtable1, 0, 0), 1.0, GSL_DBL_EPSILON,
                "%s normtable1 (0,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable1, 1, 0), 2.0, GSL_DBL_EPSILON,
                "%s normtable1 (1,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable1, 1, 1), 3.0, GSL_DBL_EPSILON,
                "%s normtable1 (1,1) value", __func__);

        gsl_vector * exact = gsl_vector_alloc(1);
        gsl_vector_set_zero(exact);
        gsl_test(suzerain_richardson_extrapolation(
                    data2, 2.0, NULL, normtable2, exact),
                "Unexpected error reported in %s");
        gsl_vector_free(exact);
        const double result2 = gsl_matrix_get(data2, 0, 0);

        for (int i = 0; i < normtable2->size1 - 1; ++i) {
            for (int j = i+1; j < normtable2->size2; ++j) {
                gsl_test(!gsl_isnan(gsl_matrix_get(normtable2, i, j)),
                        "%s expected NAN in normtable2 at (%d, %d)",
                        __func__, i, j);
            }
        }

        gsl_test_abs(result1, result2, GSL_DBL_EPSILON,
                "%s with exact == NULL", __func__);

        for (int i = 0; i < normtable1->size1; ++i) {
            for (int j = 0; j < i + 1; ++j) {
                gsl_test_abs(gsl_matrix_get(normtable1, i, j),
                             gsl_matrix_get(normtable2, i, j),
                             GSL_DBL_EPSILON,
                             "%s identical normtable results at (%d, %d)",
                             __func__, i, j);
            }
        }

        gsl_matrix_free(normtable2);
        gsl_matrix_free(normtable1);
    }

    gsl_matrix_free(data2);
    gsl_matrix_free(data1);
}

void
test_richardson_extrapolation_twolevels()
{
    gsl_matrix * data2 = gsl_matrix_alloc(1, 2);
    gsl_matrix * data3 = gsl_matrix_alloc(1, 3);
    gsl_vector_int * k1 = gsl_vector_int_alloc(1);
    gsl_vector_int * k2 = gsl_vector_int_alloc(2);

    {
        gsl_matrix_set(data2, 0, 0, 1.0);
        gsl_matrix_set(data2, 0, 1, 2.0);
        gsl_vector_int_set(k1, 0, 1);

        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data2, 0, 0), 3.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);
    }

    {
        gsl_matrix_set(data2, 0, 0, 2.0);
        gsl_matrix_set(data2, 0, 1, 3.0);
        gsl_vector_int_set(k1, 0, 1);

        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data2, 0, 0), 4.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);
    }

    {
        gsl_matrix_set(data2, 0, 0, 3.0);
        gsl_matrix_set(data2, 0, 1, 4.0);
        gsl_vector_int_set(k1, 0, 2);

        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data2, 0, 0), 13.0/3.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);
    }

    {
        gsl_matrix_set(data3, 0, 0, 1.0);
        gsl_matrix_set(data3, 0, 1, 2.0);
        gsl_matrix_set(data3, 0, 2, 3.0);
        gsl_vector_int_set(k1, 0, 1);

        gsl_test(suzerain_richardson_extrapolation(data3, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data3, 0, 0), 13.0/3.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);
    }

    {
        gsl_matrix * normtable = gsl_matrix_alloc(data3->size2, data3->size2);

        gsl_matrix_set(data3, 0, 0, 1.0);
        gsl_matrix_set(data3, 0, 1, 2.0);
        gsl_matrix_set(data3, 0, 2, 3.0);
        gsl_vector_int_set(k2, 0, 1);
        gsl_vector_int_set(k2, 1, 2);

        gsl_test(suzerain_richardson_extrapolation(data3, 2, k2, normtable, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data3, 0, 0), 13.0/3.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);

        for (int i = 0; i < normtable->size1 - 1; ++i) {
            for (int j = i+1; j < normtable->size2; ++j) {
                gsl_test(!gsl_isnan(gsl_matrix_get(normtable, i, j)),
                        "%s expected NAN in normtable at (%d, %d)",
                        __func__, i, j);
            }
        }

        gsl_matrix_free(normtable);
    }

    gsl_vector_int_free(k2);
    gsl_vector_int_free(k1);
    gsl_matrix_free(data3);
    gsl_matrix_free(data2);

}

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    test_richardson_extrapolation_step();
    test_richardson_extrapolation_defaults();
    test_richardson_extrapolation_twolevels();

    exit(gsl_test_summary());
}
