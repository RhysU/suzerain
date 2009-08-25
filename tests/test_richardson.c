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

        gsl_vector * k = gsl_vector_alloc(1);
        gsl_vector_set(k, 0, 1);
        gsl_test(suzerain_richardson_extrapolation(data2, 2, k, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_vector_free(k);
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
    gsl_vector * k1 = gsl_vector_alloc(1);
    gsl_vector * k2 = gsl_vector_alloc(2);

    {
        gsl_matrix_set(data2, 0, 0, 1.0);
        gsl_matrix_set(data2, 0, 1, 2.0);
        gsl_vector_set(k1, 0, 1);

        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data2, 0, 0), 3.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);
    }

    {
        gsl_matrix_set(data2, 0, 0, 2.0);
        gsl_matrix_set(data2, 0, 1, 3.0);
        gsl_vector_set(k1, 0, 1);

        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data2, 0, 0), 4.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);
    }

    {
        gsl_matrix_set(data2, 0, 0, 3.0);
        gsl_matrix_set(data2, 0, 1, 4.0);
        gsl_vector_set(k1, 0, 2);

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
        gsl_vector_set(k1, 0, 1);

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
        gsl_vector_set(k2, 0, 1);
        gsl_vector_set(k2, 1, 2);

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

        gsl_test_abs(gsl_matrix_get(normtable, 0, 0), 1.0, GSL_DBL_EPSILON,
                "%s normtable (0,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 1, 0), 2.0, GSL_DBL_EPSILON,
                "%s normtable (1,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 2, 0), 3.0, GSL_DBL_EPSILON,
                "%s normtable (2,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 1, 1), 3.0, GSL_DBL_EPSILON,
                "%s normtable (1,1) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 2, 1), 4.0, GSL_DBL_EPSILON,
                "%s normtable (2,1) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 2, 2), 13.0/3.0, GSL_DBL_EPSILON,
                "%s normtable (2,2) value", __func__);

        gsl_matrix_free(normtable);
    }

    {
        gsl_matrix_set(data3, 0, 0, 1.0);
        gsl_matrix_set(data3, 0, 1, 2.0);
        gsl_matrix_set(data3, 0, 2, 3.0);
        gsl_vector_set(k2, 0, 2);
        gsl_vector_set(k2, 1, 3);

        gsl_test(suzerain_richardson_extrapolation(data3, 2, k2, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data3, 0, 0), 73.0/21.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);
    }

    {
        gsl_matrix_set(data3, 0, 0, 1.0);
        gsl_matrix_set(data3, 0, 1, 2.0);
        gsl_matrix_set(data3, 0, 2, 3.0);
        gsl_vector_set(k1, 0, 2);

        gsl_test(suzerain_richardson_extrapolation(data3, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data3, 0, 0), 73.0/21.0, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);
    }


    gsl_vector_free(k2);
    gsl_vector_free(k1);
    gsl_matrix_free(data3);
    gsl_matrix_free(data2);

}

void
test_richardson_extrapolation_multiplelevels()
{
    gsl_matrix * data2 = gsl_matrix_alloc(1, 2);
    gsl_matrix * data4 = gsl_matrix_alloc(1, 4);
    gsl_vector * k1 = gsl_vector_alloc(1);
    gsl_vector * k2 = gsl_vector_alloc(2);
    gsl_vector * k3 = gsl_vector_alloc(3);

    {
        gsl_vector_set(k1, 0, 3);

        gsl_matrix_set(data2, 0, 0, 1.0);
        gsl_matrix_set(data2, 0, 1, 2.0);
        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        const double xx = gsl_matrix_get(data2, 0, 0);

        gsl_matrix_set(data2, 0, 0, 2.0);
        gsl_matrix_set(data2, 0, 1, 3.0);
        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        const double xy = gsl_matrix_get(data2, 0, 0);

        gsl_matrix_set(data2, 0, 0, 3.0);
        gsl_matrix_set(data2, 0, 1, 4.0);
        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        const double xz = gsl_matrix_get(data2, 0, 0);

        gsl_vector_set(k1, 0, 6);

        gsl_matrix_set(data2, 0, 0, xx);
        gsl_matrix_set(data2, 0, 1, xy);
        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        const double yx = gsl_matrix_get(data2, 0, 0);

        gsl_matrix_set(data2, 0, 0, xy);
        gsl_matrix_set(data2, 0, 1, xz);
        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        const double yy = gsl_matrix_get(data2, 0, 0);

        gsl_vector_set(k1, 0, 9);

        gsl_matrix_set(data2, 0, 0, yx);
        gsl_matrix_set(data2, 0, 1, yy);
        gsl_test(suzerain_richardson_extrapolation(data2, 2, k1, NULL, NULL),
                "Unexpected error reported in %s");
        const double zx = gsl_matrix_get(data2, 0, 0);

        gsl_vector_set(k3, 0, 3);
        gsl_vector_set(k3, 1, 6);
        gsl_vector_set(k3, 2, 9);

        gsl_matrix_set(data4, 0, 0, 1.0);
        gsl_matrix_set(data4, 0, 1, 2.0);
        gsl_matrix_set(data4, 0, 2, 3.0);
        gsl_matrix_set(data4, 0, 3, 4.0);
        gsl_test(suzerain_richardson_extrapolation(data4, 2, k3, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data4, 0, 0), zx, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);

        gsl_matrix * normtable = gsl_matrix_alloc(data4->size2, data4->size2);

        gsl_vector_set(k2, 0, 3);
        gsl_vector_set(k2, 1, 6);

        gsl_matrix_set(data4, 0, 0, 1.0);
        gsl_matrix_set(data4, 0, 1, 2.0);
        gsl_matrix_set(data4, 0, 2, 3.0);
        gsl_matrix_set(data4, 0, 3, 4.0);
        gsl_test(suzerain_richardson_extrapolation(data4, 2, k2, normtable, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data4, 0, 0), zx, GSL_DBL_EPSILON,
                "%s scalar correct result at %s:%d",
                __func__, __FILE__, __LINE__);

        for (int i = 0; i < normtable->size1 - 1; ++i) {
            for (int j = i+1; j < normtable->size2; ++j) {
                gsl_test(!gsl_isnan(gsl_matrix_get(normtable, i, j)),
                        "%s expected NAN in normtable at (%d, %d)",
                        __func__, i, j);
            }
        }

        gsl_test_abs(gsl_matrix_get(normtable, 0, 0), 1.0, GSL_DBL_EPSILON,
                "%s normtable (0,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 1, 0), 2.0, GSL_DBL_EPSILON,
                "%s normtable (1,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 2, 0), 3.0, GSL_DBL_EPSILON,
                "%s normtable (2,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 3, 0), 4.0, GSL_DBL_EPSILON,
                "%s normtable (3,0) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 1, 1), xx, GSL_DBL_EPSILON,
                "%s normtable (1,1) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 2, 1), xy, GSL_DBL_EPSILON,
                "%s normtable (2,1) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 3, 1), xz, GSL_DBL_EPSILON,
                "%s normtable (3,1) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 2, 2), yx, GSL_DBL_EPSILON,
                "%s normtable (2,2) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 3, 2), yy, GSL_DBL_EPSILON,
                "%s normtable (3,2) value", __func__);
        gsl_test_abs(gsl_matrix_get(normtable, 3, 3), zx, GSL_DBL_EPSILON,
                "%s normtable (3,3) value", __func__);

        gsl_matrix_free(normtable);
    }

    gsl_vector_free(k3);
    gsl_vector_free(k2);
    gsl_vector_free(k1);
    gsl_matrix_free(data4);
    gsl_matrix_free(data2);
}

void
test_richardson_extrapolation_vectorinputs()
{
    {
        gsl_matrix * data = gsl_matrix_alloc(2,2);
        gsl_matrix_set(data, 0, 0, 1.0);
        gsl_matrix_set(data, 0, 1, 2.0);
        gsl_matrix_set(data, 1, 0, 2.0);
        gsl_matrix_set(data, 1, 1, 3.0);

        gsl_matrix * normtable = gsl_matrix_alloc(data->size2, data->size2);

        gsl_vector * exact = gsl_vector_alloc(data->size2);
        gsl_vector_set(exact, 0, 3.0);
        gsl_vector_set(exact, 1, 4.0);

        gsl_test(suzerain_richardson_extrapolation(
                    data, 2, NULL, normtable, exact),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data, 0, 0), 3.0, GSL_DBL_EPSILON,
                "%s data (0,0) value at %s:%d", __func__, __FILE__, __LINE__);
        gsl_test_abs(gsl_matrix_get(data, 1, 0), 4.0, GSL_DBL_EPSILON,
                "%s data (1,0) value at %s:%d", __func__, __FILE__, __LINE__);

        gsl_test_abs(gsl_matrix_get(normtable, 1, 1), 0.0, GSL_DBL_EPSILON,
                "%s normtable (1,1) value at %s:%d",
                __func__, __FILE__, __LINE__);

        gsl_vector_free(exact);
        gsl_matrix_free(normtable);
        gsl_matrix_free(data);
    }

    {
        gsl_matrix * data = gsl_matrix_alloc(2,3);
        gsl_matrix_set(data, 0, 0, 1.0);
        gsl_matrix_set(data, 0, 1, 2.0);
        gsl_matrix_set(data, 0, 2, 3.0);
        gsl_matrix_set(data, 1, 0, 2.0);
        gsl_matrix_set(data, 1, 1, 3.0);
        gsl_matrix_set(data, 1, 2, 4.0);

        gsl_test(suzerain_richardson_extrapolation(data, 2, NULL, NULL, NULL),
                "Unexpected error reported in %s");
        gsl_test_abs(gsl_matrix_get(data, 0, 0), 13.0/3.0, GSL_DBL_EPSILON,
                "%s data (0,0) value at %s:%d", __func__, __FILE__, __LINE__);
        gsl_test_abs(gsl_matrix_get(data, 1, 0), 16.0/3.0, GSL_DBL_EPSILON,
                "%s data (1,0) value at %s:%d", __func__, __FILE__, __LINE__);

        gsl_matrix_free(data);
    }
}

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    test_richardson_extrapolation_step();
    test_richardson_extrapolation_defaults();
    test_richardson_extrapolation_twolevels();
    test_richardson_extrapolation_multiplelevels();
    test_richardson_extrapolation_vectorinputs();

    exit(gsl_test_summary());
}
