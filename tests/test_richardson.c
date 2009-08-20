#include <suzerain/config.h>

#include <stdlib.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>
#include <suzerain/richardson.h>

void
test_richardson_step()
{
    const int    ki           = 1;
    const size_t n            = 1;
    gsl_vector *Aih   = gsl_vector_alloc(n);
    gsl_vector *Aiht  = gsl_vector_alloc(n);
    gsl_vector *Aip1h = gsl_vector_alloc(n);

    gsl_vector_set(Aih,  0, 1.0);
    gsl_vector_set(Aiht, 0, 2.0);

    {
        const double t = 2.0;
        gsl_test(suzerain_richardson_step(Aih, Aiht, 1, t, Aip1h),
                "Unexpected error reported in %s for t=%f", __func__, t);
        gsl_test_abs(gsl_vector_get(Aip1h, 0), 3.0, GSL_DBL_EPSILON,
                "Extrapolation in %s for ki=%d, t=%f", __func__, ki, t);
    }

    {
        const double t = 3.0;
        gsl_test(suzerain_richardson_step(Aih, Aiht, 1, t, Aip1h),
                "Unexpected error reported in %s for t=%f", __func__, t);
        gsl_test_abs(gsl_vector_get(Aip1h, 0), 5.0/2.0, GSL_DBL_EPSILON,
                "Extrapolation in %s for ki=%d, t=%f", __func__, ki, t);
    }

    gsl_vector_free(Aip1h);
    gsl_vector_free(Aiht);
    gsl_vector_free(Aih);
}

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    test_richardson_step();

    exit(gsl_test_summary());
}
