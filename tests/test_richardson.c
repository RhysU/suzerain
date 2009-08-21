#include <suzerain/config.h>

#include <stdlib.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
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
        gsl_test(suzerain_richardson_extrapolation_step(Ah, Aht, 1, t),
                "Unexpected error reported in %s for t=%f", __func__, t);
        gsl_test_abs(gsl_vector_get(Ah, 0), 3.0, GSL_DBL_EPSILON,
                "Extrapolation in %s for ki=%d, t=%f", __func__, ki, t);
    }

    {
        const double t = 3.0;
        gsl_vector_set(Ah,  0, 1.0);
        gsl_vector_set(Aht, 0, 2.0);
        gsl_test(suzerain_richardson_extrapolation_step(Ah, Aht, 1, t),
                "Unexpected error reported in %s for t=%f", __func__, t);
        gsl_test_abs(gsl_vector_get(Ah, 0), 5.0/2.0, GSL_DBL_EPSILON,
                "Extrapolation in %s for ki=%d, t=%f", __func__, ki, t);
    }

    gsl_vector_free(Aht);
    gsl_vector_free(Ah);
}

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    test_richardson_extrapolation_step();

    exit(gsl_test_summary());
}
