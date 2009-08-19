#include <suzerain/config.h>

#include <stdlib.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <suzerain/timestepper.h>

int
main(int argc, char **argv)
{
    gsl_ieee_env_setup();

    exit(gsl_test_summary());
}
