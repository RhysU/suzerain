#define BOOST_TEST_MODULE $Id$

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

double f(double x, void *p) {
    double retval = 0.0;

    if (0.0 < x && x <= 1.0) {
        retval = x;
    } else if (1.0 < x && x <= 2.0) {
        retval = 1.0 - x;
    }
    retval *= retval;

    return retval;
}

int main(int argc, char **argv)
{
    double abserr = 0.0;
    double relerr = 1.0;
    double left = 0.0;
    double right = 2.0;

    const gsl_function gsl_f = { f, NULL };

    if (argc >= 2) abserr = atof(argv[1]);
    if (argc >= 3) relerr = atof(argv[2]);
    if (argc >= 4) left   = atof(argv[2]);
    if (argc >= 5) right  = atof(argv[2]);

    printf("abserr:           %g\n", abserr);
    printf("relerr:           %g\n", relerr);
    printf("left:             %g\n", left);
    printf("right:            %g\n", right);

    {
        double result1, result2;
        double esterr1, esterr2;
        size_t neval1, neval2;

        const double center = (left+right)/2.0;
        gsl_integration_qng(&gsl_f, left, center, abserr, relerr,
                            &result1, &esterr1, &neval1);
        gsl_integration_qng(&gsl_f, center, right, abserr, relerr,
                            &result2, &esterr2, &neval2);

        printf("\nTwo intervals...\n");
        printf("Result:           %g\n", result1 + result2);
        printf("Estimated error:  %g\n", esterr1 + esterr2);
        printf("# of evaluations: %zu\n", neval1 + neval2);
    }

    {
        double result;
        double esterr;
        size_t neval;

        gsl_integration_qng(&gsl_f, left, right, abserr, relerr,
                            &result, &esterr, &neval);
        printf("\nOne interval...\n");
        printf("Result:           %g\n", result);
        printf("Estimated error:  %g\n", esterr);
        printf("# of evaluations: %zu\n", neval);
    }

    return 0;
}

