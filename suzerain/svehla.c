#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_spline.h>
#include <gsl/gsl_const_mksa.h>
#include <suzerain/svehla.h>

// Temperature in degrees Kelvin from NASA TR R-132 page 117
static const double TableIV_T[] = {
     100,
     200,
     300,
     400,
     500,
     600,
     700,
     800,
     900,
    1000,
    //
    1100,
    1200,
    1300,
    1400,
    1500,
    1600,
    1700,
    1800,
    1900,
    2000,
    //
    2100,
    2200,
    2300,
    2400,
    2500,
    2600,
    2700,
    2800,
    2900,
    3000,
    //
    3100,
    3200,
    3300,
    3400,
    3500,
    3600,
    3700,
    3800,
    3900,
    4000,
    //
    4100,
    4200,
    4300,
    4400,
    4500,
    4600,
    4700,
    4800,
    4900,
    5000
};

static const double TableIV_Air_mu[sizeof(TableIV_T)/sizeof(TableIV_T[0])] = {
    73.8e6 * GSL_CONST_MKSA_POISE,
    136.0e6 * GSL_CONST_MKSA_POISE,
    185.2e6 * GSL_CONST_MKSA_POISE,
    227.2e6 * GSL_CONST_MKSA_POISE,
    264.7e6 * GSL_CONST_MKSA_POISE,
    299.2e6 * GSL_CONST_MKSA_POISE,
    331.3e6 * GSL_CONST_MKSA_POISE,
    361.4e6 * GSL_CONST_MKSA_POISE,
    389.8e6 * GSL_CONST_MKSA_POISE,
    417.1e6 * GSL_CONST_MKSA_POISE,
    //
    443.5e6 * GSL_CONST_MKSA_POISE,
    469.5e6 * GSL_CONST_MKSA_POISE,
    495.1e6 * GSL_CONST_MKSA_POISE,
    519.7e6 * GSL_CONST_MKSA_POISE,
    543.6e6 * GSL_CONST_MKSA_POISE,
    567.0e6 * GSL_CONST_MKSA_POISE,
    589.8e6 * GSL_CONST_MKSA_POISE,
    612.1e6 * GSL_CONST_MKSA_POISE,
    633.9e6 * GSL_CONST_MKSA_POISE,
    655.3e6 * GSL_CONST_MKSA_POISE,
    //
    676.3e6 * GSL_CONST_MKSA_POISE,
    697.0e6 * GSL_CONST_MKSA_POISE,
    717.3e6 * GSL_CONST_MKSA_POISE,
    737.3e6 * GSL_CONST_MKSA_POISE,
    757.0e6 * GSL_CONST_MKSA_POISE,
    776.5e6 * GSL_CONST_MKSA_POISE,
    795.6e6 * GSL_CONST_MKSA_POISE,
    814.5e6 * GSL_CONST_MKSA_POISE,
    833.2e6 * GSL_CONST_MKSA_POISE,
    851.6e6 * GSL_CONST_MKSA_POISE,
    //
    869.8e6 * GSL_CONST_MKSA_POISE,
    887.8e6 * GSL_CONST_MKSA_POISE,
    905.6e6 * GSL_CONST_MKSA_POISE,
    923.2e6 * GSL_CONST_MKSA_POISE,
    940.6e6 * GSL_CONST_MKSA_POISE,
    957.9e6 * GSL_CONST_MKSA_POISE,
    975.9e6 * GSL_CONST_MKSA_POISE,
    991.8e6 * GSL_CONST_MKSA_POISE,
    008.6e6 * GSL_CONST_MKSA_POISE,
    025.2e6 * GSL_CONST_MKSA_POISE,
    //
    041.6e6 * GSL_CONST_MKSA_POISE,
    058.0e6 * GSL_CONST_MKSA_POISE,
    074.1e6 * GSL_CONST_MKSA_POISE,
    090.2e6 * GSL_CONST_MKSA_POISE,
    106.1e6 * GSL_CONST_MKSA_POISE,
    121.9e6 * GSL_CONST_MKSA_POISE,
    137.5e6 * GSL_CONST_MKSA_POISE,
    153.1e6 * GSL_CONST_MKSA_POISE,
    168.5e6 * GSL_CONST_MKSA_POISE,
    183.8e6 * GSL_CONST_MKSA_POISE
};

gsl_spline * suzerain_svehla_air_mu_vs_T()
{
    // Viscosity in Pascal-seconds from NASA TR R-132 page 117
    static const size_t N = sizeof(TableIV_T)/sizeof(TableIV_T[0]);
    assert(N == sizeof(TableIV_Air_mu)/sizeof(TableIV_Air_mu[0]));
    gsl_spline * s = gsl_spline_alloc(gsl_interp_cspline, N);
    if (s) {
        if (gsl_spline_init(s, TableIV_T, TableIV_Air_mu, N)) {
            gsl_spline_free(s);
            s = NULL;
        }
    }
    return s;
}
