/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012, 2013 Rhys Ulerich
 * Copyright (C) 2012, 2013 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

/** @file
 * @copydoc blasius.c
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/blasius.h>

#include <suzerain/common.h>

// Variable $\eta$ from Table3b in http://arxiv.org/format/1006.3888v1
static const double Table3b_eta[] = {
    0.0E+00,
    2.0E-01,
    4.0E-01,
    6.0E-01,
    8.0E-01,
    1.0E+00,
    1.2E+00,
    1.4E+00,
    1.6E+00,
    1.8E+00,
    2.0E+00,
    2.2E+00,
    2.4E+00,
    2.6E+00,
    2.8E+00,
    3.0E+00,
    3.2E+00,
    3.4E+00,
    3.6E+00,
    3.8E+00,
    4.0E+00,
    4.2E+00,
    4.4E+00,
    4.6E+00,
    4.8E+00,
    5.0E+00,
    5.2E+00,
    5.4E+00,
    5.6E+00,
    5.8E+00,
    6.0E+00,
    6.2E+00,
    6.4E+00,
    6.6E+00,
    6.8E+00,
    7.0E+00,
    7.2E+00,
    7.4E+00,
    7.6E+00,
    7.8E+00,
    8.0E+00,
    8.2E+00,
    8.4E+00,
    8.6E+00,
    8.8E+00
};

// Variable $f$ from Table3b in http://arxiv.org/format/1006.3888v1
static const double Table3b_f[] = {
    0.000000000E+00,
    6.640999715E-03,
    2.655988402E-02,
    5.973463750E-02,
    1.061082208E-01,
    1.655717258E-01,
    2.379487173E-01,
    3.229815738E-01,
    4.203207655E-01,
    5.295180377E-01,
    6.500243699E-01,
    7.811933370E-01,
    9.222901256E-01,
    1.072505977E+00,
    1.230977302E+00,
    1.396808231E+00,
    1.569094960E+00,
    1.746950094E+00,
    1.929525170E+00,
    2.116029817E+00,
    2.305746418E+00,
    2.498039663E+00,
    2.692360938E+00,
    2.888247990E+00,
    3.085320655E+00,
    3.283273665E+00,
    3.481867612E+00,
    3.680919063E+00,
    3.880290678E+00,
    4.079881939E+00,
    4.279620923E+00,
    4.479457297E+00,
    4.679356615E+00,
    4.879295811E+00,
    5.079259772E+00,
    5.279238811E+00,
    5.479226847E+00,
    5.679220147E+00,
    5.879216466E+00,
    6.079214481E+00,
    6.279213431E+00,
    6.479212887E+00,
    6.679212609E+00,
    6.879212471E+00,
    7.079212403E+00
};

// Variable $f'$ from Table3b in http://arxiv.org/format/1006.3888v1
static const double Table3b_fp[] = {
    0.000000000E+00,
    6.640779210E-02,
    1.327641608E-01,
    1.989372524E-01,
    2.647091387E-01,
    3.297800312E-01,
    3.937761044E-01,
    4.562617647E-01,
    5.167567844E-01,
    5.747581439E-01,
    6.297657365E-01,
    6.813103772E-01,
    7.289819351E-01,
    7.724550211E-01,
    8.115096232E-01,
    8.460444437E-01,
    8.760814552E-01,
    9.017612214E-01,
    9.233296659E-01,
    9.411179967E-01,
    9.555182298E-01,
    9.669570738E-01,
    9.758708321E-01,
    9.826835008E-01,
    9.877895262E-01,
    9.915419002E-01,
    9.942455354E-01,
    9.961553040E-01,
    9.974777682E-01,
    9.983754937E-01,
    9.989728724E-01,
    9.993625417E-01,
    9.996117017E-01,
    9.997678702E-01,
    9.998638190E-01,
    9.999216041E-01,
    9.999557173E-01,
    9.999754577E-01,
    9.999866551E-01,
    9.999928812E-01,
    9.999962745E-01,
    9.999980875E-01,
    9.999990369E-01,
    9.999995242E-01,
    9.999997695E-01
};

// Variable $f''$ from Table3b in http://arxiv.org/format/1006.3888v1
static const double Table3b_fpp[] = {
    3.320573362E-01,
    3.319838371E-01,
    3.314698442E-01,
    3.300791276E-01,
    3.273892701E-01,
    3.230071167E-01,
    3.165891911E-01,
    3.078653918E-01,
    2.966634615E-01,
    2.829310173E-01,
    2.667515457E-01,
    2.483509132E-01,
    2.280917607E-01,
    2.064546268E-01,
    1.840065939E-01,
    1.613603195E-01,
    1.391280556E-01,
    1.178762461E-01,
    9.808627878E-02,
    8.012591814E-02,
    6.423412109E-02,
    5.051974749E-02,
    3.897261085E-02,
    2.948377201E-02,
    2.187118635E-02,
    1.590679869E-02,
    1.134178897E-02,
    7.927659815E-03,
    5.431957680E-03,
    3.648413667E-03,
    2.402039844E-03,
    1.550170691E-03,
    9.806151170E-04,
    6.080442648E-04,
    3.695625701E-04,
    2.201689553E-04,
    1.285698072E-04,
    7.359298339E-05,
    4.129031111E-05,
    2.270775140E-05,
    1.224092624E-05,
    6.467978611E-06,
    3.349939753E-06,
    1.700667989E-06,
    8.462841214E-07
};

gsl_spline * suzerain_blasius_spline_fit()
{
    static const size_t N = sizeof(Table3b_eta)/sizeof(Table3b_eta[0]);
    assert(N == sizeof(Table3b_f  )/sizeof(Table3b_f  [0]));
    assert(N == sizeof(Table3b_fp )/sizeof(Table3b_fp [0]));
    assert(N == sizeof(Table3b_fpp)/sizeof(Table3b_fpp[0]));
    gsl_spline * s = gsl_spline_alloc(gsl_interp_cspline, N);
    if (s) {
        if (gsl_spline_init(s, Table3b_eta, Table3b_fp, N)) {
            gsl_spline_free(s);
            s = NULL;
        }
    }
    return s;
}
