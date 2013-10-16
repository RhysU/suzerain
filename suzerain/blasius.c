/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013 Rhys Ulerich
 * Copyright (C) 2013 The PECOS Development Team
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
 * @copydoc blasius.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/blasius.h>

#include <gsl/gsl_nan.h>

#include <suzerain/common.h>

// C99 extern declaration for inlined function in blasius.h
extern
double suzerain_blasius_eta(double y, double Re_x);

// Variable $\eta$ from Table3b in http://arxiv.org/format/1006.3888v1
const double suzerain_blasius_ganapol_eta[] = {
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

enum {
    Nganapol = sizeof(suzerain_blasius_ganapol_eta)
             / sizeof(suzerain_blasius_ganapol_eta[0])
};

// Variable $f$ from Table3b in http://arxiv.org/format/1006.3888v1
const double suzerain_blasius_ganapol_f[Nganapol] = {
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
const double suzerain_blasius_ganapol_fp[Nganapol] = {
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
const double suzerain_blasius_ganapol_fpp[Nganapol] = {
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

// See "Compile Time Assertions" by Ralf Holly (http://drdobbs.com/184401873)
// Sanity check that the data arrays are all of equal size.
enum {
    assert1 = 1 / (    sizeof(suzerain_blasius_ganapol_f  )
                    == sizeof(suzerain_blasius_ganapol_eta)),
    assert2 = 1 / (    sizeof(suzerain_blasius_ganapol_fp )
                    == sizeof(suzerain_blasius_ganapol_eta)),
    assert3 = 1 / (    sizeof(suzerain_blasius_ganapol_fpp)
                    == sizeof(suzerain_blasius_ganapol_eta))
};

// See http://agentzlerich.blogspot.com/2013/10/generating-blasius-boundary-layer.html
// and https://github.com/RhysU/blasius for background on this data.
const double suzerain_blasius_extended_eta[] = {
    0.0000000000000000e+00,
    2.0000000000000001e-01,
    4.0000000000000002e-01,
    6.0000000000000009e-01,
    8.0000000000000004e-01,
    1.0000000000000000e+00,
    1.2000000000000002e+00,
    1.4000000000000001e+00,
    1.6000000000000001e+00,
    1.8000000000000000e+00,
    2.0000000000000000e+00,
    2.2000000000000002e+00,
    2.4000000000000004e+00,
    2.6000000000000001e+00,
    2.8000000000000003e+00,
    3.0000000000000000e+00,
    3.2000000000000002e+00,
    3.4000000000000004e+00,
    3.6000000000000001e+00,
    3.8000000000000003e+00,
    4.0000000000000000e+00,
    4.2000000000000002e+00,
    4.4000000000000004e+00,
    4.6000000000000005e+00,
    4.8000000000000007e+00,
    5.0000000000000000e+00,
    5.2000000000000002e+00,
    5.4000000000000004e+00,
    5.6000000000000005e+00,
    5.8000000000000007e+00,
    6.0000000000000000e+00,
    6.2000000000000002e+00,
    6.4000000000000004e+00,
    6.6000000000000005e+00,
    6.8000000000000007e+00,
    7.0000000000000000e+00,
    7.2000000000000002e+00,
    7.4000000000000004e+00,
    7.6000000000000005e+00,
    7.8000000000000007e+00,
    8.0000000000000000e+00,
    8.2000000000000011e+00,
    8.4000000000000004e+00,
    8.5999999999999996e+00,
    8.8000000000000007e+00,
    9.0000000000000000e+00,
    9.2000000000000011e+00,
    9.4000000000000004e+00,
    9.6000000000000014e+00,
    9.8000000000000007e+00,
    1.0000000000000000e+01,
    1.0200000000000001e+01,
    1.0400000000000000e+01,
    1.0600000000000001e+01,
    1.0800000000000001e+01,
    1.1000000000000000e+01,
    1.1200000000000001e+01,
    1.1400000000000000e+01,
    1.1600000000000001e+01,
    1.1800000000000001e+01,
    1.2000000000000000e+01,
    1.2200000000000001e+01,
    1.2400000000000000e+01,
    1.2600000000000001e+01,
    1.2800000000000001e+01,
    1.3000000000000000e+01,
    1.3200000000000001e+01,
    1.3400000000000000e+01,
    1.3600000000000001e+01,
    1.3800000000000001e+01,
    1.4000000000000000e+01,
    1.4200000000000001e+01,
    1.4400000000000000e+01,
    1.4600000000000001e+01,
    1.4800000000000001e+01,
    1.5000000000000000e+01,
    1.5200000000000001e+01,
    1.5400000000000000e+01,
    1.5600000000000001e+01,
    1.5800000000000001e+01,
    1.6000000000000000e+01,
    1.6199999999999999e+01,
    1.6400000000000002e+01,
    1.6600000000000001e+01,
    1.6800000000000001e+01,
    1.7000000000000000e+01,
    1.7199999999999999e+01,
    1.7400000000000002e+01,
    1.7600000000000001e+01,
    1.7800000000000001e+01,
    1.8000000000000000e+01,
    1.8199999999999999e+01,
    1.8400000000000002e+01,
    1.8600000000000001e+01,
    1.8800000000000001e+01,
    1.9000000000000000e+01,
    1.9200000000000003e+01,
    1.9400000000000002e+01,
    1.9600000000000001e+01,
    1.9800000000000001e+01,
    2.0000000000000000e+01,
    2.0200000000000003e+01,
    2.0400000000000002e+01,
    2.0600000000000001e+01,
    2.0800000000000001e+01,
    2.1000000000000000e+01,
    2.1200000000000003e+01,
    2.1400000000000002e+01,
    2.1600000000000001e+01,
    2.1800000000000001e+01,
    2.2000000000000000e+01,
    2.2200000000000003e+01,
    2.2400000000000002e+01,
    2.2600000000000001e+01,
    2.2800000000000001e+01,
    2.3000000000000000e+01,
    2.3200000000000003e+01,
    2.3400000000000002e+01,
    2.3600000000000001e+01,
    2.3800000000000001e+01,
    2.4000000000000000e+01,
    2.4200000000000003e+01,
    2.4400000000000002e+01,
    2.4600000000000001e+01,
    2.4800000000000001e+01,
    2.5000000000000000e+01,
    2.5200000000000003e+01,
    2.5400000000000002e+01,
    2.5600000000000001e+01,
    2.5800000000000001e+01,
    2.6000000000000000e+01,
    2.6200000000000003e+01,
    2.6400000000000002e+01,
    2.6600000000000001e+01,
    2.6800000000000001e+01,
    2.7000000000000000e+01,
    2.7200000000000003e+01,
    2.7400000000000002e+01,
    2.7600000000000001e+01,
    2.7800000000000001e+01,
    2.8000000000000000e+01,
    2.8200000000000003e+01,
    2.8400000000000002e+01,
    2.8600000000000001e+01,
    2.8800000000000001e+01,
    2.9000000000000000e+01,
    2.9200000000000003e+01,
    2.9400000000000002e+01,
    2.9600000000000001e+01,
    2.9800000000000001e+01,
    3.0000000000000000e+01,
    3.0200000000000003e+01,
    3.0400000000000002e+01,
    3.0600000000000001e+01,
    3.0800000000000001e+01,
    3.1000000000000000e+01,
    3.1200000000000003e+01,
    3.1400000000000002e+01,
    3.1600000000000001e+01,
    3.1800000000000001e+01,
    3.2000000000000000e+01,
    3.2200000000000003e+01,
    3.2399999999999999e+01,
    3.2600000000000001e+01,
    3.2800000000000004e+01,
    3.3000000000000000e+01,
    3.3200000000000003e+01,
    3.3399999999999999e+01,
    3.3600000000000001e+01,
    3.3800000000000004e+01,
    3.4000000000000000e+01,
    3.4200000000000003e+01,
    3.4399999999999999e+01,
    3.4600000000000001e+01,
    3.4800000000000004e+01,
    3.5000000000000000e+01,
    3.5200000000000003e+01,
    3.5399999999999999e+01,
    3.5600000000000001e+01,
    3.5800000000000004e+01,
    3.6000000000000000e+01,
    3.6200000000000003e+01,
    3.6399999999999999e+01,
    3.6600000000000001e+01,
    3.6800000000000004e+01,
    3.7000000000000000e+01,
    3.7200000000000003e+01,
    3.7399999999999999e+01,
    3.7600000000000001e+01,
    3.7800000000000004e+01,
    3.8000000000000000e+01,
    3.8200000000000003e+01,
    3.8400000000000006e+01,
    3.8600000000000001e+01
};

enum {
    Nextended = sizeof(suzerain_blasius_extended_eta)
              / sizeof(suzerain_blasius_extended_eta[0])
};

// See comments at suzerain_blasius_extended_eta
const double suzerain_blasius_extended_f[Nextended] = {
    0.0000000000000000e+00,
    6.6409997145970585e-03,
    2.6559884017994733e-02,
    5.9734637498038382e-02,
    1.0610822083904069e-01,
    1.6557172578927976e-01,
    2.3794871728892703e-01,
    3.2298157382953580e-01,
    4.2032076550162573e-01,
    5.2951803774381023e-01,
    6.5002436993528900e-01,
    7.8119333701057347e-01,
    9.2229012563454160e-01,
    1.0725059767878351e+00,
    1.2309773023143091e+00,
    1.3968082308703464e+00,
    1.5690949600067619e+00,
    1.7469500939488432e+00,
    1.9295251697605302e+00,
    2.1160298171720697e+00,
    2.3057464184620726e+00,
    2.4980396627069701e+00,
    2.6923609375430866e+00,
    2.8882479900178160e+00,
    3.0853206551774823e+00,
    3.2832736651563144e+00,
    3.4818676115086635e+00,
    3.6809190628975084e+00,
    3.8802906776333757e+00,
    4.0798819392396917e+00,
    4.2796209225138426e+00,
    4.4794572972826767e+00,
    4.6793566154314172e+00,
    4.8792958110602429e+00,
    5.0792597724490030e+00,
    5.2792388110290958e+00,
    5.4792268473431696e+00,
    5.6792201473338295e+00,
    5.8792164658040509e+00,
    6.0792144810719346e+00,
    6.2792134313460615e+00,
    6.4792128866784848e+00,
    6.6792126094414330e+00,
    6.8792124710150580e+00,
    7.0792124032167214e+00,
    7.2792123706452294e+00,
    7.4792123552968919e+00,
    7.6792123482031105e+00,
    7.8792123449874136e+00,
    8.0792123435577192e+00,
    8.2792123429343132e+00,
    8.4792123426677222e+00,
    8.6792123425559158e+00,
    8.8792123425099305e+00,
    9.0792123424913829e+00,
    9.2792123424840440e+00,
    9.4792123424811994e+00,
    9.6792123424801169e+00,
    9.8792123424797147e+00,
    1.0079212342479567e+01,
    1.0279212342479513e+01,
    1.0479212342479494e+01,
    1.0679212342479486e+01,
    1.0879212342479486e+01,
    1.1079212342479485e+01,
    1.1279212342479484e+01,
    1.1479212342479485e+01,
    1.1679212342479484e+01,
    1.1879212342479486e+01,
    1.2079212342479485e+01,
    1.2279212342479484e+01,
    1.2479212342479485e+01,
    1.2679212342479484e+01,
    1.2879212342479486e+01,
    1.3079212342479485e+01,
    1.3279212342479484e+01,
    1.3479212342479485e+01,
    1.3679212342479484e+01,
    1.3879212342479486e+01,
    1.4079212342479485e+01,
    1.4279212342479484e+01,
    1.4479212342479483e+01,
    1.4679212342479486e+01,
    1.4879212342479486e+01,
    1.5079212342479485e+01,
    1.5279212342479484e+01,
    1.5479212342479483e+01,
    1.5679212342479486e+01,
    1.5879212342479486e+01,
    1.6079212342479483e+01,
    1.6279212342479482e+01,
    1.6479212342479482e+01,
    1.6679212342479484e+01,
    1.6879212342479484e+01,
    1.7079212342479483e+01,
    1.7279212342479482e+01,
    1.7479212342479485e+01,
    1.7679212342479484e+01,
    1.7879212342479484e+01,
    1.8079212342479483e+01,
    1.8279212342479482e+01,
    1.8479212342479485e+01,
    1.8679212342479484e+01,
    1.8879212342479484e+01,
    1.9079212342479483e+01,
    1.9279212342479482e+01,
    1.9479212342479485e+01,
    1.9679212342479484e+01,
    1.9879212342479484e+01,
    2.0079212342479483e+01,
    2.0279212342479482e+01,
    2.0479212342479485e+01,
    2.0679212342479484e+01,
    2.0879212342479484e+01,
    2.1079212342479483e+01,
    2.1279212342479482e+01,
    2.1479212342479485e+01,
    2.1679212342479484e+01,
    2.1879212342479484e+01,
    2.2079212342479483e+01,
    2.2279212342479482e+01,
    2.2479212342479485e+01,
    2.2679212342479484e+01,
    2.2879212342479484e+01,
    2.3079212342479483e+01,
    2.3279212342479482e+01,
    2.3479212342479485e+01,
    2.3679212342479484e+01,
    2.3879212342479484e+01,
    2.4079212342479483e+01,
    2.4279212342479482e+01,
    2.4479212342479485e+01,
    2.4679212342479484e+01,
    2.4879212342479484e+01,
    2.5079212342479483e+01,
    2.5279212342479482e+01,
    2.5479212342479485e+01,
    2.5679212342479484e+01,
    2.5879212342479484e+01,
    2.6079212342479483e+01,
    2.6279212342479482e+01,
    2.6479212342479485e+01,
    2.6679212342479484e+01,
    2.6879212342479484e+01,
    2.7079212342479483e+01,
    2.7279212342479482e+01,
    2.7479212342479485e+01,
    2.7679212342479484e+01,
    2.7879212342479484e+01,
    2.8079212342479483e+01,
    2.8279212342479482e+01,
    2.8479212342479485e+01,
    2.8679212342479484e+01,
    2.8879212342479484e+01,
    2.9079212342479483e+01,
    2.9279212342479482e+01,
    2.9479212342479485e+01,
    2.9679212342479484e+01,
    2.9879212342479484e+01,
    3.0079212342479483e+01,
    3.0279212342479482e+01,
    3.0479212342479485e+01,
    3.0679212342479481e+01,
    3.0879212342479484e+01,
    3.1079212342479487e+01,
    3.1279212342479482e+01,
    3.1479212342479485e+01,
    3.1679212342479481e+01,
    3.1879212342479484e+01,
    3.2079212342479487e+01,
    3.2279212342479482e+01,
    3.2479212342479485e+01,
    3.2679212342479481e+01,
    3.2879212342479484e+01,
    3.3079212342479487e+01,
    3.3279212342479482e+01,
    3.3479212342479485e+01,
    3.3679212342479481e+01,
    3.3879212342479484e+01,
    3.4079212342479487e+01,
    3.4279212342479482e+01,
    3.4479212342479485e+01,
    3.4679212342479481e+01,
    3.4879212342479484e+01,
    3.5079212342479487e+01,
    3.5279212342479482e+01,
    3.5479212342479485e+01,
    3.5679212342479481e+01,
    3.5879212342479484e+01,
    3.6079212342479487e+01,
    3.6279212342479482e+01,
    3.6479212342479485e+01,
    3.6679212342479488e+01,
    3.6879212342479484e+01
};

const double suzerain_blasius_extended_fp[Nextended] = {
    0.0000000000000000e+00,
    6.6407792096250848e-02,
    1.3276416076102221e-01,
    1.9893725242221899e-01,
    2.6470913872311741e-01,
    3.2978003124966698e-01,
    3.9377610443395650e-01,
    4.5626176470513380e-01,
    5.1675678442261486e-01,
    5.7475814388945690e-01,
    6.2976573650238588e-01,
    6.8131037723602750e-01,
    7.2898193506257503e-01,
    7.7245502114856535e-01,
    8.1150962319910902e-01,
    8.4604444365799358e-01,
    8.7608145517247749e-01,
    9.0176122143676474e-01,
    9.2332966588164878e-01,
    9.4111799672931340e-01,
    9.5551822981069434e-01,
    9.6695707375360584e-01,
    9.7587083213647785e-01,
    9.8268350076070032e-01,
    9.8778952620078742e-01,
    9.9154190016439347e-01,
    9.9424553535992810e-01,
    9.9615530396275309e-01,
    9.9747776821291501e-01,
    9.9837549365019085e-01,
    9.9897287243586053e-01,
    9.9936254171909067e-01,
    9.9961170171767588e-01,
    9.9976787021003244e-01,
    9.9986381903703803e-01,
    9.9992160414794995e-01,
    9.9995571727792665e-01,
    9.9997545768489782e-01,
    9.9998665513912843e-01,
    9.9999288116610119e-01,
    9.9999627453530060e-01,
    9.9999808745923757e-01,
    9.9999903687484326e-01,
    9.9999952424780447e-01,
    9.9999976948975089e-01,
    9.9999989045379900e-01,
    9.9999994893901312e-01,
    9.9999997665715090e-01,
    9.9999998953401714e-01,
    9.9999999539788897e-01,
    9.9999999801539019e-01,
    9.9999999916068538e-01,
    9.9999999965190545e-01,
    9.9999999985842569e-01,
    9.9999999994353506e-01,
    9.9999999997791622e-01,
    9.9999999999153044e-01,
    9.9999999999681477e-01,
    9.9999999999882527e-01,
    9.9999999999957512e-01,
    9.9999999999984923e-01,
    9.9999999999994749e-01,
    9.9999999999998201e-01,
    9.9999999999999389e-01,
    9.9999999999999789e-01,
    9.9999999999999922e-01,
    9.9999999999999967e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01,
    9.9999999999999978e-01
};

const double suzerain_blasius_extended_fpp[Nextended] = {
    3.3205733621519629e-01,
    3.3198383711462798e-01,
    3.3146984420145503e-01,
    3.3007912757429625e-01,
    3.2738927014925229e-01,
    3.2300711668694282e-01,
    3.1658919106111544e-01,
    3.0786539179016786e-01,
    2.9666346145571931e-01,
    2.8293101725975589e-01,
    2.6675154569727827e-01,
    2.4835091319037131e-01,
    2.2809176068668360e-01,
    2.0645462679942550e-01,
    1.8400659386536070e-01,
    1.6136031954087834e-01,
    1.3912805557242577e-01,
    1.1787624608752340e-01,
    9.8086278784280084e-02,
    8.0125918139197061e-02,
    6.4234121091690577e-02,
    5.0519747486645686e-02,
    3.8972610853629498e-02,
    2.9483772011648642e-02,
    2.1871186347443238e-02,
    1.5906798685318153e-02,
    1.1341788968929135e-02,
    7.9276598147061551e-03,
    5.4319576799275147e-03,
    3.6484136667473462e-03,
    2.4020398437572996e-03,
    1.5501706906563828e-03,
    9.8061511700977232e-04,
    6.0804426478449074e-04,
    3.6956257014031428e-04,
    2.2016895527113460e-04,
    1.2856980723515674e-04,
    7.3592983389223794e-05,
    4.1290311113698544e-05,
    2.2707751402803766e-05,
    1.2240926243249920e-05,
    6.4679786108453658e-06,
    3.3499397531974408e-06,
    1.7006679885717248e-06,
    8.4628412130585306e-07,
    4.1278790156475232e-07,
    1.9735668380412326e-07,
    9.2489158677487331e-08,
    4.2485812620375726e-08,
    1.9129831304430804e-08,
    8.4429158651022967e-09,
    3.6524803913092198e-09,
    1.5488074704694815e-09,
    6.4375569794427733e-10,
    2.6227617757128940e-10,
    1.0473955128094568e-10,
    4.0999322474265197e-11,
    1.5731015212516414e-11,
    5.9163099772618016e-12,
    2.1810177063481187e-12,
    7.8810042675859475e-13,
    2.7913740190928823e-13,
    9.6910000670083129e-14,
    3.2978678842745192e-14,
    1.1000489193601443e-14,
    3.5967050786090774e-15,
    1.1526879056075572e-15,
    3.6210349563133137e-16,
    1.1149817624085032e-16,
    3.3652463895846418e-17,
    9.9558883380119627e-18,
    2.8870692572896716e-18,
    8.2063192574246674e-19,
    2.2864074180800616e-19,
    6.2441426902306141e-20,
    1.6714984659791484e-20,
    4.3858430443215050e-21,
    1.1280129861182341e-21,
    2.8437341294364077e-22,
    7.0271256502345964e-23,
    1.7020810680574499e-23,
    4.0410712223087798e-24,
    9.4042994774277148e-25,
    2.1452109779172718e-25,
    4.7965293589945852e-26,
    1.0512297679161997e-26,
    2.2582995021569423e-27,
    4.7553080499735291e-28,
    9.8149675730042055e-29,
    1.9856926465778935e-29,
    3.9377483627391731e-30,
    7.6541412156787733e-31,
    1.4583351355181007e-31,
    2.7235184131024432e-32,
    4.9855734246948443e-33,
    8.9456401044497370e-34,
    1.5733265449215635e-34,
    2.7122958127887408e-35,
    4.5831652130786345e-36,
    7.5910828469873827e-37,
    1.2323984036641570e-37,
    1.9611329988618580e-38,
    3.0589385941978318e-39,
    4.6767207619714147e-40,
    7.0083885099190416e-41,
    1.0294373363811153e-41,
    1.4821273674923266e-42,
    2.0915763086607784e-43,
    2.8930972939667386e-44,
    3.9224001552355360e-45,
    5.2124091708031195e-46,
    6.7892352135297409e-47,
    8.6675577626222035e-48,
    1.0845852714395082e-48,
    1.3302062620067328e-49,
    1.5990408451522715e-50,
    1.8840049200109247e-51,
    2.1756181226916025e-52,
    2.4623926789265991e-53,
    2.7315010074974225e-54,
    2.9696805873471981e-55,
    3.1642938498112364e-56,
    3.3044283832872253e-57,
    3.3819055610426380e-58,
    3.3920684938547889e-59,
    3.3342440151369026e-60,
    3.2118151386126913e-61,
    3.0318933015271105e-62,
    2.8046345683889967e-63,
    2.5422912249132635e-64,
    2.2581218340582574e-65,
    1.9652940874035787e-66,
    1.6759049914493029e-67,
    1.4002153927008813e-68,
    1.1461568792521680e-69,
    9.1912644543839046e-71,
    7.2204543462141436e-72,
    5.5562992217299953e-73,
    4.1880300588302606e-74,
    3.0917579003585325e-75,
    2.2353118803015678e-76,
    1.5825952580418952e-77,
    1.0971328044714600e-78,
    7.4466463321320530e-80,
    4.9479433619150357e-81,
    3.2180793936666466e-82,
    2.0484118520586418e-83,
    1.2759132252044612e-84,
    7.7756562082842033e-86,
    4.6353981261175121e-87,
    2.7026251888740908e-88,
    1.5407780980219590e-89,
    8.5871257260844221e-91,
    4.6773229889612407e-92,
    2.4892388818274261e-93,
    1.2939628648472887e-94,
    6.5677525065747211e-96,
    3.2538033326029564e-97,
    1.5727870853228035e-98,
    7.4141340404797473e-100,
    3.4068421671381657e-101,
    1.5251530497264114e-102,
    6.6479965692561246e-104,
    2.8197143445244591e-105,
    1.1629192637092119e-106,
    4.6600038355110267e-108,
    1.8127702002052962e-109,
    6.8392386592324925e-111,
    2.4999224286723917e-112,
    8.8429116394024642e-114,
    3.0231156782214689e-115,
    9.9743067437888638e-117,
    3.1709178257332913e-118,
    9.6958279486040461e-120,
    2.8458406138630829e-121,
    7.9998070326770773e-123,
    2.1482251580718956e-124,
    5.4947361028691160e-126,
    1.3342409449725092e-127,
    3.0639171774547037e-129,
    6.6243587249864036e-131,
    1.3414513702652775e-132,
    2.5287078430355464e-134,
    4.4047117308988946e-136,
    7.0266393543409267e-138,
    1.0152807419752380e-139,
    1.3102634354511120e-141,
    1.4830676323939600e-143,
    1.4365356668016854e-145,
    1.1498674445293797e-147,
    7.2098563079040367e-150,
    3.2299875759095029e-152,
    8.4957414420568161e-155,
    6.0946204192285745e-158
};

// See "Compile Time Assertions" by Ralf Holly (http://drdobbs.com/184401873)
// Sanity check that the data arrays are all of equal size.
enum {
    assert4 = 1 / (    sizeof(suzerain_blasius_extended_f  )
                    == sizeof(suzerain_blasius_extended_eta)),
    assert5 = 1 / (    sizeof(suzerain_blasius_extended_fp )
                    == sizeof(suzerain_blasius_extended_eta)),
    assert6 = 1 / (    sizeof(suzerain_blasius_extended_fpp)
                    == sizeof(suzerain_blasius_extended_eta))
};

// Handles the boilerplate aspects of preparing a spline fit
// Also isolates spline-type selections in one place
static
gsl_spline * prepare_fit(
        const size_t N,
        const double * const x,
        const double * const y)
{
    gsl_spline * s = gsl_spline_alloc(gsl_interp_akima, N);
    if (s) {
        if (gsl_spline_init(s, x, y, N)) {
            gsl_spline_free(s);
            s = NULL;
        }
    }
    return s;
}

// TODO Use non-natural AKIMA conditions to set known derivatives
// That is, suzerain_blasius_extended_fpp and f'''(eta) are
// known and provide important information wrt skin friction.
gsl_spline * suzerain_blasius_u(const double Re_x)
{
    double y[Nextended];
    const double invSqrtRe = sqrt(1 / Re_x);
    for (size_t i = 0; i < Nextended; ++i) {
        y[i] = invSqrtRe * suzerain_blasius_extended_eta[i];
    }
    return prepare_fit(Nextended, y, suzerain_blasius_extended_fp);
}

gsl_spline * suzerain_blasius_v(const double Re_x)
{
    double y[Nextended], v[Nextended];
    const double invSqrtRe = sqrt(1 / Re_x);
    for (size_t i = 0; i < Nextended; ++i) {
        y[i] = invSqrtRe * suzerain_blasius_extended_eta[i];
        v[i] = invSqrtRe/M_SQRT2 * (  suzerain_blasius_extended_f  [i]
                                    + suzerain_blasius_extended_eta[i]
                                    * suzerain_blasius_extended_fp [i]);
    }
    return prepare_fit(Nextended, y, v);
}

gsl_spline * suzerain_blasius_T(const double Re_x, const double Pr)
{
    return suzerain_blasius_u(Re_x * Pr);
}
