//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc cev.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/cev.hpp>

#include <gsl/gsl_spline.h>

#include <suzerain/common.hpp>

namespace suzerain {

namespace cev {

// See method Doxygen for origin of these magic numbers.
int iss_laminar(
    const double dstag,
    double& gammae,
    double& Mae,
    double& pexi,
    double& T_ratio)
{
    enum { N = 71 };
    static const double data_dstag[N] = {
        -0.6143492210355781,    -0.59390333664534101,  -0.57179300778174724,
        -0.54789092718131283,   -0.52206122007465172,  -0.49415938719521613,
        -0.46403232768676173,   -0.43151718546930518,  -0.3964421451945217,
        -0.35862567789533895,   -0.31787735472377765,  -0.27399655136109757,
        -0.2267747407134193,    -0.17599395375299576,  -0.12142773350753355,
        -0.062843090903680121,  -0,                    0.067346704177323247,
        0.13944580580883326,    0.21654744952740934,   0.29890258933871827,
        0.38675857037813088,    0.48035747367467119,   0.57993312475275527,
        0.68103062660465374,    0.78212886925651937,   0.88322732787220115,
        0.9843255166409417,     1.0854233790807597,    1.186520826236158,
        1.2876186095690088,     1.3887167424250881,    1.4898149297790302,
        1.5909130106137055,     1.6920108449633222,    1.7931083547689379,
        1.8942058645744817,     1.9953036989240984,    2.0964017797587737,
        2.1974999671127158,     2.2985980999687952,    2.399695883301646,
        2.5007933304570438,     2.6018911928968622,    2.7029893816656028,
        2.8040878402812845,     2.9051860829331497,    3.0062835847850486,
        3.1058592358631327,     3.199458139159673,     3.2873141201990856,
        3.3696692600103946,     3.4467709037289707,    3.5188700053604807,
        3.5862167095378044,     3.6490598004414849,    3.7076444430453384,
        3.7622106632908006,     3.8129914502512245,    3.8602132608989028,
        3.9040940642615829,     3.9448423874331442,    3.9826588547323265,
        4.01773389500711,       4.0502490372245665,    4.0803760967330209,
        4.1082779296124556,     4.1341076367191167,    4.1580097173195512,
        4.1801200461831449,     4.2005659305733829
    };
    static const double data_gammae[N] = {
        1.4099815429390452,  1.4098104536887972,  1.4096908057455511,
        1.4096108296401819,  1.4095556146292538,  1.4095189225317255,
        1.4094954447208019,  1.4094812278433841,  1.4094741105376656,
        1.4094728832908761,  1.4094750964867506,  1.409479332470402,
        1.4094845509118477,  1.4094896108287147,  1.409493739436265,
        1.409496087805638,   1.4094965858307458,  1.4094949070017295,
        1.4094907750622041,  1.4094836020905861,  1.4094729365300411,
        1.4094581868477518,  1.409438670326961,   1.4094139594507757,
        1.4093838947568316,  1.409349694990228,   1.4093108365063949,
        1.4092681785505701,  1.4092214742059146,  1.4091712712490843,
        1.4091178420751405,  1.4090617036849444,  1.4090029840768459,
        1.4089424630791352,  1.4088808377294937,  1.4088184299707727,
        1.4087553064055358,  1.4086914743946093,  1.4086264139081992,
        1.4085601149390066,  1.4084929150788372,  1.4084256297807161,
        1.4083595370046604,  1.4082956118387255,  1.4082345113997261,
        1.4081769445147463,  1.4081242617376561,  1.4080790926320617,
        1.4080536636759231,  1.4080294579155366,  1.4080078090065777,
        1.4079951744388131,  1.407991112426707,   1.4079948516205061,
        1.4080058983225989,  1.4080234319050933,  1.4080461078617779,
        1.4080733203251381,  1.4081050590057502,  1.4081403498844409,
        1.4081787684846814,  1.4082204150965079,  1.408265818053513,
        1.4083167289094889,  1.4083781550996435,  1.4084445755187867,
        1.4085010844389505,  1.4085542572945138,  1.4086045117223955,
        1.4086554301854781,  1.4087083877256554
    };
    static const double data_Mae[N] = {
        0.42190694027972758,   0.38261248638894363,    0.3477875752073547,
        0.31640737083790416,   0.2876995764882363,     0.25985405007086548,
        0.23489666492710973,   0.21197219841991385,    0.18901267536866218,
        0.16560121027857475,   0.14356407281720557,    0.12173776923151167,
        0.099465572885036532,  0.076881172916196644,   0.054136052992762411,
        0.030073499145841599,  0.0054429205954211607,  0.020625533435922171,
        0.04710689731459318,   0.074540336218711409,   0.10229648360013355,
        0.13005033212091255,   0.15869200176946574,    0.18735787088001371,
        0.21665239131864536,   0.24514188073129567,    0.27350049613561794,
        0.30155122032277992,   0.32914547794748539,    0.35667346504868797,
        0.38391132699969333,   0.41120414466739352,    0.43862203023951135,
        0.46611985089770513,   0.49373079729648561,    0.52144481034749257,
        0.54927217716869747,   0.57704944407135428,    0.60480443623717972,
        0.63230228827088442,   0.65970364956076455,    0.68695701225763917,
        0.71394455827465408,   0.7409908147908203,     0.76817051568966199,
        0.79517940351262706,   0.82219640239697489,    0.84900156497529833,
        0.87490773103516306,   0.89858115487825818,    0.92143086386963002,
        0.94344288692974954,   0.9637540241995538,     0.98245494114359588,
        1.000847126522,        1.0177358093237081,     1.0332195780131908,
        1.0481561893730225,    1.0623982037546464,     1.075287642592152,
        1.0874605127903292,    1.0985922749827015,     1.1093731426053473,
        1.118775131627562,     1.1272786975850249,     1.1348406174944963,
        1.1433080447195736,    1.1521944684829317,     1.1621859200981024,
        1.1744408936560817,    1.1905934935259677
    };
    static const double data_pexi[N] = {
        -0.10101086946499484,    -0.095023446053762964,   -0.089182343095845121,
        -0.083447045126352787,   -0.079212639987524411,   -0.077883638614280201,
        -0.076400383502228372,   -0.075017607543867923,   -0.076596710412302649,
        -0.081926940030949899,   -0.087216018089626909,   -0.095805208115930338,
        -0.11274893621616011,    -0.14376217440601388,    -0.2019320381281717,
        -0.39128487154052993,    -3.5106047841435903,     -0.32211414005226363,
        -0.17260017496292279,    -0.10940150511034639,    -0.077877994719461052,
        -0.059829843422088434,   -0.047537526508757005,   -0.038927185529409138,
        -0.033207371957848009,   -0.029286208180788198,   -0.026241315356759633,
        -0.023850353644973017,   -0.021862794866939709,   -0.0202783088695449,
        -0.018972065957037879,   -0.017931635525217017,   -0.017105486287050446,
        -0.016420139593348454,   -0.01582059702106817,    -0.01527493730587331,
        -0.014755313619465062,   -0.014217427339774808,   -0.013690304564595569,
        -0.013159861199299431,   -0.012694892400477699,   -0.012285393056196152,
        -0.011892537649319856,   -0.011571158658678039,   -0.011314885044750676,
        -0.01104037110704819,    -0.010840033105601166,   -0.010830296722694014,
        -0.010573127574118263,   -0.010191783208898577,   -0.010111390042733736,
        -0.010087225149769299,   -0.010004726557591544,   -0.0098711124687766924,
        -0.0098858336057264413,  -0.0098399105377029502,  -0.009736713144315564,
        -0.0097314199000317486,  -0.009786887911991721,   -0.0097198961918436976,
        -0.0096911626691443126,  -0.0096502759056394972,  -0.0097415965750651881,
        -0.0097736338699843522,  -0.010045362685481096,   -0.010248414045932849,
        -0.011066019395791476,   -0.012350269196758374,   -0.014683079037586257,
        -0.018803670726498776,   -0.025438605531063266
    };
    static const double data_T_ratio[N] = {
        3.4872999452040814,  3.6479923945555432,  3.7199291691044016,
        3.7644092644152076,  3.7999017268666924,  3.8313715997967619,
        3.8602181515453204,  3.8857189776984429,  3.9097235134249479,
        3.9313475230762478,  3.9505823320748172,  3.9669976109161964,
        3.9807965955454621,  3.9918326739643462,  4.0007546597793091,
        4.0081656310024787,  4.015032259169498,   4.0221316905250157,
        4.0308078701497205,  4.0407977228460457,  4.0515352658071633,
        4.0623404117234481,  4.0725715903486446,  4.0830141612741606,
        4.0933729464407564,  4.1010991108409218,  4.1072300067106804,
        4.1125933915687005,  4.1174438736460042,  4.1217382515014434,
        4.1254436508821914,  4.1291175669260829,  4.1334874468323957,
        4.1379994390070687,  4.142982497780304,   4.1482700949021858,
        4.1541354238844006,  4.1611072377993317,  4.1677384625970166,
        4.1748141000911261,  4.181692718061174,   4.1885944038008942,
        4.1960728537236802,  4.2035455920447395,  4.2112397768750833,
        4.2197680032275438,  4.2295387519018597,  4.2400864198458468,
        4.2510425457759693,  4.2619128750651454,  4.2714687975204821,
        4.2805456956631689,  4.2880414432341372,  4.2952355514753942,
        4.3010590003939635,  4.3057514460810973,  4.3102253879254713,
        4.3134371768417648,  4.3161547307280461,  4.3181456381433216,
        4.3196888662022648,  4.3203828218121583,  4.3205184682086744,
        4.3194370500077701,  4.3168567736825629,  4.3119612190476326,
        4.302792717238483,   4.2853609220662268,  4.2545325614550382,
        4.1864251808993052,  4.0039925486698449
    };

    if (dstag < data_dstag[0]) {
        std::ostringstream oss;
        oss << "Requested dstag " << dstag
            << " is less than minimum " << data_dstag[0];
        throw std::domain_error(oss.str());
    }

    if (dstag > data_dstag[N-1]) {
        std::ostringstream oss;
        oss << "Requested dstag " << dstag
            << " is greater than maximum " << data_dstag[N-1];
        throw std::domain_error(oss.str());
    }

    shared_ptr<gsl_interp_accel> a(gsl_interp_accel_alloc(),
                                   gsl_interp_accel_free);
    shared_ptr<gsl_interp>       f(gsl_interp_alloc(gsl_interp_cspline, N),
                                   gsl_interp_free);

    gsl_interp_init(f.get(), data_dstag, data_gammae, N);
    gammae = gsl_interp_eval(f.get(), data_dstag, data_gammae, dstag, a.get());

    gsl_interp_init(f.get(), data_dstag, data_Mae, N);
    Mae = gsl_interp_eval(f.get(), data_dstag, data_Mae, dstag, a.get());

    gsl_interp_init(f.get(), data_dstag, data_pexi, N);
    pexi = gsl_interp_eval(f.get(), data_dstag, data_pexi, dstag, a.get());

    gsl_interp_init(f.get(), data_dstag, data_T_ratio, N);
    T_ratio = gsl_interp_eval(f.get(), data_dstag, data_T_ratio, dstag, a.get());

    return 0;
}

} // namespace cev

} // namespace suzerain
