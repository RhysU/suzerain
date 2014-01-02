//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#include <suzerain/rholut.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

#pragma warning(disable:1599)

using suzerain::Vector3r;
using suzerain::Matrix3r;

// Provides data for the test field
//      rho = 2*(x^2)*y*z + 3*x*(y^2)*z + 5*x*y*(z^2)
//      m   = {
//               7*sin( 7*x) + 11*sin(11*y) + 13*sin(13*z),
//              17*sin(17*x) + 19*sin(19*y) + 23*sin(23*z),
//              31*sin(31*x) + 37*sin(37*y) + 41*sin(41*z)
//      }
//      e   = 311*(x^2)*y*z + 313*x*(y^2)*z + 317*x*y*(z^2)
// symbolically evaluated at (x,y,z) = (1, 2, 3)
//
// Floating point values evaluated using test_rholut.sage
// to 200 bits of precision.
static
void rholut_test_data(
        double   &rho,
        Vector3r &grad_rho,
        double   &div_grad_rho,
        Matrix3r &grad_grad_rho,
        Vector3r &m,
        double   &div_m,
        Matrix3r &grad_m,
        Vector3r &div_grad_m,
        Vector3r &grad_div_m,
        double   &e,
        Vector3r &grad_e,
        double   &div_grad_e)
{
    rho = 138.;

    grad_rho(0) = 150.;
    grad_rho(1) =  87.;
    grad_rho(2) =  76.;

    div_grad_rho = 62.;

    grad_grad_rho(0,0) = 24.;
    grad_grad_rho(0,1) = 93.;
    grad_grad_rho(0,2) = 80.;
    grad_grad_rho(1,0) = 93.;
    grad_grad_rho(1,1) = 18.;
    grad_grad_rho(1,2) = 44.;
    grad_grad_rho(2,0) = 80.;
    grad_grad_rho(2,1) = 44.;
    grad_grad_rho(2,2) = 20.;

    m(0) =  17.030881810530221790029045449138418084725998922603314999607L;
    m(1) = -13.352805083487451615530843887500866367157365945228654781328L;
    m(2) = -67.831621760613409212021683467570521921918364223433827929255L;

    div_m = -1110.9529361851136052549019796872846553765997751760446179283L;

    grad_m(0,0) =    36.941210462821927268918678564239918598697356181546215137637L;
    grad_m(0,1) =  -120.99525999375109230095514344473594667045460456500973446019L;
    grad_m(0,2) =    45.062655568829395508033559350881935096816196590558655110591L;
    grad_m(1,0) =  - 79.522204696911510521679027743488377694525395918294315600644L;
    grad_m(1,1) =   344.78158550107344358873958098152158053944887665250279502393L;
    grad_m(1,2) =   525.50351087308169627005787890310679646720537555710570892294L;
    grad_m(2,0) =   879.06740585015455908290082893049093028384492918816862062545L;
    grad_m(2,1) =   235.08104096633447429877919230313280932462835892676775285483L;
    grad_m(2,2) = -1492.6757321490089761125602392330461545147460080100936280899L;

    div_grad_m(0) = - 2331.0237743611578930680818995644668983923397087883721673025L;
    div_grad_m(1) =   4087.1406255366278601045485877398647158902558157286863079530L;
    div_grad_m(2) =  93634.307505134882301285667998181793756271588535303593390142L;

    grad_div_m(0) = -  225.34640336054465800617068841661686603234641360063909560764L;
    grad_div_m(1) = - 2032.7920813676738919937627213144098066291866038214905756451L;
    grad_div_m(2) =  31697.008481817318630325961933724829632738184041308835720121L;

    e = 11328.;

    grad_e(0) = 13194.;
    grad_e(1) =  7542.;
    grad_e(2) =  5678.;

    div_grad_e = 6878.;
}

// Helper macro for loading test data via rhoult_test_data
#define ADD_TEST_DATA_TO_ENCLOSING_SCOPE             \
    double          rho;                             \
    Vector3r grad_rho;                               \
    double          div_grad_rho;                    \
    Matrix3r grad_grad_rho;                          \
    Vector3r m;                                      \
    double          div_m;                           \
    Matrix3r grad_m;                                 \
    Vector3r div_grad_m;                             \
    Vector3r grad_div_m;                             \
    double          e;                               \
    Vector3r grad_e;                                 \
    double          div_grad_e;                      \
                                                     \
    rholut_test_data(                                \
        rho, grad_rho, div_grad_rho, grad_grad_rho,  \
        m, div_m, grad_m, div_grad_m, grad_div_m,    \
        e, grad_e, div_grad_e)

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const Vector3r u = suzerain::rholut::u(rho, m);

    /* Expected results found using test_rholut.sage */
    const Vector3r ans(
             0.12341218703282769413064525687781462380236231103335735506961L,
            -0.096759457126720663880658289039861350486647579313251121603828L,
            -0.49153349101893774791320060483746755015882872625676686905257L);

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e2;
    BOOST_CHECK_CLOSE(u(0), ans(0), close_enough);
    BOOST_CHECK_CLOSE(u(1), ans(1), close_enough);
    BOOST_CHECK_CLOSE(u(2), ans(2), close_enough);
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_p_T_mu_lambda )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    /* Variant that provides derivatives */
    suzerain::rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* Expected results found using test_rholut.sage */
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;
    BOOST_CHECK_CLOSE(p,
            4441.1984111498830681828969449447964342868418205971441637701L,
            close_enough);
    BOOST_CHECK_CLOSE(T,
            45.055636055143743363613040114696980466792183873003314822257L,
            close_enough);
    BOOST_CHECK_CLOSE(mu,
            12.661915668902014668188456691480749088816699629988322914800L,
            close_enough);
    BOOST_CHECK_CLOSE(lambda,
            54.868301231908730228816645663083246051539031729949399297467L,
            close_enough);

    const Vector3r grad_p_ans(
            7432.6318085149443001747520526182772856422787520432873299439L,
            3876.3726619344082666965507166633629625480052934837862921105L,
            -1052.4623117340587311081836624513934718452258172103981784331L);
    BOOST_CHECK_CLOSE(grad_p(0), grad_p_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_p(1), grad_p_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_p(2), grad_p_ans(2), close_enough);

    const Vector3r grad_T_ans(
            26.429993649633053564012507465458338360389321254158523061160L,
            10.920879637033810058447318598785230160203637285265515837835L,
            -35.490402729120338418589989472089860721815316114333655060589L);
    BOOST_CHECK_CLOSE(grad_T(0), grad_T_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_T(1), grad_T_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_T(2), grad_T_ans(2), close_enough);

    const Vector3r grad_mu_ans(
            4.9517201401349647231331890756487079296179499785864492804454L,
            2.0460519349175446824606520620070658148800973612906173152405L,
            -6.6492086341354822416604185512988943122060233671691387576565L);
    BOOST_CHECK_CLOSE(grad_mu(0), grad_mu_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_mu(1), grad_mu_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_mu(2), grad_mu_ans(2), close_enough);

    const Vector3r grad_lambda_ans(
            21.457453940584847133577152661144401028344449907207946881930L,
            8.8662250513093602906628256020306185311470885655926750327087L,
            -28.813237414587089713861813722295208686226101257732934616511L);
    BOOST_CHECK_CLOSE(grad_lambda(0), grad_lambda_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_lambda(1), grad_lambda_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_lambda(2), grad_lambda_ans(2), close_enough);

    /* Check computation of grad_e from gradient of internal, kinetic */
    const Vector3r alt_grad_e
        = suzerain::rholut::energy_internal_gradient(gamma, grad_p)
        + suzerain::rholut::energy_kinetic_gradient(Ma,rho,grad_rho,m,grad_m);
    BOOST_CHECK_CLOSE(grad_e(0), alt_grad_e(0), close_enough);
    BOOST_CHECK_CLOSE(grad_e(1), alt_grad_e(1), close_enough);
    BOOST_CHECK_CLOSE(grad_e(2), alt_grad_e(2), close_enough);

    /* Variant that does not provide derivatives */
    p = T = mu = lambda = 555;
    suzerain::rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, m, e, p, T, mu, lambda);

    BOOST_CHECK_CLOSE(p,
            4441.1984111498830681828969449447964342868418205971441637701L,
            close_enough);
    BOOST_CHECK_CLOSE(T,
            45.055636055143743363613040114696980466792183873003314822257L,
            close_enough);
    BOOST_CHECK_CLOSE(mu,
            12.661915668902014668188456691480749088816699629988322914800L,
            close_enough);
    BOOST_CHECK_CLOSE(lambda,
            54.868301231908730228816645663083246051539031729949399297467L,
            close_enough);

    /* Variant that does not provide viscosities */
    p = T = 555;
    suzerain::rholut::p_T(alpha, beta, gamma, Ma, rho, m, e, p, T);

    BOOST_CHECK_CLOSE(p,
            4441.1984111498830681828969449447964342868418205971441637701L,
            close_enough);
    BOOST_CHECK_CLOSE(T,
            45.055636055143743363613040114696980466792183873003314822257L,
            close_enough);

    /* Variant that does not provide viscosities but does give derivatives */
    p = T = 555; grad_p.setConstant(555); grad_T.setConstant(555);
    suzerain::rholut::p_T(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T);

    BOOST_CHECK_CLOSE(p,
            4441.1984111498830681828969449447964342868418205971441637701L,
            close_enough);
    BOOST_CHECK_CLOSE(grad_p(0), grad_p_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_p(1), grad_p_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_p(2), grad_p_ans(2), close_enough);
    BOOST_CHECK_CLOSE(T,
            45.055636055143743363613040114696980466792183873003314822257L,
            close_enough);
    BOOST_CHECK_CLOSE(grad_T(0), grad_T_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_T(1), grad_T_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_T(2), grad_T_ans(2), close_enough);

    /* Variant computing pressure from density and temperature */
    p = 555;
    suzerain::rholut::p(gamma, rho, T, p);

    BOOST_CHECK_CLOSE(p,
            4441.1984111498830681828969449447964342868418205971441637701L,
            close_enough);

    /* Variant providing only pressure */
    p = 555;
    suzerain::rholut::p(alpha, beta, gamma, Ma, rho, m, e, p);

    BOOST_CHECK_CLOSE(p,
            4441.1984111498830681828969449447964342868418205971441637701L,
            close_enough);
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_grad_p_and_div_grad_T )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);


    /* Expected results found using test_rholut.sage */
    {
        const double div_grad_p = suzerain::rholut::div_grad_p(
                gamma, Ma,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, grad_e, div_grad_e);

        BOOST_CHECK_CLOSE(div_grad_p,
                106145.68434338822542886306883326636335276520086897118283389L,
                close_enough);

        const double div_grad_T = suzerain::rholut::div_grad_T(
                gamma,
                rho, grad_rho, div_grad_rho,
                p, grad_p, div_grad_p);

        BOOST_CHECK_CLOSE(div_grad_T,
                1024.4624544088198099733042779277734936507131500839456171059L,
                close_enough);
    }

    /* With zero refcoeffs, explicit operator matches the full one */
    {
        const double mu_div_grad_T = suzerain::rholut::explicit_mu_div_grad_T(
                gamma, Ma, mu,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, div_grad_e,
                p, grad_p,
                0, Vector3r::Zero(), 0);

        BOOST_CHECK_CLOSE(mu_div_grad_T,
            12971.657203680849666061961791022666712780901207909760791631L,
            close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from the full one */
    {
        const double refcoeff_div_grad_e = 34;
        const double mu_div_grad_T = suzerain::rholut::explicit_mu_div_grad_T(
                gamma, Ma, mu,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, div_grad_e,
                p, grad_p,
                refcoeff_div_grad_e, Vector3r::Zero(), 0);

        BOOST_CHECK_CLOSE(mu_div_grad_T,
            -117985.46279631915033393803820897733328721909879209023920837L,
            close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from the full one */
    {
        const Vector3r refcoeff_div_grad_m(8, 13, 21);
        const double mu_div_grad_T = suzerain::rholut::explicit_mu_div_grad_T(
                gamma, Ma, mu,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, div_grad_e,
                p, grad_p,
                0, refcoeff_div_grad_m, 0);

        BOOST_CHECK_CLOSE(mu_div_grad_T,
            1.3738494612641828659176265335876502663867247605783926741805e7L,
            close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from the full one */
    {
        const double refcoeff_div_grad_rho = 987;
        const double mu_div_grad_T = suzerain::rholut::explicit_mu_div_grad_T(
                gamma, Ma, mu,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, div_grad_e,
                p, grad_p,
                0, Vector3r::Zero(), refcoeff_div_grad_rho);

        BOOST_CHECK_CLOSE(mu_div_grad_T,
            -72699.942796319150333938038208977333287219098792090239208369L,
            close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double mu     = 4181.0;
        const double rho    = 67.0;
        const Vector3r m(144.0, 233.0, 377.0);
        const double p      = 55.0;
        const double e      = p/(gamma-1)+Ma*Ma*m.squaredNorm()/(2*rho);

        BOOST_CHECK_CLOSE(
            suzerain::rholut::explicit_mu_div_grad_T_refcoeff_div_grad_e(
                    mu, rho),
                mu/rho,
                close_enough);
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(
                suzerain::rholut::explicit_mu_div_grad_T_refcoeff_div_grad_m(
                        mu, rho, m)[i],
                    mu/rho/rho*m[i],
                    close_enough);
        }
        BOOST_CHECK_CLOSE(
            suzerain::rholut::explicit_mu_div_grad_T_refcoeff_div_grad_rho(
                    gamma, mu, rho, e, p),
                mu/rho/rho*((gamma-1)*e-2*p),
                close_enough);
    }
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_grad_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const Matrix3r grad_u
        = suzerain::rholut::grad_u(
                rho, grad_rho, m, grad_m);

    /* Expected results found using test_rholut.sage */
    Matrix3r ans;
    ans(0,0) =   0.13354624933259255905305717414904148571263050381552617302315L;
    ans(0,1) = - 0.95458058163483407021971942603699868798014583786167988660320L;
    ans(0,2) =   0.25857485039372819387032260745049292527417870255089489945870L;
    ans(1,0) = - 0.47107453715872036912739336512687808059078448566164237217442L;
    ans(1,1) =   2.5594178135586821836692525516520979567520812757446061058222L;
    ans(1,2) =   3.8612842725703801936593326729720018775665984897457448852524L;
    ans(2,0) =   6.9043291992970668207962385482326888609251394067150989201691L;
    ans(2,1) =   2.0133656136592902780233887313332788854235250587761338439305L;
    ans(2,2) = -10.545791208924418168631572414966656671758514672569415551028L;

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;
    BOOST_CHECK_CLOSE(grad_u(0,0), ans(0,0), close_enough);
    BOOST_CHECK_CLOSE(grad_u(0,1), ans(0,1), close_enough);
    BOOST_CHECK_CLOSE(grad_u(0,2), ans(0,2), close_enough);
    BOOST_CHECK_CLOSE(grad_u(1,0), ans(1,0), close_enough);
    BOOST_CHECK_CLOSE(grad_u(1,1), ans(1,1), close_enough);
    BOOST_CHECK_CLOSE(grad_u(1,2), ans(1,2), close_enough);
    BOOST_CHECK_CLOSE(grad_u(2,0), ans(2,0), close_enough);
    BOOST_CHECK_CLOSE(grad_u(2,1), ans(2,1), close_enough);
    BOOST_CHECK_CLOSE(grad_u(2,2), ans(2,2), close_enough);
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double div_u = suzerain::rholut::div_u(
            rho, grad_rho, m, div_m);

    const double ans = -7.8528271460331434259092626891655172292938028930092832721830L;
    const double close_enough = std::numeric_limits<double>::epsilon() *1e2;
    BOOST_CHECK_CLOSE(div_u, ans, close_enough);
}

// Checks derived formula and computed result against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_grad_div_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    /* Expected results found using test_rholut.sage */
    {
        const Vector3r grad_div_u
            = suzerain::rholut::grad_div_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);

        const Vector3r ans(
            3.5808667611324763961641377901365615487146733886413982972779L,
            - 11.378277959392865631969701464497046416747061061234520455900L,
            237.13623643835318300159939437852909295361695206836632310924L);

        BOOST_CHECK_CLOSE(grad_div_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(grad_div_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(grad_div_u(2), ans(2), close_enough);
    }

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* With zero refcoeffs, explicit operator matches the full one */
    {
        const Vector3r mu_plus_lambda_grad_div_u
            = suzerain::rholut::explicit_mu_plus_lambda_grad_div_u(
                    mu, lambda, rho, grad_rho, grad_grad_rho, m,
                    div_m, grad_m, grad_div_m, 0, Vector3r::Zero());

        const Vector3r ans(
            241.81670907217979012053789700043342987966984221576324700315L,
            -768.37757855551448923175290916353582295349427995015552075162L,
            16013.861481723930920490790400651244064018877204543331179645L);

        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const double refcoeff_grad_div_m = 89.0;

        const Vector3r mu_plus_lambda_grad_div_u
            = suzerain::rholut::explicit_mu_plus_lambda_grad_div_u(
                    mu, lambda, rho, grad_rho, grad_grad_rho, m,
                    div_m, grad_m, grad_div_m,
                    refcoeff_grad_div_m, Vector3r::Zero());

        const Vector3r ans(
              20297.646608160654352669729166079334506758500652672642756083L,
              180150.11766316746189821312928781893696704411346016250571166L,
              -2.8050198934000174271785198217008585932496795024719430479111e6L);

        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const Vector3r refcoeff_grad_grad_rho(144.0, 233.0, 377.0);

        const Vector3r mu_plus_lambda_grad_div_u
            = suzerain::rholut::explicit_mu_plus_lambda_grad_div_u(
                    mu, lambda, rho, grad_rho, grad_grad_rho, m,
                    div_m, grad_m, grad_div_m,
                    0, refcoeff_grad_grad_rho);

        const Vector3r ans(
            55526.816709072179790120537897000433429879669842215763247003L,
            33405.622421444485510768247090836464177046505720049844479248L,
            45325.861481723930920490790400651244064018877204543331179645L);

        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        BOOST_CHECK_CLOSE(
            suzerain::rholut
                ::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m(
                    mu, lambda, rho),
                (mu+lambda)/rho,
                close_enough);
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(
                suzerain::rholut
                    ::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho(
                        mu, lambda, rho, m)[i],
                    (mu+lambda)/rho/rho*m[i],
                    close_enough);
        }
    }
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_grad_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    /* Expected results found using test_rholut.sage */
    {
      const Vector3r div_grad_u
          = suzerain::rholut::div_grad_u(
                  rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);
      const Vector3r ans(
        - 16.318446092843163297609832709693050751558008045050013765743L,
          23.204406985769508424700716673327465318352740571910767421546L,
         672.79795991861399600487930607996269111701083307478613865859L);

      BOOST_CHECK_CLOSE(div_grad_u(0), ans(0), close_enough);
      BOOST_CHECK_CLOSE(div_grad_u(1), ans(1), close_enough);
      BOOST_CHECK_CLOSE(div_grad_u(2), ans(2), close_enough);
    }

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* With zero refcoeffs, explicit operator matches the full one */
    {
        const Vector3r ans(
            -206.62278827510370976247563121644883265411981528196528710676L,
            293.81224440069430722715479229137439391998698236497510955195L,
            8518.9110306988081894321298782906014223456866685196821822478L);

        const Vector3r mu_div_grad_u
            = suzerain::rholut::explicit_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                0, Vector3r::Zero());

        BOOST_CHECK_CLOSE(mu_div_grad_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(2), ans(2), close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const double refcoeff_div_grad_m = 2584.0;

        const Vector3r ans(
            6.0231588101609568919781611528433660166131516876938717150225e6L,
            -1.0560877564142245696202926395927519051466501040860560444641e7L,
            -2.4194253168223783705833273397742346446478343908855596563795e8L);

        const Vector3r mu_div_grad_u
            = suzerain::rholut::explicit_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                refcoeff_div_grad_m, Vector3r::Zero());

        BOOST_CHECK_CLOSE(mu_div_grad_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(2), ans(2), close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const Vector3r refcoeff_div_grad_rho(610.0, 987.0, 1597.0);
        const Vector3r ans(
            37613.377211724896290237524368783551167345880184718034712893L,
            61487.812244400694307227154792291374393919986982364975109552L,
            107532.91103069880818943212987829060142234568666851968218225L);

        const Vector3r mu_div_grad_u
            = suzerain::rholut::explicit_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                0, refcoeff_div_grad_rho);

        BOOST_CHECK_CLOSE(mu_div_grad_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(2), ans(2), close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const Vector3r m(3.0, 5.0, 7.0);

        BOOST_CHECK_CLOSE(suzerain::rholut
            ::explicit_mu_div_grad_u_refcoeff_div_grad_m(
                mu, rho),
            mu/rho, close_enough);
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(suzerain::rholut
                ::explicit_mu_div_grad_u_refcoeff_div_grad_rho(
                    mu, rho, m)[i],
                (mu/rho)/rho * m[i], close_enough);
        }
    }
}

// Checks computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_tau_and_div_tau )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
    const double          div_u  = suzerain::rholut::div_u(rho, grad_rho, m, div_m);
    const Matrix3r grad_u = suzerain::rholut::grad_u(rho, grad_rho, m, grad_m);
    const Vector3r grad_div_u = suzerain::rholut::grad_div_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);
    const Vector3r div_grad_u = suzerain::rholut::div_grad_u(
                rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e4;
    {
        using namespace suzerain;
        const Matrix3r tau = rholut::tau(mu, lambda, div_u, grad_u);

        /* Expected results found using test_rholut.sage */
        Matrix3r ans;
        ans(0,0) = -427.48938267676174372726167773557639521947426572548406990997L;
        ans(0,1) = -18.051524887102470041909905956029142873764292311782332127087L;
        ans(0,2) = 90.696087021621572471230415659234366391108287726762864275979L;
        ans(1,0) = ans(0,1);
        ans(1,1) = -366.05702033712541523451517498439074337943946797646351939919L;
        ans(1,2) = 74.384321443764902359727211537482893239920338590287217190636L;
        ans(2,0) = ans(0,2);
        ans(2,1) = ans(1,2);
        ans(2,2) = -697.93112326915506200245437123956619312724596885631534519240L;

        // Check for acceptability of entries
        BOOST_CHECK_CLOSE(tau(0,0), ans(0,0), close_enough);
        BOOST_CHECK_CLOSE(tau(0,1), ans(0,1), close_enough);
        BOOST_CHECK_CLOSE(tau(0,2), ans(0,2), close_enough);
        BOOST_CHECK_CLOSE(tau(1,0), ans(1,0), close_enough);
        BOOST_CHECK_CLOSE(tau(1,1), ans(1,1), close_enough);
        BOOST_CHECK_CLOSE(tau(1,2), ans(1,2), close_enough);
        BOOST_CHECK_CLOSE(tau(2,0), ans(2,0), close_enough);
        BOOST_CHECK_CLOSE(tau(2,1), ans(2,1), close_enough);
        BOOST_CHECK_CLOSE(tau(2,2), ans(2,2), close_enough);

        // Check for exact symmetry of the computed result
        BOOST_CHECK_EQUAL(tau(0,1), tau(1,0));
        BOOST_CHECK_EQUAL(tau(0,2), tau(2,0));
        BOOST_CHECK_EQUAL(tau(2,1), tau(1,2));
    }

    {
        using namespace suzerain;
        const Vector3r div_tau = rholut::div_tau(
                mu, grad_mu, lambda, grad_lambda,
                div_u, grad_u, div_grad_u, grad_div_u);

        /* Expected results found using test_rholut.sage */
        Vector3r ans;
        ans(0) = -182.52979655440578315227418884504311513720860612160033737778L;
        ans(1) = -579.83808129185098230400222995823008543250988843416411315663L;
        ans(2) = 24946.768752288838607765076090999038104032108245970482384226L;

        BOOST_CHECK_CLOSE(div_tau(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_tau(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_tau(2), ans(2), close_enough);
    }
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_e_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    using namespace suzerain;
    const Vector3r u = rholut::u(rho, m);
    const double div_u = rholut::div_u(
            rho, grad_rho, m, div_m);
    const double div_e_u = rholut::div_e_u(e, grad_e, u, div_u);

    /* Expected results found using test_rholut.sage */
    const double ans
        = -90849.212502207575911979472073826868082163956391101506206117L;

    const double close_enough = std::numeric_limits<double>::epsilon() * 1e2;
    BOOST_CHECK_CLOSE(div_e_u, ans, close_enough);
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_p_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon()*1e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    using namespace suzerain;
    rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* Expected results found using test_rholut.sage */
    {
        const Vector3r u = rholut::u(rho, m);
        const double div_u = rholut::div_u(
                rho, grad_rho, m, div_m);
        const double div_p_u = rholut::div_p_u(p, grad_p, u, div_u);
        const double ans
            = -33816.441337235609376147444530387100036324237536209394475325L;

        BOOST_CHECK_CLOSE(div_p_u, ans, close_enough);
    }

    /* Explicit operator should coincide when refcoeffs are zero */
    {
        const double explicit_div_p_u = rholut::explicit_div_p_u(
                rho, grad_rho, m, div_m, p, grad_p,
                0, Vector3r::Zero());
        const double ans
            = -33816.441337235609376147444530387100036324237536209394475325L;
        BOOST_CHECK_CLOSE(explicit_div_p_u, ans, close_enough);
    }

    /* Explicit operator differs for nonzero refcoeff */
    {
        const double refcoeff_div_m = 77.0;
        const double explicit_div_p_u = rholut::explicit_div_p_u(
                rho, grad_rho, m, div_m, p, grad_p,
                refcoeff_div_m, Vector3r::Zero());
        const double ans
            = 51726.934749018138228480007905533818427673945152346041105157L;
        BOOST_CHECK_CLOSE(explicit_div_p_u, ans, close_enough);
    }

    /* Explicit operator differs for nonzero refcoeff */
    {
        const Vector3r refcoeff_grad_rho(77.0, 88.0, 99.0);
        const double explicit_div_p_u = rholut::explicit_div_p_u(
                rho, grad_rho, m, div_m, p, grad_p,
                0, refcoeff_grad_rho);
        const double ans
            = -7086.4413372356093761474445303871000363242375362093944753249L;
        BOOST_CHECK_CLOSE(explicit_div_p_u, ans, close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double rho    = 67.0;
        const Vector3r m(144.0, 233.0, 377.0);
        const double p      = 55.0;

        BOOST_CHECK_CLOSE(rholut::explicit_div_p_u_refcoeff_div_m(
            rho, p), p/rho, close_enough);
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(rholut::explicit_div_p_u_refcoeff_grad_rho(
                rho, m, p)[i], p/rho/rho*m[i], close_enough);
        }
    }
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_grad_p )
{
    const double close_enough = std::numeric_limits<double>::epsilon()*1.1e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double gamma = 1.4;
    const double Ma    = 3.5;

    using namespace suzerain;

    /* Should recover the full operator when no refcoeffs are used */
    /* and the result is appropriately augmented by a grad_e term. */
    {
        const Vector3r ans = (gamma - 1)*grad_e
            + rholut::explicit_grad_p(gamma, Ma, rho, grad_rho, m, grad_m,
                                     0, Vector3r::Zero());

        const Vector3r grad_p_ans(
                7432.6318085149443001747520526182772856422787520432873299439L,
                3876.3726619344082666965507166633629625480052934837862921105L,
                -1052.4623117340587311081836624513934718452258172103981784331L);
        BOOST_CHECK_CLOSE(ans(0), grad_p_ans(0), close_enough);
        BOOST_CHECK_CLOSE(ans(1), grad_p_ans(1), close_enough);
        BOOST_CHECK_CLOSE(ans(2), grad_p_ans(2), close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    {
        const double refcoeff_grad_rho = 7;
        const Vector3r ans = (gamma - 1)*grad_e
            + (gamma - 1)/2 * Ma * Ma * refcoeff_grad_rho * grad_rho
            + rholut::explicit_grad_p(gamma, Ma, rho, grad_rho, m, grad_m,
                                     refcoeff_grad_rho,
                                     Vector3r::Zero());

        const Vector3r grad_p_ans(
                7432.6318085149443001747520526182772856422787520432873299439L,
                3876.3726619344082666965507166633629625480052934837862921105L,
                -1052.4623117340587311081836624513934718452258172103981784331L);
        BOOST_CHECK_CLOSE(ans(0), grad_p_ans(0), close_enough);
        BOOST_CHECK_CLOSE(ans(1), grad_p_ans(1), close_enough);
        BOOST_CHECK_CLOSE(ans(2), grad_p_ans(2), close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    {
        const Vector3r refcoeff_grad_m(2, 3, 5);
        const Vector3r ans = (gamma - 1)*grad_e
            - (gamma - 1) * Ma * Ma *grad_m.transpose()*refcoeff_grad_m
            + rholut::explicit_grad_p(gamma, Ma, rho, grad_rho, m, grad_m,
                                     0, refcoeff_grad_m);

        const Vector3r grad_p_ans(
                7432.6318085149443001747520526182772856422787520432873299439L,
                3876.3726619344082666965507166633629625480052934837862921105L,
                -1052.4623117340587311081836624513934718452258172103981784331L);
        BOOST_CHECK_CLOSE(ans(0), grad_p_ans(0),     close_enough);
        BOOST_CHECK_CLOSE(ans(1), grad_p_ans(1),     close_enough);
        BOOST_CHECK_CLOSE(ans(2), grad_p_ans(2), 2.5*close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double rho    = 67.0;
        const Vector3r m(144.0, 233.0, 377.0);

        BOOST_CHECK_CLOSE(rholut::explicit_grad_p_refcoeff_grad_rho(rho, m),
            217154.0L/4489.0L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_grad_p_refcoeff_grad_m(rho, m)[0],
            144.0L/67.0L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_grad_p_refcoeff_grad_m(rho, m)[1],
            233.0L/67.0L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_grad_p_refcoeff_grad_m(rho, m)[2],
            377.0L/67.0L,
            close_enough);
    }
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_e_plus_p_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon()*5.0e2;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    using namespace suzerain;
    rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* Should recover the full operator when no refcoeffs used */
    {
        const double ans
            = -124665.65383944318528812691660421396811848819392731090068144L;

        const double div_e_plus_p_u = rholut::explicit_div_e_plus_p_u(
                gamma, Ma, rho, grad_rho, m, div_m, grad_m, e, grad_e, p,
                0, Vector3r::Zero(), Vector3r::Zero());
        BOOST_CHECK_CLOSE(div_e_plus_p_u, ans, close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    {
        const double refcoeff_div_m = 55.0;
        const double div_e_plus_p_u = rholut::explicit_div_e_plus_p_u(
                gamma, Ma, rho, grad_rho, m, div_m, grad_m, e, grad_e, p,
                refcoeff_div_m, Vector3r::Zero(), Vector3r::Zero());
        const double ans
            = -63563.242349261936999107307721413312072775206292628446695383L;
        BOOST_CHECK_CLOSE(div_e_plus_p_u, ans, close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    {
        const Vector3r refcoeff_grad_rho(55.0, 66.0, 77.0);
        const double div_e_plus_p_u = rholut::explicit_div_e_plus_p_u(
                gamma, Ma, rho, grad_rho, m, div_m, grad_m, e, grad_e, p,
                0, refcoeff_grad_rho, Vector3r::Zero());
        const double ans
            = -144509.65383944318528812691660421396811848819392731090068144L;
        BOOST_CHECK_CLOSE(div_e_plus_p_u, ans, close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    /* Odd looking 1/gamma factors due to removing gamma from refcoeff */
    {
        const Vector3r refcoeff_grad_e(55/gamma, 66/gamma, 77/gamma);
        const double div_e_plus_p_u = rholut::explicit_div_e_plus_p_u(
                gamma, Ma, rho, grad_rho, m, div_m, grad_m, e, grad_e, p,
                0, Vector3r::Zero(), refcoeff_grad_e);
        const double ans
            = -1.7853136538394431852881269166042139681184881939273109006814e6L;
        BOOST_CHECK_CLOSE(div_e_plus_p_u, ans, close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    /* Odd looking 1/gamma factors due to removing gamma from refcoeff */
    {
        BOOST_CHECK_CLOSE(rholut::explicit_div_e_plus_p_u_refcoeff_div_m(
            rho, e, p),
            114.26955370398466560143539493966043062887674127089522487304L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_div_e_plus_p_u_refcoeff_grad_rho(
            gamma, rho, m, e, p)[0],
            -14.021767904044756739151691226277061209937266437049794157090L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_div_e_plus_p_u_refcoeff_grad_rho(
            gamma, rho, m, e, p)[1],
            10.993554874700935706806005403121127980699831032811586707400L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_div_e_plus_p_u_refcoeff_grad_rho(
            gamma, rho, m, e, p)[2],
            55.846741669840866302567971498244085862638828473481931303727L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_div_e_plus_p_u_refcoeff_grad_e(
            rho, m)[0], (1 / gamma) *
            0.17277706184595877178290335962894047332330723544670029709746L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_div_e_plus_p_u_refcoeff_grad_e(
            rho, m)[1], (1 / gamma) *
            -0.13546323997740892943292160465580589068130661103855157024536L,
            close_enough);
        BOOST_CHECK_CLOSE(rholut::explicit_div_e_plus_p_u_refcoeff_grad_e(
            rho, m)[2], (1 / gamma) *
            -0.68814688742651284707848084677245457022236021675947361667360L,
            close_enough);
    }
}


// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_tau_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    using namespace suzerain;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    const Vector3r u = rholut::u(rho, m);
    const double div_u = rholut::div_u(
            rho, grad_rho, m, div_m);
    const Matrix3r grad_u = rholut::grad_u(
            rho, grad_rho, m, grad_m);
    const Vector3r div_grad_u = rholut::div_grad_u(
            rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);
    const Vector3r grad_div_u = rholut::grad_div_u(
            rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);
    const Matrix3r tau = rholut::tau(mu, lambda, div_u, grad_u);
    const Vector3r div_tau = rholut::div_tau(
            mu, grad_mu, lambda, grad_lambda,
            div_u, grad_u, div_grad_u, grad_div_u);

    const double div_tau_u = rholut::div_tau_u<double>(
            u, grad_u, tau, div_tau);

    /* Expected results found using test_rholut.sage */
    const double ans
        = -4749.9760126048652783688963154973108497305613678261015403062L;

    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e3;
    BOOST_CHECK_CLOSE(div_tau_u, ans, close_enough);
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_mu_grad_T )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double Ma    = 3.5;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholut::p_T_mu_lambda(
            alpha, beta, gamma, Ma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    const double div_grad_p = suzerain::rholut::div_grad_p(
        gamma, Ma,
        rho, grad_rho, div_grad_rho,
        m, grad_m, div_grad_m,
        e, grad_e, div_grad_e);
    const double div_grad_T = suzerain::rholut::div_grad_T(
        gamma,
        rho, grad_rho, div_grad_rho,
        p, grad_p, div_grad_p);

    const double div_mu_grad_T = suzerain::rholut::div_mu_grad_T(
            grad_T, div_grad_T, mu, grad_mu);

    /* Expected results found using test_rholut.sage */
    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e3;
    BOOST_CHECK_CLOSE(div_mu_grad_T,
            13360.858914707143946925379545454624752172684036897854744501L,
            close_enough);
}

// Checks derived formula and computation against rholut_test_data()
BOOST_AUTO_TEST_CASE( rholut_div_rho_inverse_m_outer_m )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    /* Expected results found using test_rholut.sage */
    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e3;
    const Vector3r ans(
            - 139.62394416218589511382305374934119093270491143228219517924L,
            - 196.62019324654648600809920433959614291380888927326388185922L,
             1352.1114315042037165058005031972730396574030986356606221464L );

    /* Compute using div_rho_inverse_m_outer_m */
    {
        const Vector3r div_rho_inverse_m_outer_m
            = suzerain::rholut::div_rho_inverse_m_outer_m(
                    rho, grad_rho, m, div_m, grad_m);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(2), ans(2), close_enough);
    }

    /* Compute using div_u_outer_m */
    {
        const Vector3r u = suzerain::rholut::u(
                rho, m);
        const double div_u = suzerain::rholut::div_u(
                rho, grad_rho, m, div_m);

        const Vector3r div_u_outer_m
            = suzerain::rholut::div_u_outer_m(
                    m, grad_m, u, div_u);

        BOOST_CHECK_CLOSE(div_u_outer_m(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_u_outer_m(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_u_outer_m(2), ans(2), close_enough);
    }
}
