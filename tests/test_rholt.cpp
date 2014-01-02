//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

#include <suzerain/rholt.hpp>

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
// Floating point values evaluated using test_rholt.sage
// to 200 bits of precision.
static
void rholt_test_data(
        double   &rho,
        Vector3r &grad_rho,
        double   &div_grad_rho,
        Matrix3r &grad_grad_rho,
        Vector3r &m,
        double   &div_m,
        Matrix3r &grad_m,
        Vector3r &div_grad_m,
        Vector3r &grad_div_m,
        Vector3r &curl_curl_m,
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

    // Uses curl_curl_m = grad_div_m - div_grad_m identity
    curl_curl_m(0) = grad_div_m(0) - div_grad_m(0);
    curl_curl_m(1) = grad_div_m(1) - div_grad_m(1);
    curl_curl_m(2) = grad_div_m(2) - div_grad_m(2);

    e = 11328.;

    grad_e(0) = 13194.;
    grad_e(1) =  7542.;
    grad_e(2) =  5678.;

    div_grad_e = 6878.;
}

// Helper macro for loading test data via rholt_test_data
#define ADD_TEST_DATA_TO_ENCLOSING_SCOPE                        \
    double          rho;                                        \
    Vector3r grad_rho;                                          \
    double          div_grad_rho;                               \
    Matrix3r grad_grad_rho;                                     \
    Vector3r m;                                                 \
    double          div_m;                                      \
    Matrix3r grad_m;                                            \
    Vector3r div_grad_m;                                        \
    Vector3r grad_div_m;                                        \
    Vector3r curl_curl_m;                                       \
    double          e;                                          \
    Vector3r grad_e;                                            \
    double          div_grad_e;                                 \
                                                                \
    rholt_test_data(                                            \
        rho, grad_rho, div_grad_rho, grad_grad_rho,             \
        m, div_m, grad_m, div_grad_m, grad_div_m, curl_curl_m,  \
        e, grad_e, div_grad_e)

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const Vector3r u = suzerain::rholt::u(rho, m);

    /* Expected results found using test_rholt.sage */
    const Vector3r ans(
             0.12341218703282769413064525687781462380236231103335735506961L,
            -0.096759457126720663880658289039861350486647579313251121603828L,
            -0.49153349101893774791320060483746755015882872625676686905257L);

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e2;
    BOOST_CHECK_CLOSE(u(0), ans(0), close_enough);
    BOOST_CHECK_CLOSE(u(1), ans(1), close_enough);
    BOOST_CHECK_CLOSE(u(2), ans(2), close_enough);
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_p_T_mu_lambda )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    /* Variant that provides derivatives */
    suzerain::rholt::p_T_mu_lambda(
            alpha, beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* Expected results found using test_rholt.sage */
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e4;
    BOOST_CHECK_CLOSE(p,
            4523.8529315224394491652244924378285975908038303556842666509L,
            close_enough);
    BOOST_CHECK_CLOSE(T,
            45.894160174865327745154451372557681424834241757231579516748L,
            close_enough);
    BOOST_CHECK_CLOSE(mu,
            12.818531787494714755711740676149526369900526768964830775976L,
            close_enough);
    BOOST_CHECK_CLOSE(lambda,
            55.546971079143763941417542929981280936235615998847600029231L,
            close_enough);

    const Vector3r grad_p_ans(
            5453.5209639604035869584122627487425816813409506343630929293L,
            3086.9691968926047335425739267829758568639555587282859345292L,
            1999.8806276135462562543621611330864774644251338499615213427L);
    BOOST_CHECK_CLOSE(grad_p(0), grad_p_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_p(1), grad_p_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_p(2), grad_p_ans(2), close_enough);

    const Vector3r grad_T_ans(
             5.4406182848896076809319526229317927581792830964012420477440L,
             2.3838039162055298052983060006061443162968025314525571766942L,
            -4.9864006857304358686639947733917588394000520736206805318339L);
    BOOST_CHECK_CLOSE(grad_T(0), grad_T_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_T(1), grad_T_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_T(2), grad_T_ans(2), close_enough);

    const Vector3r grad_mu_ans(
             1.0130662690381109382304916850142818302534744683039668368759L,
             0.44387442989997855519790991486035061354804176216270679061604L,
            -0.92848901983288259118625731085663474679900828620195051533844L);
    BOOST_CHECK_CLOSE(grad_mu(0), grad_mu_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_mu(1), grad_mu_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_mu(2), grad_mu_ans(2), close_enough);

    const Vector3r grad_lambda_ans(
            4.3899538324984807323321306350618879310983893626505229597954L,
            1.9234558628999070725242762977281859920415143027050627593362L,
            -4.0234524192758245618071150137120839027957025735417855664666L);
    BOOST_CHECK_CLOSE(grad_lambda(0), grad_lambda_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_lambda(1), grad_lambda_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_lambda(2), grad_lambda_ans(2), close_enough);

    /* Check computation of grad_e from gradient of internal, kinetic */
    const Vector3r alt_grad_e
        = suzerain::rholt::energy_internal_gradient(gamma, grad_p)
        + suzerain::rholt::energy_kinetic_gradient(rho, grad_rho, m, grad_m);
    BOOST_CHECK_CLOSE(grad_e(0), alt_grad_e(0), close_enough);
    BOOST_CHECK_CLOSE(grad_e(1), alt_grad_e(1), close_enough);
    BOOST_CHECK_CLOSE(grad_e(2), alt_grad_e(2), close_enough);

    /* Variant that does not provide derivatives */
    p = T = mu = lambda = 555;
    suzerain::rholt::p_T_mu_lambda(
            alpha, beta, gamma, rho, m, e, p, T, mu, lambda);

    BOOST_CHECK_CLOSE(p,
            4523.8529315224394491652244924378285975908038303556842666509L,
            close_enough);
    BOOST_CHECK_CLOSE(T,
            45.894160174865327745154451372557681424834241757231579516748L,
            close_enough);
    BOOST_CHECK_CLOSE(mu,
            12.818531787494714755711740676149526369900526768964830775976L,
            close_enough);
    BOOST_CHECK_CLOSE(lambda,
            55.546971079143763941417542929981280936235615998847600029231L,
            close_enough);

    /* Variant that does not provide viscosities */
    p = T = 555;
    suzerain::rholt::p_T(alpha, beta, gamma, rho, m, e, p, T);

    BOOST_CHECK_CLOSE(p,
            4523.8529315224394491652244924378285975908038303556842666509L,
            close_enough);
    BOOST_CHECK_CLOSE(T,
            45.894160174865327745154451372557681424834241757231579516748L,
            close_enough);

    /* Variant that does not provide viscosities but does derivatives */
    p = T = 555; grad_p.setConstant(555); grad_T.setConstant(555);
    suzerain::rholt::p_T(
            alpha, beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T);

    BOOST_CHECK_CLOSE(p,
            4523.8529315224394491652244924378285975908038303556842666509L,
            close_enough);
    BOOST_CHECK_CLOSE(grad_p(0), grad_p_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_p(1), grad_p_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_p(2), grad_p_ans(2), close_enough);
    BOOST_CHECK_CLOSE(T,
            45.894160174865327745154451372557681424834241757231579516748L,
            close_enough);
    BOOST_CHECK_CLOSE(grad_T(0), grad_T_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_T(1), grad_T_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_T(2), grad_T_ans(2), close_enough);

    /* Variant computing pressure from density and temperature */
    p = 555;
    suzerain::rholt::p(gamma, rho, T, p);

    BOOST_CHECK_CLOSE(p,
            4523.8529315224394491652244924378285975908038303556842666509L,
            close_enough);

    /* Variant providing only pressure */
    p = 555;
    suzerain::rholt::p(alpha, beta, gamma, rho, m, e, p);

    BOOST_CHECK_CLOSE(p,
            4523.8529315224394491652244924378285975908038303556842666509L,
            close_enough);
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_grad_p_and_div_grad_T )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholt::p_T_mu_lambda(
            alpha, beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);


    /* Expected results found using test_rholt.sage */
    {
        const double div_grad_p = suzerain::rholt::div_grad_p(
                gamma,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, grad_e, div_grad_e);

        BOOST_CHECK_CLOSE(div_grad_p,
                11191.566068848019028203711722987453542724808010087855921969L,
                close_enough);

        const double div_grad_T = suzerain::rholt::div_grad_T(
                gamma,
                rho, grad_rho, div_grad_rho,
                p, grad_p, div_grad_p);

        BOOST_CHECK_CLOSE(div_grad_T,
                83.577681904999696238558381171408689250539040369361349554271L,
                close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from the full one */
    {
        const double mu = 4184.0;
        const double refcoeff_div_grad_e = 34;
        const double mu_div_grad_T = suzerain::rholt::explicit_mu_div_grad_T(
                gamma, mu,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, div_grad_e,
                p, grad_p,
                refcoeff_div_grad_e, Vector3r::Zero(), 0);

        BOOST_CHECK_CLOSE(mu_div_grad_T,
            218731.90109051878978113376152599227382356957602236358637891L,
            close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from the full one */
    {
        const double mu = 4184.0;
        const Vector3r refcoeff_div_grad_m(8, 13, 21);
        const double mu_div_grad_T = suzerain::rholt::explicit_mu_div_grad_T(
                gamma, mu,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, div_grad_e,
                p, grad_p,
                0, refcoeff_div_grad_m, 0);

        BOOST_CHECK_CLOSE(mu_div_grad_T,
            1.4701398745956734774018503638258185883274160880096467829055e6L,
            close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from the full one */
    {
        const double mu = 4184.0;
        const double refcoeff_div_grad_rho = 987;
        const double mu_div_grad_T = suzerain::rholt::explicit_mu_div_grad_T(
                gamma, mu,
                rho, grad_rho, div_grad_rho,
                m, grad_m, div_grad_m,
                e, div_grad_e,
                p, grad_p,
                0, Vector3r::Zero(), refcoeff_div_grad_rho);

        BOOST_CHECK_CLOSE(mu_div_grad_T,
            264017.42109051872109609388363445174647981957602236358637891L,
            close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double mu     = 4181.0;
        const double rho    = 67.0;
        const Vector3r m(144.0, 233.0, 377.0);
        const double p      = 55.0;
        const double e      = p/(gamma-1)+m.squaredNorm()/(2*rho);

        BOOST_CHECK_CLOSE(
            suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_e(
                    mu, rho),
                mu/rho,
                close_enough);
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(
                suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_m(
                        mu, rho, m)[i],
                    mu/rho/rho*m[i],
                    close_enough);
        }
        BOOST_CHECK_CLOSE(
            suzerain::rholt::explicit_mu_div_grad_T_refcoeff_div_grad_rho(
                    gamma, mu, rho, e, p),
                mu/rho/rho*((gamma-1)*e-2*p),
                close_enough);
    }
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_grad_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const Matrix3r grad_u = suzerain::rholt::grad_u(rho, grad_rho, m, grad_m);

    /* Expected results found using test_rholt.sage */
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

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double div_u = suzerain::rholt::div_u(
            rho, grad_rho, m, div_m);

    const double ans = -7.8528271460331434259092626891655172292938028930092832721830L;
    const double close_enough = std::numeric_limits<double>::epsilon() *1e2;
    BOOST_CHECK_CLOSE(div_u, ans, close_enough);
}

// Checks derived formula and computed result against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_grad_div_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    /* Expected results found using test_rholt.sage */
    {
        const Vector3r grad_div_u
            = suzerain::rholt::grad_div_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);

        const Vector3r ans(
            3.5808667611324763961641377901365615487146733886413982972779L,
            - 11.378277959392865631969701464497046416747061061234520455900L,
            237.13623643835318300159939437852909295361695206836632310924L);

        BOOST_CHECK_CLOSE(grad_div_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(grad_div_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(grad_div_u(2), ans(2), close_enough);
    }

    /* With zero refcoeffs, explicit operator matches the full one */
    {
        const double mu     = 4181.0;
        const double lambda = 6765.0;

        const Vector3r mu_plus_lambda_grad_div_u
            = suzerain::rholt::explicit_mu_plus_lambda_grad_div_u(
                    mu, lambda, rho, grad_rho, grad_grad_rho, m,
                    div_m, grad_m, grad_div_m, 0, Vector3r::Zero());

        const Vector3r ans(
            (mu+lambda)
            *3.5808667611324763961641377901365615487146733886413982972779L,
            (mu+lambda)
            *-11.378277959392865631969701464497046416747061061234520455900L,
            (mu+lambda)
            *237.13623643835318300159939437852909295361695206836632310924L);

        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const double mu     = 4181.0;
        const double lambda = 6765.0;
        const double refcoeff_grad_div_m = 89.0;

        const Vector3r mu_plus_lambda_grad_div_u
            = suzerain::rholt::explicit_mu_plus_lambda_grad_div_u(
                    mu, lambda, rho, grad_rho, grad_grad_rho, m,
                    div_m, grad_m, grad_div_m,
                    refcoeff_grad_div_m, Vector3r::Zero());

        const Vector3r ans(
                59251.997466444561194961843519913703789109645722525625271084L,
                56371.864698208669179904529966597802712284277363839600322130L,
              -225340.51082752741696350364123413038584340722233614860633701L);

        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const double mu     = 4181.0;
        const double lambda = 6765.0;
        const Vector3r refcoeff_grad_grad_rho(144.0, 233.0, 377.0);

        const Vector3r mu_plus_lambda_grad_div_u
            = suzerain::rholt::explicit_mu_plus_lambda_grad_div_u(
                    mu, lambda, rho, grad_rho, grad_grad_rho, m,
                    div_m, grad_m, grad_div_m,
                    0, refcoeff_grad_grad_rho);

        const Vector3r ans(
             94481.167567356086632412652250834802712230814912068745762004L,
            -90372.630543514307207540352230384670077713330376273060910286L,
                 2.6250052440542139411355069708673794514702911573403377727538e6L);

        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double mu     = 4181.0;
        const double lambda = 6765.0;
        const double rho    = 67.0;
        const Vector3r m(144.0, 233.0, 377.0);

        BOOST_CHECK_CLOSE(
            suzerain::rholt
                ::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m(
                    mu, lambda, rho),
                (mu+lambda)/rho,
                close_enough);
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(
                suzerain::rholt
                    ::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho(
                        mu, lambda, rho, m)[i],
                    (mu+lambda)/rho/rho*m[i],
                    close_enough);
        }
    }
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_grad_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    /* Expected results found using test_rholt.sage */
    {
      const Vector3r div_grad_u = suzerain::rholt::div_grad_u(
                  rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);
      const Vector3r ans(
        - 16.318446092843163297609832709693050751558008045050013765743L,
          23.204406985769508424700716673327465318352740571910767421546L,
         672.79795991861399600487930607996269111701083307478613865859L);

      BOOST_CHECK_CLOSE(div_grad_u(0), ans(0), close_enough);
      BOOST_CHECK_CLOSE(div_grad_u(1), ans(1), close_enough);
      BOOST_CHECK_CLOSE(div_grad_u(2), ans(2), close_enough);
    }

    /* With zero refcoeffs, explicit operator matches the full one */
    {
        const double mu = 4181.0;

        const Vector3r ans(
            mu * - 16.318446092843163297609832709693050751558008045050013765743L,
            mu *   23.204406985769508424700716673327465318352740571910767421546L,
            mu *  672.79795991861399600487930607996269111701083307478613865859L);

        const Vector3r mu_div_grad_u
            = suzerain::rholt::explicit_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                0, Vector3r::Zero());

        BOOST_CHECK_CLOSE(mu_div_grad_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(2), ans(2), close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const double mu = 4181.0;
        const double refcoeff_div_grad_m = 2584.0;

        const Vector3r ans(
            5.9551380098350547299406169179153558202535417758727995727551e6L,
            -1.0464153750779144075786479854308628293364388219511766501161e7L,
            -2.3913808232284881074922576572858143105464556248213880447440e8L);

        const Vector3r mu_div_grad_u
            = suzerain::rholt::explicit_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                refcoeff_div_grad_m, Vector3r::Zero());

        BOOST_CHECK_CLOSE(mu_div_grad_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(2), ans(2), close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const double mu = 4181.0;
        const Vector3r refcoeff_div_grad_rho(610.0, 987.0, 1597.0);
        const Vector3r ans(
            - 30407.423114177265747306710559226645192264031636354107554571L,
            158211.62560750231472367369641118213249603280833115891858948L,
                2.9119822704197251172964003787203240115602222930856808457316e6L);

        const Vector3r mu_div_grad_u
            = suzerain::rholt::explicit_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                0, refcoeff_div_grad_rho);

        BOOST_CHECK_CLOSE(mu_div_grad_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(mu_div_grad_u(2), ans(2), close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double mu  = 4181.0;
        const double rho = 5678.0;
        const Vector3r m(3.0, 5.0, 7.0);

        BOOST_CHECK_CLOSE(suzerain::rholt
            ::explicit_mu_div_grad_u_refcoeff_div_grad_m(
                mu, rho),
            mu/rho, close_enough);
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(suzerain::rholt
                ::explicit_mu_div_grad_u_refcoeff_div_grad_rho(
                    mu, rho, m)[i],
                (mu/rho)/rho * m[i], close_enough);
        }
    }
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_u_dot_explicit_mu_div_grad_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;
    const double mu
        = 12.818531787494714755711740676149526369900526768964830775976L;

    /* With zero refcoeffs, explicit operator matches the full one */
    {
        const double ans
            = -4293.7193901585066389168068570176066143148258347128441480446L;

        const double u_dot_mu_div_grad_u
            = suzerain::rholt::explicit_u_dot_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                Vector3r::Zero(), 0);

        BOOST_CHECK_CLOSE(u_dot_mu_div_grad_u, ans, close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const Vector3r refcoeff_div_grad_m(3, 5, 7);
        const double ans
            = -673176.50373070234836923498008429608579249020553411631291690L;

        const double u_dot_mu_div_grad_u
            = suzerain::rholt::explicit_u_dot_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                refcoeff_div_grad_m, 0);

        BOOST_CHECK_CLOSE(u_dot_mu_div_grad_u, ans, close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const double refcoeff_div_grad_rho
            = 123;
        const double ans
            = 3332.2806098414933610831931429823933856851741652871558519554L;

        const double u_dot_mu_div_grad_u
            = suzerain::rholt::explicit_u_dot_mu_div_grad_u(
                mu, rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m,
                Vector3r::Zero(), refcoeff_div_grad_rho);

        BOOST_CHECK_CLOSE(u_dot_mu_div_grad_u, ans, close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double mu  = 4181.0;
        const double rho = 5678.0;
        const Vector3r m(3.0, 5.0, 7.0);

        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(suzerain::rholt
                ::explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_m(
                    mu, rho, m)[i],
                (mu/rho)/rho * m[i], close_enough);
        }
        BOOST_CHECK_CLOSE(suzerain::rholt
            ::explicit_u_dot_mu_div_grad_u_refcoeff_div_grad_rho(
                mu, rho, m),
            mu/(rho*rho*rho)*m.squaredNorm(), close_enough);
    }
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_u_dot_explicit_mu_plus_lambda_grad_div_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;
    const double mu
        = 12.818531787494714755711740676149526369900526768964830775976L;
    const double lambda
        = 55.546971079143763941417542929981280936235615998847600029231L;

    /* With zero refcoeffs, explicit operator matches the full one */
    {
        const double ans
            = -7863.2308357512977715114493722454403366670208895633912928559L;

        const double u_dot_mu_plus_lambda_grad_div_u
            = suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u(
                mu, lambda, rho, grad_rho, grad_grad_rho,
                m, div_m, grad_m, grad_div_m,
                Vector3r::Zero(), Matrix3r::Zero());

        BOOST_CHECK_CLOSE(u_dot_mu_plus_lambda_grad_div_u, ans, close_enough);
    }

    /* With nonzero refcoeffs, explicit operator differs from full one */
    {
        const Vector3r refcoeff_grad_div_m(3, 5, 7);
        const double ans
            = -218902.29059155252474980585723649734813459133691881587116865L;

        const double u_dot_mu_plus_lambda_grad_div_u
            = suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u(
                mu, lambda, rho, grad_rho, grad_grad_rho,
                m, div_m, grad_m, grad_div_m,
                refcoeff_grad_div_m, Matrix3r::Zero());

        BOOST_CHECK_CLOSE(u_dot_mu_plus_lambda_grad_div_u, ans, close_enough);
    }

    /* With zero refcoeffs, explicit operator matches the full one */
    {
        Matrix3r refcoeff_grad_grad_rho; // should be symmetric
        refcoeff_grad_grad_rho << 3,  5,  7,
                                  5, 11, 13,
                                  7, 13, 17;
        const double ans
            = -4059.2308357512977715114493722454403366670208895633912928559L;

        const double u_dot_mu_plus_lambda_grad_div_u
            = suzerain::rholt::explicit_u_dot_mu_plus_lambda_grad_div_u(
                mu, lambda, rho, grad_rho, grad_grad_rho,
                m, div_m, grad_m, grad_div_m,
                Vector3r::Zero(), refcoeff_grad_grad_rho);

        BOOST_CHECK_CLOSE(u_dot_mu_plus_lambda_grad_div_u, ans, close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double mu     = 4181.0;
        const double lambda = 1234.5;
        const double rho    = 5678.0;
        const Vector3r m(3.0, 5.0, 7.0);

        const Vector3r refcoeff_grad_div_m = suzerain::rholt
                ::explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m
                (mu, lambda, rho, m);

        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(refcoeff_grad_div_m[i],
                    (mu+lambda)/(rho*rho)*m[i], close_enough);
        }

        const Matrix3r refcoeff_grad_grad_rho = suzerain::rholt
                ::explicit_u_dot_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho<
                double, Vector3r, Matrix3r>(mu, lambda, rho, m);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                BOOST_CHECK_CLOSE(refcoeff_grad_grad_rho(i,j),
                        (mu+lambda)/(rho*rho*rho)*m[i]*m[j], close_enough);
            }
        }
    }
}

// Checks derived formula and computed result against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_curl_curl_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    /* Expected results found using curl_curl_u = grad_div_u - div_grad_u */
    {
        const Vector3r curl_curl_u
            = suzerain::rholt::curl_curl_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, curl_curl_m);

        const Vector3r grad_div_u
            = suzerain::rholt::grad_div_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);
        const Vector3r div_grad_u
            = suzerain::rholt::div_grad_u(
                    rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);
        const Vector3r ans = grad_div_u - div_grad_u;

        BOOST_CHECK_CLOSE(curl_curl_u(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(curl_curl_u(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(curl_curl_u(2), ans(2), close_enough);
    }

// TODO Implement mu_plus_lambda_curl_curl_u kernels (Ticket #2183)
//  /* With zero refcoeffs, explicit operator matches the full one */
//  {
//      const double mu     = 4181.0;
//      const double lambda = 6765.0;
//
//      const Vector3r mu_plus_lambda_grad_div_u
//          = suzerain::rholt::explicit_mu_plus_lambda_grad_div_u(
//                  mu, lambda, rho, grad_rho, grad_grad_rho, m,
//                  div_m, grad_m, grad_div_m, 0, Vector3r::Zero());
//
//      const Vector3r ans(
//          (mu+lambda)
//          *3.5808667611324763961641377901365615487146733886413982972779L,
//          (mu+lambda)
//          *-11.378277959392865631969701464497046416747061061234520455900L,
//          (mu+lambda)
//          *237.13623643835318300159939437852909295361695206836632310924L);
//
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
//  }
//
//  /* With nonzero refcoeffs, explicit operator differs from full one */
//  {
//      const double mu     = 4181.0;
//      const double lambda = 6765.0;
//      const double refcoeff_grad_div_m = 89.0;
//
//      const Vector3r mu_plus_lambda_grad_div_u
//          = suzerain::rholt::explicit_mu_plus_lambda_grad_div_u(
//                  mu, lambda, rho, grad_rho, grad_grad_rho, m,
//                  div_m, grad_m, grad_div_m,
//                  refcoeff_grad_div_m, Vector3r::Zero());
//
//      const Vector3r ans(
//              59251.997466444561194961843519913703789109645722525625271084L,
//              56371.864698208669179904529966597802712284277363839600322130L,
//            -225340.51082752741696350364123413038584340722233614860633701L);
//
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
//  }
//
//  /* With nonzero refcoeffs, explicit operator differs from full one */
//  {
//      const double mu     = 4181.0;
//      const double lambda = 6765.0;
//      const Vector3r refcoeff_grad_grad_rho(144.0, 233.0, 377.0);
//
//      const Vector3r mu_plus_lambda_grad_div_u
//          = suzerain::rholt::explicit_mu_plus_lambda_grad_div_u(
//                  mu, lambda, rho, grad_rho, grad_grad_rho, m,
//                  div_m, grad_m, grad_div_m,
//                  0, refcoeff_grad_grad_rho);
//
//      const Vector3r ans(
//           94481.167567356086632412652250834802712230814912068745762004L,
//          -90372.630543514307207540352230384670077713330376273060910286L,
//               2.6250052440542139411355069708673794514702911573403377727538e6L);
//
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(0), ans(0), close_enough);
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(1), ans(1), close_enough);
//      BOOST_CHECK_CLOSE(mu_plus_lambda_grad_div_u(2), ans(2), close_enough);
//  }
//
//  /* Ensure the coefficient calculations are correct */
//  {
//      const double mu     = 4181.0;
//      const double lambda = 6765.0;
//      const double rho    = 67.0;
//      const Vector3r m(144.0, 233.0, 377.0);
//
//      BOOST_CHECK_CLOSE(
//          suzerain::rholt
//              ::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_div_m(
//                  mu, lambda, rho),
//              (mu+lambda)/rho,
//              close_enough);
//      for (int i = 0; i < 3; ++i) {
//          BOOST_CHECK_CLOSE(
//              suzerain::rholt
//                  ::explicit_mu_plus_lambda_grad_div_u_refcoeff_grad_grad_rho(
//                      mu, lambda, rho, m)[i],
//                  (mu+lambda)/rho/rho*m[i],
//                  close_enough);
//      }
//  }
}

// Checks computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_tau_and_div_tau )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholt::p_T_mu_lambda(
            alpha, beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
    const double          div_u  = suzerain::rholt::div_u(rho, grad_rho, m, div_m);
    const Matrix3r grad_u = suzerain::rholt::grad_u(rho, grad_rho, m, grad_m);
    const Vector3r grad_div_u = suzerain::rholt::grad_div_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);
    const Vector3r div_grad_u = suzerain::rholt::div_grad_u(
                rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.02e4;
    {
        using namespace suzerain;
        const Matrix3r tau = rholt::tau(mu, lambda, div_u, grad_u);

        /* Expected results found using test_rholt.sage */
        Matrix3r ans;
        ans(0,0) = -432.77702868587701572931303205296522463376744126783255181324L;
        ans(0,1) = -18.274805458259731295309074758455374124262335886995620262487L;
        ans(0,2) = 91.817913251736077743933865033061417711733974379411774042718L;
        ans(1,0) = ans(0,1);
        ans(1,1) = -370.58480516905370299703749948866597794578072398454903444316L;
        ans(1,2) = 75.304386307037240507158364283048615894309818988016971749851L;
        ans(2,0) = ans(0,2);
        ans(2,1) = ans(1,2);
        ans(2,2) = -706.56388204197802262981359995686124154741489842153295213821L;

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
        const Vector3r div_tau = rholt::div_tau(
                mu, grad_mu, lambda, grad_lambda,
                div_u, grad_u, div_grad_u, grad_div_u);

        /* Expected results found using test_rholt.sage */
        Vector3r ans;
        ans(0) = -5.8572189782846240727678259765123413043664342815994539555325L;
        ans(1) = -500.16654308872810027481518191178470854090357788902597100115L;
        ans(2) = 24897.262970203388363153941644146159218248689630111325247307L;

        BOOST_CHECK_CLOSE(div_tau(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_tau(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_tau(2), ans(2), close_enough);
    }
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_e_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    using namespace suzerain;
    const Vector3r u = rholt::u(rho, m);
    const double div_u = rholt::div_u(
            rho, grad_rho, m, div_m);
    const double div_e_u = rholt::div_e_u(e, grad_e, u, div_u);

    /* Expected results found using test_rholt.sage */
    const double ans
        = -90849.212502207575911979472073826868082163956391101506206117L;

    const double close_enough = std::numeric_limits<double>::epsilon() * 1e2;
    BOOST_CHECK_CLOSE(div_e_u, ans, close_enough);
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_p_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon()*5.0e2;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    using namespace suzerain;
    rholt::p_T_mu_lambda(
            alpha, beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* Expected results found using test_rholt.sage */
    {
        const Vector3r u = rholt::u(rho, m);
        const double div_u = rholt::div_u(
                rho, grad_rho, m, div_m);
        const double div_p_u = rholt::div_p_u(p, grad_p, u, div_u);
        const double ans
            = -36133.705926299567425548649114294679973238338337791365684538L;

        BOOST_CHECK_CLOSE(div_p_u, ans, close_enough);
    }

    /* Explicit operator should coincide when refcoeffs are zero */
    {
        const double explicit_div_p_u = rholt::explicit_div_p_u(
                rho, grad_rho, m, div_m, p, grad_p,
                0, Vector3r::Zero());
        const double ans
            = -36133.705926299567425548649114294679973238338337791365684538L;
        BOOST_CHECK_CLOSE(explicit_div_p_u, ans, close_enough);
    }

    /* Explicit operator differs for nonzero refcoeff */
    {
        const double refcoeff_div_m = 77.0;
        const double explicit_div_p_u = rholt::explicit_div_p_u(
                rho, grad_rho, m, div_m, p, grad_p,
                refcoeff_div_m, Vector3r::Zero());
        const double ans
            = 49409.670159954180179078803321626238490759844350764069895944L;
        BOOST_CHECK_CLOSE(explicit_div_p_u, ans, close_enough);
    }

    /* Explicit operator differs for nonzero refcoeff */
    {
        const Vector3r refcoeff_grad_rho(77.0, 88.0, 99.0);
        const double explicit_div_p_u = rholt::explicit_div_p_u(
                rho, grad_rho, m, div_m, p, grad_p,
                0, refcoeff_grad_rho);
        const double ans
            = -9403.7059262995674255486491142946799732383383377913656845385L;
        BOOST_CHECK_CLOSE(explicit_div_p_u, ans, close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double rho    = 67.0;
        const Vector3r m(144.0, 233.0, 377.0);
        const double p      = 55.0;

        BOOST_CHECK_CLOSE(rholt::explicit_div_p_u_refcoeff_div_m(
            rho, p), p/rho, close_enough);
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_CLOSE(rholt::explicit_div_p_u_refcoeff_grad_rho(
                rho, m, p)[i], p/rho/rho*m[i], close_enough);
        }
    }
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_grad_p )
{
    const double close_enough = std::numeric_limits<double>::epsilon()*5.0e2;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double gamma = 1.4;

    using namespace suzerain;

    /* Should recover the full operator when no refcoeffs are used */
    /* and the result is appropriately augmented by a grad_e term. */
    {
        const Vector3r ans = (gamma - 1)*grad_e
            + rholt::explicit_grad_p(gamma, rho, grad_rho, m, grad_m,
                                     0, Vector3r::Zero());

        const Vector3r grad_p_ans(
                5453.5209639604035869584122627487425816813409506343630929293L,
                3086.9691968926047335425739267829758568639555587282859345292L,
                1999.8806276135462562543621611330864774644251338499615213427L);
        BOOST_CHECK_CLOSE(ans(0), grad_p_ans(0), close_enough);
        BOOST_CHECK_CLOSE(ans(1), grad_p_ans(1), close_enough);
        BOOST_CHECK_CLOSE(ans(2), grad_p_ans(2), close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    {
        const double refcoeff_grad_rho = 7;
        const Vector3r ans = (gamma - 1)*grad_e
            + (gamma - 1)/2 * refcoeff_grad_rho * grad_rho
            + rholt::explicit_grad_p(gamma, rho, grad_rho, m, grad_m,
                                     refcoeff_grad_rho,
                                     Vector3r::Zero());

        const Vector3r grad_p_ans(
                5453.5209639604035869584122627487425816813409506343630929293L,
                3086.9691968926047335425739267829758568639555587282859345292L,
                1999.8806276135462562543621611330864774644251338499615213427L);
        BOOST_CHECK_CLOSE(ans(0), grad_p_ans(0), close_enough);
        BOOST_CHECK_CLOSE(ans(1), grad_p_ans(1), close_enough);
        BOOST_CHECK_CLOSE(ans(2), grad_p_ans(2), close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    {
        const Vector3r refcoeff_grad_m(2, 3, 5);
        const Vector3r ans = (gamma - 1)*grad_e
            - (gamma - 1)*grad_m.transpose()*refcoeff_grad_m
            + rholt::explicit_grad_p(gamma, rho, grad_rho, m, grad_m,
                                     0, refcoeff_grad_m);

        const Vector3r grad_p_ans(
                5453.5209639604035869584122627487425816813409506343630929293L,
                3086.9691968926047335425739267829758568639555587282859345292L,
                1999.8806276135462562543621611330864774644251338499615213427L);
        BOOST_CHECK_CLOSE(ans(0), grad_p_ans(0), close_enough);
        BOOST_CHECK_CLOSE(ans(1), grad_p_ans(1), close_enough);
        BOOST_CHECK_CLOSE(ans(2), grad_p_ans(2), close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    {
        const double rho    = 67.0;
        const Vector3r m(144.0, 233.0, 377.0);

        BOOST_CHECK_CLOSE(rholt::explicit_grad_p_refcoeff_grad_rho(rho, m),
            217154.0L/4489.0L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_grad_p_refcoeff_grad_m(rho, m)[0],
            144.0L/67.0L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_grad_p_refcoeff_grad_m(rho, m)[1],
            233.0L/67.0L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_grad_p_refcoeff_grad_m(rho, m)[2],
            377.0L/67.0L,
            close_enough);
    }
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_e_plus_p_u )
{
    const double close_enough = std::numeric_limits<double>::epsilon()*5.0e2;

    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    using namespace suzerain;
    rholt::p_T_mu_lambda(
            alpha, beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* Should recover the full operator when no refcoeffs used */
    {
        const double ans
            = -126982.91842850714333752812118812154805540229472889287189066L;

        const double div_e_plus_p_u = rholt::explicit_div_e_plus_p_u(
                gamma, rho, grad_rho, m, div_m, grad_m, e, grad_e, p,
                0, Vector3r::Zero(), Vector3r::Zero());
        BOOST_CHECK_CLOSE(div_e_plus_p_u, ans, close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    {
        const double refcoeff_div_m = 55.0;
        const double div_e_plus_p_u = rholt::explicit_div_e_plus_p_u(
                gamma, rho, grad_rho, m, div_m, grad_m, e, grad_e, p,
                refcoeff_div_m, Vector3r::Zero(), Vector3r::Zero());
        const double ans
            = -65880.506938325895048508512305320892009689307094210417904597L;
        BOOST_CHECK_CLOSE(div_e_plus_p_u, ans, close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    {
        const Vector3r refcoeff_grad_rho(55.0, 66.0, 77.0);
        const double div_e_plus_p_u = rholt::explicit_div_e_plus_p_u(
                gamma, rho, grad_rho, m, div_m, grad_m, e, grad_e, p,
                0, refcoeff_grad_rho, Vector3r::Zero());
        const double ans
            = -146826.91842850714333752812118812154805540229472889287189066L;
        BOOST_CHECK_CLOSE(div_e_plus_p_u, ans, close_enough);
    }

    /* Explicit operator differs when refcoeffs in use */
    /* Odd looking 1/gamma factors due to removing gamma from refcoeff */
    {
        const Vector3r refcoeff_grad_e(55/gamma, 66/gamma, 77/gamma);
        const double div_e_plus_p_u = rholt::explicit_div_e_plus_p_u(
                gamma, rho, grad_rho, m, div_m, grad_m, e, grad_e, p,
                0, Vector3r::Zero(), refcoeff_grad_e);
        const double ans
            = -1.7876309184285071433375281211881215480554022947288928718907e6L;
        BOOST_CHECK_CLOSE(div_e_plus_p_u, ans, close_enough);
    }

    /* Ensure the coefficient calculations are correct */
    /* Odd looking 1/gamma factors due to removing gamma from refcoeff */
    {
        const double rho    = 67.0;
        const Vector3r m(144.0, 233.0, 377.0);
        const double e      = 55.0;
        const double p      = (gamma-1)*(e - m.squaredNorm()/(2*rho));

        BOOST_CHECK_CLOSE(rholt::explicit_div_e_plus_p_u_refcoeff_div_m(
            rho, e, p),
            -8.5256850077968367119625751837825796391178436177322343506349L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_div_e_plus_p_u_refcoeff_grad_rho(
            gamma, rho, m, e, p)[0],
            39.117758500879429983076375751006606530723526497607750953408L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_div_e_plus_p_u_refcoeff_grad_rho(
            gamma, rho, m, e, p)[1],
            63.294706463228522125394413541559300844851261624601430362112L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_div_e_plus_p_u_refcoeff_grad_rho(
            gamma, rho, m, e, p)[2],
            102.41246496410795210847078929256590737557478812220918131552L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_div_e_plus_p_u_refcoeff_grad_e(
            rho, m)[0], (1 / gamma) *
            3.0089552238805970149253731343283582089552238805970149253731L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_div_e_plus_p_u_refcoeff_grad_e(
            rho, m)[1], (1 / gamma) *
            4.8686567164179104477611940298507462686567164179104477611940L,
            close_enough);
        BOOST_CHECK_CLOSE(rholt::explicit_div_e_plus_p_u_refcoeff_grad_e(
            rho, m)[2], (1 / gamma) *
            7.8776119402985074626865671641791044776119402985074626865672L,
            close_enough);
    }
}


// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_tau_u )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    using namespace suzerain;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholt::p_T_mu_lambda(
            alpha, beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    const Vector3r u = rholt::u(rho, m);
    const double div_u = rholt::div_u(
            rho, grad_rho, m, div_m);
    const Matrix3r grad_u = rholt::grad_u(
            rho, grad_rho, m, grad_m);
    const Vector3r div_grad_u = rholt::div_grad_u(
            rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);
    const Vector3r grad_div_u = rholt::grad_div_u(
            rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);
    const Matrix3r tau = rholt::tau(mu, lambda, div_u, grad_u);
    const Vector3r div_tau = rholt::div_tau(
            mu, grad_mu, lambda, grad_lambda,
            div_u, grad_u, div_grad_u, grad_div_u);

    const double div_tau_u = rholt::div_tau_u<double>(
            u, grad_u, tau, div_tau);

    /* Expected results found using test_rholt.sage */
    const double ans
        = -4619.0441415000015916988672295614662059743675820473242574398L;

    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e3;
    BOOST_CHECK_CLOSE(div_tau_u, ans, close_enough);
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_mu_grad_T )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    const double alpha = 5.0;
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Vector3r grad_p, grad_T, grad_mu, grad_lambda;

    suzerain::rholt::p_T_mu_lambda(
            alpha, beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    const double div_grad_p = suzerain::rholt::div_grad_p(
        gamma,
        rho, grad_rho, div_grad_rho,
        m, grad_m, div_grad_m,
        e, grad_e, div_grad_e);
    const double div_grad_T = suzerain::rholt::div_grad_T(
        gamma,
        rho, grad_rho, div_grad_rho,
        p, grad_p, div_grad_p);

    const double div_mu_grad_T = suzerain::rholt::div_mu_grad_T(
            grad_T, div_grad_T, mu, grad_mu);

    /* Expected results found using test_rholt.sage */
    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e3;
    BOOST_CHECK_CLOSE(div_mu_grad_T,
            1082.5428069809810057587904051307120521551327361396488157143L,
            close_enough);
}

// Checks derived formula and computation against rholt_test_data()
BOOST_AUTO_TEST_CASE( rholt_div_rho_inverse_m_outer_m )
{
    ADD_TEST_DATA_TO_ENCLOSING_SCOPE;

    /* Expected results found using test_rholt.sage */
    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e3;
    const Vector3r ans(
            - 139.62394416218589511382305374934119093270491143228219517924L,
            - 196.62019324654648600809920433959614291380888927326388185922L,
             1352.1114315042037165058005031972730396574030986356606221464L );

    /* Compute using div_rho_inverse_m_outer_m */
    {
        const Vector3r div_rho_inverse_m_outer_m
            = suzerain::rholt::div_rho_inverse_m_outer_m(
                    rho, grad_rho, m, div_m, grad_m);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(2), ans(2), close_enough);
    }

    /* Compute using div_u_outer_m */
    {
        const Vector3r u = suzerain::rholt::u(
                rho, m);
        const double div_u = suzerain::rholt::div_u(
                rho, grad_rho, m, div_m);

        const Vector3r div_u_outer_m
            = suzerain::rholt::div_u_outer_m(
                    m, grad_m, u, div_u);

        BOOST_CHECK_CLOSE(div_u_outer_m(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_u_outer_m(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_u_outer_m(2), ans(2), close_enough);
    }

    /* Compute using explicit_div_rho_inverse_m_outer_m without refs */
    {
        const Vector3r u = suzerain::rholt::u(
                rho, m);
        const Vector3r div_rho_inverse_m_outer_m
            = suzerain::rholt::explicit_div_rho_inverse_m_outer_m(
                    grad_rho, div_m, grad_m, u,
                    Vector3r::Zero(), Matrix3r::Zero());
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(2), ans(2), close_enough);
    }

    /* Compute using explicit_div_rho_inverse_m_outer_m with vector ref */
    {
        const Vector3r u = suzerain::rholt::u(
                rho, m);
        const Vector3r uref(7, 11, 13);
        const Vector3r div_rho_inverse_m_outer_m
            = suzerain::rholt::explicit_div_rho_inverse_m_outer_m(
                    grad_rho, div_m, grad_m, u,
                    uref, Matrix3r::Zero());
        const Vector3r refadjust = grad_m*uref + div_m*uref;
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(0),
                          ans(0) - refadjust(0), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(1),
                          ans(1) - refadjust(1), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(2),
                          ans(2) - refadjust(2), close_enough);
    }

    /* Compute using explicit_div_rho_inverse_m_outer_m with tensor ref */
    {
        const Vector3r u = suzerain::rholt::u(
                rho, m);
        const Vector3r uref(7, 11, 13);
        const Vector3r div_rho_inverse_m_outer_m
            = suzerain::rholt::explicit_div_rho_inverse_m_outer_m(
                    grad_rho, div_m, grad_m, u,
                    Vector3r::Zero(), uref*uref.transpose());
        const Vector3r refadjust = -uref*uref.transpose()*grad_rho;
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(0),
                          ans(0) - refadjust(0), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(1),
                          ans(1) - refadjust(1), close_enough);
        BOOST_CHECK_CLOSE(div_rho_inverse_m_outer_m(2),
                          ans(2) - refadjust(2), close_enough);
    }
}
