#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <cmath>
#include <limits>
#include <Eigen/Core>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <suzerain/orthonormal.hpp>

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
// Floating point values evaluated using test_orthonormal.sage
// to 200 bits of precision.
void orthonormal_rhome_test_data(
        double          &rho,
        Eigen::Vector3d &grad_rho,
        double          &div_grad_rho,
        Eigen::Matrix3d &grad_grad_rho,
        Eigen::Vector3d &m,
        double          &div_m,
        Eigen::Matrix3d &grad_m,
        Eigen::Vector3d &div_grad_m,
        Eigen::Vector3d &grad_div_m,
        double          &e,
        Eigen::Vector3d &grad_e,
        double          &div_grad_e)
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

    m(0) =  17.030881810530221790029045449138418084725998922603314999607;
    m(1) = -13.352805083487451615530843887500866367157365945228654781328;
    m(2) = -67.831621760613409212021683467570521921918364223433827929255;

    div_m = -1110.9529361851136052549019796872846553765997751760446179283;

    grad_m(0,0) =    36.941210462821927268918678564239918598697356181546215137637;
    grad_m(0,1) =  -120.99525999375109230095514344473594667045460456500973446019;
    grad_m(0,2) =    45.062655568829395508033559350881935096816196590558655110591;
    grad_m(1,0) =  - 79.522204696911510521679027743488377694525395918294315600644;
    grad_m(1,1) =   344.78158550107344358873958098152158053944887665250279502393;
    grad_m(1,2) =   525.50351087308169627005787890310679646720537555710570892294;
    grad_m(2,0) =   879.06740585015455908290082893049093028384492918816862062545;
    grad_m(2,1) =   235.08104096633447429877919230313280932462835892676775285483;
    grad_m(2,2) = -1492.6757321490089761125602392330461545147460080100936280899;

    div_grad_m(0) = - 2331.0237743611578930680818995644668983923397087883721673025;
    div_grad_m(1) =   4087.1406255366278601045485877398647158902558157286863079530;
    div_grad_m(2) =  93634.307505134882301285667998181793756271588535303593390142;

    grad_div_m(0) = -  225.34640336054465800617068841661686603234641360063909560764;
    grad_div_m(1) = - 2032.7920813676738919937627213144098066291866038214905756451;
    grad_div_m(2) =  31697.008481817318630325961933724829632738184041308835720121;

    e = 11328.;

    grad_e(0) = 13194.;
    grad_e(1) =  7542.;
    grad_e(2) =  5678.;

    div_grad_e = 6878.;
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_rhome_u )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const Eigen::Vector3d u = pecos::suzerain::orthonormal::rhome::u(rho, m);

    /* Expected results found using test_orthonormal.sage */
    const Eigen::Vector3d ans(
             0.12341218703282769413064525687781462380236231103335735506961,
            -0.096759457126720663880658289039861350486647579313251121603828,
            -0.49153349101893774791320060483746755015882872625676686905257);

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e2;
    BOOST_CHECK_CLOSE(u(0), ans(0), close_enough);
    BOOST_CHECK_CLOSE(u(1), ans(1), close_enough);
    BOOST_CHECK_CLOSE(u(2), ans(2), close_enough);
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_rhome_p_T_mu_lambda )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Eigen::Vector3d grad_p, grad_T, grad_mu, grad_lambda;

    pecos::suzerain::orthonormal::rhome::p_T_mu_lambda(
            beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* Expected results found using test_orthonormal.sage */
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e4;
    BOOST_CHECK_CLOSE(p,
            4523.8529315224394491652244924378285975908038303556842666509,
            close_enough);
    BOOST_CHECK_CLOSE(T,
            45.894160174865327745154451372557681424834241757231579516748,
            close_enough);
    BOOST_CHECK_CLOSE(mu,
            12.818531787494714755711740676149526369900526768964830775976,
            close_enough);
    BOOST_CHECK_CLOSE(lambda,
            -8.5456878583298098371411604507663509132670178459765538506509,
            close_enough);

    const Eigen::Vector3d grad_p_ans(
            5453.5209639604035869584122627487425816813409506343630929293,
            3086.9691968926047335425739267829758568639555587282859345292,
            1999.8806276135462562543621611330864774644251338499615213427);
    BOOST_CHECK_CLOSE(grad_p(0), grad_p_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_p(1), grad_p_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_p(2), grad_p_ans(2), close_enough);

    const Eigen::Vector3d grad_T_ans(
             5.4406182848896076809319526229317927581792830964012420477440,
             2.3838039162055298052983060006061443162968025314525571766942,
            -4.9864006857304358686639947733917588394000520736206805318339);
    BOOST_CHECK_CLOSE(grad_T(0), grad_T_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_T(1), grad_T_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_T(2), grad_T_ans(2), close_enough);

    const Eigen::Vector3d grad_mu_ans(
             1.0130662690381109382304916850142818302534744683039668368759,
             0.44387442989997855519790991486035061354804176216270679061604,
            -0.92848901983288259118625731085663474679900828620195051533844);
    BOOST_CHECK_CLOSE(grad_mu(0), grad_mu_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_mu(1), grad_mu_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_mu(2), grad_mu_ans(2), close_enough);

    const Eigen::Vector3d grad_lambda_ans(
            -0.67537751269207395882032779000952122016898297886931122458391,
            -0.29591628659998570346527327657356707569869450810847119374402,
             0.61899267988858839412417154057108983119933885746796701022562);
    BOOST_CHECK_CLOSE(grad_lambda(0), grad_lambda_ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_lambda(1), grad_lambda_ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_lambda(2), grad_lambda_ans(2), close_enough);
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_rhome_div_grad_p_and_div_grad_T )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Eigen::Vector3d grad_p, grad_T, grad_mu, grad_lambda;

    pecos::suzerain::orthonormal::rhome::p_T_mu_lambda(
            beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    const double div_grad_p = pecos::suzerain::orthonormal::rhome::div_grad_p(
        gamma,
        rho, grad_rho, div_grad_rho,
        m, grad_m, div_grad_m,
        e, grad_e, div_grad_e);

    /* Expected results found using test_orthonormal.sage */
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;
    BOOST_CHECK_CLOSE(div_grad_p,
            11191.566068848019028203711722987453542724808010087855921969,
            close_enough);

    const double div_grad_T = pecos::suzerain::orthonormal::rhome::div_grad_T(
        gamma,
        rho, grad_rho, div_grad_rho,
        p, grad_p, div_grad_p);

    /* Expected results found using test_orthonormal.sage */
    BOOST_CHECK_CLOSE(div_grad_T,
            83.577681904999696238558381171408689250539040369361349554271,
            close_enough);
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_rhome_grad_u )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const Eigen::Matrix3d grad_u
        = pecos::suzerain::orthonormal::rhome::grad_u(
                rho, grad_rho, m, grad_m);

    /* Expected results found using test_orthonormal.sage */
    Eigen::Matrix3d ans;
    ans(0,0) =   0.13354624933259255905305717414904148571263050381552617302315;
    ans(0,1) = - 0.95458058163483407021971942603699868798014583786167988660320;
    ans(0,2) =   0.25857485039372819387032260745049292527417870255089489945870;
    ans(1,0) = - 0.47107453715872036912739336512687808059078448566164237217442;
    ans(1,1) =   2.5594178135586821836692525516520979567520812757446061058222;
    ans(1,2) =   3.8612842725703801936593326729720018775665984897457448852524;
    ans(2,0) =   6.9043291992970668207962385482326888609251394067150989201691;
    ans(2,1) =   2.0133656136592902780233887313332788854235250587761338439305;
    ans(2,2) = -10.545791208924418168631572414966656671758514672569415551028;

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

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_rhome_div_u )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const double div_u = pecos::suzerain::orthonormal::rhome::div_u(
            rho, grad_rho, m, div_m);

    const double ans = -7.8528271460331434259092626891655172292938028930092832721830;
    const double close_enough = std::numeric_limits<double>::epsilon() *1e2;
    BOOST_CHECK_CLOSE(div_u, ans, close_enough);
}

// Checks derived formula and computed result against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_rhome_grad_div_u )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const Eigen::Vector3d grad_div_u
        = pecos::suzerain::orthonormal::rhome::grad_div_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);

    /* Expected results found using test_orthonormal.sage */
    Eigen::Vector3d ans;
    ans(0) =   3.5808667611324763961641377901365615487146733886413982972779;
    ans(1) = -11.378277959392865631969701464497046416747061061234520455900;
    ans(2) = 237.13623643835318300159939437852909295361695206836632310924;

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;
    BOOST_CHECK_CLOSE(grad_div_u(0), ans(0), close_enough);
    BOOST_CHECK_CLOSE(grad_div_u(1), ans(1), close_enough);
    BOOST_CHECK_CLOSE(grad_div_u(2), ans(2), close_enough);
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_rhome_div_grad_u )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const Eigen::Vector3d div_grad_u
        = pecos::suzerain::orthonormal::rhome::div_grad_u(
                rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);

    /* Expected results found using test_orthonormal.sage */
    Eigen::Vector3d ans;
    ans(0) = - 16.318446092843163297609832709693050751558008045050013765743;
    ans(1) =   23.204406985769508424700716673327465318352740571910767421546;
    ans(2) =  672.79795991861399600487930607996269111701083307478613865859;

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;
    BOOST_CHECK_CLOSE(div_grad_u(0), ans(0), close_enough);
    BOOST_CHECK_CLOSE(div_grad_u(1), ans(1), close_enough);
    BOOST_CHECK_CLOSE(div_grad_u(2), ans(2), close_enough);
}

// Checks computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_tau_and_div_tau )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Eigen::Vector3d grad_p, grad_T, grad_mu, grad_lambda;

    using namespace pecos::suzerain::orthonormal;
    rhome::p_T_mu_lambda(beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
    const double          div_u  = rhome::div_u(rho, grad_rho, m, div_m);
    const Eigen::Matrix3d grad_u = rhome::grad_u(rho, grad_rho, m, grad_m);
    const Eigen::Vector3d grad_div_u = rhome::grad_div_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);
    const Eigen::Vector3d div_grad_u = rhome::div_grad_u(
                rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;
    {
        using namespace pecos::suzerain;
        const Eigen::Matrix3d tau = orthonormal::tau(mu, lambda, div_u, grad_u);

        /* Expected results found using test_orthonormal.sage */
        Eigen::Matrix3d ans;
        ans(0,0) =   70.531543279759231389408345113198923408553579956805627651633;
        ans(0,1) = - 18.274805458259731295309074758455374124262335886995620262487;
        ans(0,2) =   91.817913251736077743933865033061417711733974379411774042718;
        ans(1,0) = ans(0,1);
        ans(1,1) =  132.72376679658254412168387767749817009654029724008914502171;
        ans(1,2) =   75.304386307037240507158364283048615894309818988016971749851;
        ans(2,0) = ans(0,2);
        ans(2,1) = ans(1,2);
        ans(2,2) = -203.25531007634177551109222279069709350509387719689477267335;

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
        using namespace pecos::suzerain;
        const Eigen::Vector3d div_tau = orthonormal::div_tau(
                mu, grad_mu, lambda, grad_lambda,
                div_u, grad_u, div_grad_u, grad_div_u);

        /* Expected results found using test_orthonormal.sage */
        Eigen::Vector3d ans;
        ans(0) = - 195.58731950891911935327871659569254694467509625810898354185;
        ans(1) =   246.52589132115431226459081257854313932950474889459145467105;
        ans(2) =  9662.1147275452378450481114526062453379951497410467223391699;

        BOOST_CHECK_CLOSE(div_tau(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_tau(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_tau(2), ans(2), close_enough);
    }
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_div_e_u )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    using namespace pecos::suzerain;
    const Eigen::Vector3d u = orthonormal::rhome::u(rho, m);
    const double div_u = orthonormal::rhome::div_u(
            rho, grad_rho, m, div_m);
    const double div_e_u = orthonormal::div_e_u(e, grad_e, u, div_u);

    /* Expected results found using test_orthonormal.sage */
    const double ans
        = -90849.212502207575911979472073826868082163956391101506206117;

    const double close_enough = std::numeric_limits<double>::epsilon() * 1e2;
    BOOST_CHECK_CLOSE(div_e_u, ans, close_enough);
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_div_p_u )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Eigen::Vector3d grad_p, grad_T, grad_mu, grad_lambda;

    using namespace pecos::suzerain;
    orthonormal::rhome::p_T_mu_lambda(
            beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    const Eigen::Vector3d u = orthonormal::rhome::u(rho, m);
    const double div_u = orthonormal::rhome::div_u(
            rho, grad_rho, m, div_m);
    const double div_p_u = orthonormal::div_p_u(p, grad_p, u, div_u);

    /* Expected results found using test_orthonormal.sage */
    const double ans
        = -36133.705926299567425548649114294679973238338337791365684538;

    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e2;
    BOOST_CHECK_CLOSE(div_p_u, ans, close_enough);
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_div_tau_u )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    using namespace pecos::suzerain;

    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Eigen::Vector3d grad_p, grad_T, grad_mu, grad_lambda;

    pecos::suzerain::orthonormal::rhome::p_T_mu_lambda(
            beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    const Eigen::Vector3d u = orthonormal::rhome::u(rho, m);
    const double div_u = orthonormal::rhome::div_u(
            rho, grad_rho, m, div_m);
    const Eigen::Matrix3d grad_u = orthonormal::rhome::grad_u(
            rho, grad_rho, m, grad_m);
    const Eigen::Vector3d div_grad_u = orthonormal::rhome::div_grad_u(
            rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);
    const Eigen::Vector3d grad_div_u = orthonormal::rhome::grad_div_u(
            rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);
    const Eigen::Matrix3d tau = orthonormal::tau(mu, lambda, div_u, grad_u);
    const Eigen::Vector3d div_tau = orthonormal::div_tau(
            mu, grad_mu, lambda, grad_lambda,
            div_u, grad_u, div_grad_u, grad_div_u);

    const double div_tau_u = orthonormal::div_tau_u<double>(
            u, grad_u, tau, div_tau);

    /* Expected results found using test_orthonormal.sage */
    const double ans
        = -1178.5183176047041859331013266287927977178206769790673702902;

    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e3;
    BOOST_CHECK_CLOSE(div_tau_u, ans, close_enough);
}

// Checks derived formula and computation against orthonormal_rhome_test_data()
BOOST_AUTO_TEST_CASE( orthonormal_rhome_div_mu_grad_T )
{
    double          rho;
    Eigen::Vector3d grad_rho;
    double          div_grad_rho;
    Eigen::Matrix3d grad_grad_rho;
    Eigen::Vector3d m;
    double          div_m;
    Eigen::Matrix3d grad_m;
    Eigen::Vector3d div_grad_m;
    Eigen::Vector3d grad_div_m;
    double          e;
    Eigen::Vector3d grad_e;
    double          div_grad_e;

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e, div_grad_e);

    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Eigen::Vector3d grad_p, grad_T, grad_mu, grad_lambda;

    pecos::suzerain::orthonormal::rhome::p_T_mu_lambda(
            beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    const double div_grad_p = pecos::suzerain::orthonormal::rhome::div_grad_p(
        gamma,
        rho, grad_rho, div_grad_rho,
        m, grad_m, div_grad_m,
        e, grad_e, div_grad_e);
    const double div_grad_T = pecos::suzerain::orthonormal::rhome::div_grad_T(
        gamma,
        rho, grad_rho, div_grad_rho,
        p, grad_p, div_grad_p);

    const double div_mu_grad_T = pecos::suzerain::orthonormal::div_mu_grad_T(
            grad_T, div_grad_T, mu, grad_mu);

    /* Expected results found using test_orthonormal.sage */
    const double close_enough = std::numeric_limits<double>::epsilon()*1.0e3;
    BOOST_CHECK_CLOSE(div_mu_grad_T,
            1082.5428069809810057587904051307120521551327361396488157143,
            close_enough);
}

