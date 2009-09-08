#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <cmath>
#include <limits>
#include <Eigen/Core>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <suzerain/classical.hpp>

// Provides data for the test field
//      rho = 2*(x^2)*y*z + 3*x*(y^2)*z + 5*x*y*(z^2)
//      m   = {
//               7*sin( 7*x) + 11*sin(11*y) + 13*sin(13*z),
//              17*sin(17*x) + 19*sin(19*y) + 23*sin(23*z),
//              31*sin(31*x) + 37*sin(37*y) + 41*sin(41*z)
//      }
//      e   = 311*(x^2)*y*z + 313*x*(y^2)*z + 317*x*y*(z^2)
// symbolically evaluated at (x,y,z) = (1, 2, 3)
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
        Eigen::Vector3d &grad_e)
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

    m(0) =  7.*sin( 7.) + 11.*sin(22.) + 13.*sin(39.);
    m(1) = 17.*sin(17.) + 19.*sin(38.) + 23.*sin(69.);
    m(2) = 31.*sin(31.) + 37.*sin(74.) + 41.*sin(123.);

    div_m = 49.*cos(7.) + 361.*cos(38.) + 1681.*cos(123.);

    grad_m(0,0) =   49.*cos(  7.);
    grad_m(0,1) =  121.*cos( 22.);
    grad_m(0,2) =  169.*cos( 39.);
    grad_m(1,0) =  289.*cos( 17.);
    grad_m(1,1) =  361.*cos( 38.);
    grad_m(1,2) =  529.*cos( 69.);
    grad_m(2,0) =  961.*cos( 31.);
    grad_m(2,1) = 1369.*cos( 74.);
    grad_m(2,2) = 1681.*cos(123.);

    div_grad_m(0) = -  343.*sin( 7.) -  1331.*sin(22.) -  2197.*sin( 39.);
    div_grad_m(1) = - 4913.*sin(17.) -  6859.*sin(38.) - 12167.*sin( 69.);
    div_grad_m(2) = -29791.*sin(31.) - 50653.*sin(74.) - 68921.*sin(123.);

    grad_div_m(0) = -  343.*sin(  7.);
    grad_div_m(1) = - 6859.*sin( 38.);
    grad_div_m(2) = -68921.*sin(123.);

    e = 11328.;

    grad_e(0) = 13194.;
    grad_e(1) =  7542.;
    grad_e(2) =  5678.;
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

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e);

    const double beta  = 2.0/3.0;
    const double gamma = 1.4;

    double p, T, mu, lambda;
    Eigen::Vector3d grad_p, grad_T, grad_mu, grad_lambda;

    pecos::suzerain::orthonormal::rhome::p_T_mu_lambda(
            beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

    /* Expected results found using sage's RealField(200) */
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

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e);

    const Eigen::Matrix3d grad_u
        = pecos::suzerain::orthonormal::rhome::grad_u(
                rho, grad_rho, m, grad_m);

    Eigen::Matrix3d ans; /* Found using sage's RealField(200) */
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

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e);

    const double div_u = pecos::suzerain::orthonormal::rhome::div_u(
            rho, grad_rho, m, div_m);

    const double ans = -7.8528271460331434259092626891655172292938028930092832721830;
    const double close_enough = std::numeric_limits<double>::epsilon();
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

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e);

    const Eigen::Vector3d grad_div_u
        = pecos::suzerain::orthonormal::rhome::grad_div_u(
                rho, grad_rho, grad_grad_rho, m, div_m, grad_m, grad_div_m);

    Eigen::Vector3d ans; /* Found using sage's RealField(200) */
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

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e);

    const Eigen::Vector3d div_grad_u
        = pecos::suzerain::orthonormal::rhome::div_grad_u(
                rho, grad_rho, div_grad_rho, m, grad_m, div_grad_m);

    Eigen::Vector3d ans; /* Found using sage's RealField(200) */
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

    orthonormal_rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e);

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

        Eigen::Matrix3d ans; /* Found using sage's RealField(200) */
        ans(0,0) =   70.531543279759231389408345113198923408553579956805627651633;
        ans(0,1) = - 18.274805458259731295309074758455374124262335886995620262487;
        ans(0,2) =   91.817913251736077743933865033061417711733974379411774042718;
        ans(1,0) = ans(0,1);
        ans(1,1) =  132.72376679658254412168387767749817009654029724008914502171;
        ans(1,2) =   75.304386307037240507158364283048615894309818988016971749851;
        ans(2,0) = ans(0,2);
        ans(2,1) = ans(1,2);
        ans(2,2) = -203.25531007634177551109222279069709350509387719689477267335;

        BOOST_CHECK_CLOSE(tau(0,0), ans(0,0), close_enough);
        BOOST_CHECK_CLOSE(tau(0,1), ans(0,1), close_enough);
        BOOST_CHECK_CLOSE(tau(0,2), ans(0,2), close_enough);
        BOOST_CHECK_CLOSE(tau(1,0), ans(1,0), close_enough);
        BOOST_CHECK_CLOSE(tau(1,1), ans(1,1), close_enough);
        BOOST_CHECK_CLOSE(tau(1,2), ans(1,2), close_enough);
        BOOST_CHECK_CLOSE(tau(2,0), ans(2,0), close_enough);
        BOOST_CHECK_CLOSE(tau(2,1), ans(2,1), close_enough);
        BOOST_CHECK_CLOSE(tau(2,2), ans(2,2), close_enough);
    }

    {
        using namespace pecos::suzerain;
        const Eigen::Vector3d div_tau = orthonormal::div_tau(
                mu, grad_mu, lambda, grad_lambda,
                div_u, grad_u, div_grad_u, grad_div_u);

        Eigen::Vector3d ans; /* Found using sage's RealField(200) */
        ans(0) = - 195.58731950891911935327871659569254694467509625810898354185;
        ans(1) =   246.52589132115431226459081257854313932950474889459145467105;
        ans(2) =  9662.1147275452378450481114526062453379951497410467223391699;

        BOOST_CHECK_CLOSE(div_tau(0), ans(0), close_enough);
        BOOST_CHECK_CLOSE(div_tau(1), ans(1), close_enough);
        BOOST_CHECK_CLOSE(div_tau(2), ans(2), close_enough);
    }
}
