#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <cmath>
#include <limits>
#include <Eigen/Core>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <suzerain/classical.hpp>

BOOST_AUTO_TEST_CASE( rhome_p_T_mu_lambda )
{
    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double rho   = 2.0;
    Eigen::Vector3d m;
    m << 3.0, 5.0, 7.0;
    const double e = 163.0;
    double p, T, mu, lambda;

    pecos::suzerain::orthonormal::rhome::p_T_mu_lambda(
            beta, gamma, rho, m, e, p, T, mu, lambda);

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e4;
    BOOST_CHECK_CLOSE(p,       56.9,             close_enough);
    BOOST_CHECK_CLOSE(T,       39.83,            close_enough);
    BOOST_CHECK_CLOSE(mu,      11.6629085673383, close_enough);
    BOOST_CHECK_CLOSE(lambda,  -7.7752723782255, close_enough);
}

// Provides data for the test field
//      rho = 2*(x^2)*y*z + 3*x*(y^2)*z + 5*x*y*(z^2)
//      m   = {
//               7*sin( 7*x) + 11*sin(11*y) + 13*sin(13*z),
//              17*sin(17*x) + 19*sin(19*y) + 23*sin(23*z),
//              31*sin(31*x) + 37*sin(37*y) + 41*sin(41*z)
//      }
//      e   = 311*(x^2)*y*z + 313*x*(y^2)*z + 317*x*y*(z^2)
// symbolically evaluated at (x,y,z) = (1, 2, 3)
void rhome_test_data(
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

// Checks derived formula and computed result against rhome_test_data()
BOOST_AUTO_TEST_CASE( rhome_grad_u )
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

    rhome_test_data(
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

// Checks derived formula and computed result against rhome_test_data()
BOOST_AUTO_TEST_CASE( rhome_div_u )
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

    rhome_test_data(
        rho, grad_rho, div_grad_rho, grad_grad_rho,
        m, div_m, grad_m, div_grad_m, grad_div_m,
        e, grad_e);

    const double div_u = pecos::suzerain::orthonormal::rhome::div_u(
            rho, grad_rho, m, div_m);

    const double ans = -7.8528271460331434259092626891655172292938028930092832721830;
    const double close_enough = std::numeric_limits<double>::epsilon();
    BOOST_CHECK_CLOSE(div_u, ans, close_enough);
}

// Checks derived formula and computed result against rhome_test_data()
BOOST_AUTO_TEST_CASE( rhome_grad_div_u )
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

    rhome_test_data(
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
