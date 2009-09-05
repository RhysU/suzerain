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

    pecos::suzerain::cartesian::rhome::p_T_mu_lambda(
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
    rho = 138.0;

    grad_rho(0) = 150.0;
    grad_rho(1) =  87.0;
    grad_rho(2) =  76.0;

    div_grad_rho = 62.0;

    grad_grad_rho(0,0) = 24.0;
    grad_grad_rho(0,1) = 93.0;
    grad_grad_rho(0,2) = 80.0;
    grad_grad_rho(1,0) = 93.0;
    grad_grad_rho(1,1) = 18.0;
    grad_grad_rho(1,2) = 44.0;
    grad_grad_rho(2,0) = 80.0;
    grad_grad_rho(2,1) = 44.0;
    grad_grad_rho(2,2) = 20.0;

    m(0) =  7.0*sin( 7.0) + 11.0*sin(22.0) + 13.0*sin(39.0);
    m(1) = 17.0*sin(17.0) + 19.0*sin(38.0) + 23.0*sin(69.0);
    m(2) = 31.0*sin(31.0) + 37.0*sin(74.0) + 41.0*sin(123.0);

    div_m = 49.0*cos(7.0) + 361.0*cos(38.0) + 1681.0*cos(123.0);

    grad_m(0,0) =   49.0*cos(  7.0);
    grad_m(0,1) =  121.0*cos( 22.0);
    grad_m(0,2) =  169.0*cos( 39.0);
    grad_m(1,0) =  289.0*cos( 17.0);
    grad_m(1,1) =  361.0*cos( 38.0);
    grad_m(1,2) =  529.0*cos( 69.0);
    grad_m(2,0) =  961.0*cos( 31.0);
    grad_m(2,1) = 1369.0*cos( 74.0);
    grad_m(2,2) = 1681.0*cos(123.0);

    div_grad_m(0) = -  343.0*sin( 7.0) -  1331.0*sin(22.0) -  2197.0*sin( 39.0);
    div_grad_m(1) = - 4913.0*sin(17.0) -  6859.0*sin(38.0) - 12167.0*sin( 69.0);
    div_grad_m(2) = -29791.0*sin(31.0) - 50653.0*sin(74.0) - 68921.0*sin(123.0);

    grad_div_m(0) = -  343.0*sin(  7.0);
    grad_div_m(1) = - 6859.0*sin( 38.0);
    grad_div_m(2) = -68921.0*sin(123.0);

    e = 11328.0;

    grad_e(0) = 13194.0;
    grad_e(1) =  7542.0;
    grad_e(2) =  5678.0;
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
        = pecos::suzerain::cartesian::rhome::grad_u(rho, grad_rho, m, grad_m);

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

    const double div_u = pecos::suzerain::cartesian::rhome::div_u(
            rho, grad_rho, m, div_m);

    const double ans = -7.8528271460331434259092626891655172292938028930092832721830;
    const double close_enough = std::numeric_limits<double>::epsilon();
    BOOST_CHECK_CLOSE(div_u, ans, close_enough);
}
