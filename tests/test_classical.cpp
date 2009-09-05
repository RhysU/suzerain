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

// Check that the derived formula and computed result are correct against
//      rho = 2*x + 3*y + 5*z
//      m   = {
//               7*sin( 7*x) + 11*sin(11*y) + 13*sin(13*z),
//              17*sin(17*x) + 19*sin(19*y) + 23*sin(23*z),
//              31*sin(31*x) + 37*sin(37*y) + 41*sin(41*z)
//      }
// symbolically evaluated at (x,y,z) = (1, 2, 3)
BOOST_AUTO_TEST_CASE( rhome_grad_u )
{
    const double rho = 23;
    Eigen::Vector3d grad_rho;
    grad_rho << 2.0, 3.0, 5.0;
    Eigen::Vector3d m;
    m <<  7.0*sin( 7.0) + 11.0*sin(22.0) + 13.0*sin(39.0),
         17.0*sin(17.0) + 19.0*sin(38.0) + 23.0*sin(69.0),
         31.0*sin(31.0) + 37.0*sin(74.0) + 41.0*sin(123.0);
    Eigen::Matrix3d grad_m;
    grad_m(0,0) =   49.0*cos(  7.0);
    grad_m(0,1) =  121.0*cos( 22.0);
    grad_m(0,2) =  169.0*cos( 39.0);
    grad_m(1,0) =  289.0*cos( 17.0);
    grad_m(1,1) =  361.0*cos( 38.0);
    grad_m(1,2) =  529.0*cos( 69.0);
    grad_m(2,0) =  961.0*cos( 31.0);
    grad_m(2,1) = 1369.0*cos( 74.0);
    grad_m(2,2) = 1681.0*cos(123.0);

    const Eigen::Matrix3d grad_u
        = pecos::suzerain::cartesian::rhome::grad_u(rho, grad_rho, m, grad_m);

    Eigen::Matrix3d expected;
    expected(0,0) = -14./529.*sin(7.)-22./529.*sin(22.)-26./529.*sin(39.)+49./23.*cos(7.);
    expected(0,1) = -21./529.*sin(7.)-33./529.*sin(22.)-39./529.*sin(39.)+121./23.*cos(22.);
    expected(0,2) = -35./529.*sin(7.)-55./529.*sin(22.)-65./529.*sin(39.)+169./23.*cos(39.);
    expected(1,0) = -34./529.*sin(17.)-38./529.*sin(38.)-2./23.*sin(69.)+289./23.*cos(17.);
    expected(1,1) = -51./529.*sin(17.)-57./529.*sin(38.)-3./23.*sin(69.)+361./23.*cos(38.);
    expected(1,2) = -85./529.*sin(17.)-95./529.*sin(38.)-5./23.*sin(69.)+23.*cos(69.);
    expected(2,0) = -62./529.*sin(31.)-74./529.*sin(74.)-82./529.*sin(123.)+961./23.*cos(31.);
    expected(2,1) = -93./529.*sin(31.)-111./529.*sin(74.)-123./529.*sin(123.)+1369./23.*cos(74.);
    expected(2,2) = -155./529.*sin(31.)-185./529.*sin(74.)-205./529.*sin(123.)+1681./23.*cos(123.);

    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e2;
    BOOST_CHECK_CLOSE(grad_u(0,0), expected(0,0), close_enough);
    BOOST_CHECK_CLOSE(grad_u(0,1), expected(0,1), close_enough);
    BOOST_CHECK_CLOSE(grad_u(0,2), expected(0,2), close_enough);
    BOOST_CHECK_CLOSE(grad_u(1,0), expected(1,0), close_enough);
    BOOST_CHECK_CLOSE(grad_u(1,1), expected(1,1), close_enough);
    BOOST_CHECK_CLOSE(grad_u(1,2), expected(1,2), close_enough);
    BOOST_CHECK_CLOSE(grad_u(2,0), expected(2,0), close_enough);
    BOOST_CHECK_CLOSE(grad_u(2,1), expected(2,1), close_enough);
    BOOST_CHECK_CLOSE(grad_u(2,2), expected(2,2), close_enough);
}

BOOST_AUTO_TEST_CASE( rhome_div_u )
{
    const double rho = 2.0;
    Eigen::Vector3d grad_rho;
    grad_rho << 3.0, 5.0, 7.0;
    Eigen::Vector3d m;
    m << 11.0, 13.0, 17.0;
    const double div_m = 19.0;

    const double div_u = pecos::suzerain::cartesian::rhome::div_u(
            rho, grad_rho, m, div_m);

    const double expected = -0.25*(3.0*11.0+5.0*13.0+7.0*17.0) + 0.5*19.0;
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e-3;
    BOOST_CHECK_CLOSE(div_u, expected, close_enough);
}
