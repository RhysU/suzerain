#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>
#include <Eigen/Core>
#include <limits>

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

BOOST_AUTO_TEST_CASE( rhome_grad_u )
{
    const double rho = 2.0;
    Eigen::Vector3d grad_rho;
    grad_rho << 3.0, 5.0, 7.0;
    Eigen::Vector3d m;
    m << 11.0, 13.0, 17.0;
    Eigen::Matrix3d grad_m;
    grad_m << 19.0, 23.0, 29.0,
              31.0, 37.0, 41.0,
              43.0, 47.0, 53.0;
    Eigen::Matrix3d grad_u;

    pecos::suzerain::cartesian::rhome::grad_u(
            rho, grad_rho, m, grad_m, grad_u);

    Eigen::Matrix3d expected;
    expected(0,0) = 0.5*19.0 - 0.25*3.0*11.0;
    expected(0,1) = 0.5*23.0 - 0.25*3.0*13.0;
    expected(0,2) = 0.5*29.0 - 0.25*3.0*17.0;
    expected(1,0) = 0.5*31.0 - 0.25*5.0*11.0;
    expected(1,1) = 0.5*37.0 - 0.25*5.0*13.0;
    expected(1,2) = 0.5*41.0 - 0.25*5.0*17.0;
    expected(2,0) = 0.5*43.0 - 0.25*7.0*11.0;
    expected(2,1) = 0.5*47.0 - 0.25*7.0*13.0;
    expected(2,2) = 0.5*53.0 - 0.25*7.0*17.0;
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e3;

    BOOST_CHECK_SMALL((grad_u - expected).norm(), close_enough);
}
