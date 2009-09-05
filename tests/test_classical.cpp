#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <limits>
#include <Eigen/Core>
#include <boost/test/included/unit_test.hpp>
#include <suzerain/classical.hpp>


BOOST_AUTO_TEST_CASE( rhome_p_T_mu_lambda )
{
    const double close_enough = std::numeric_limits<double>::epsilon() * 1.0e4;
    using namespace pecos::suzerain::cartesian::rhome;

    const double beta  = 2.0/3.0;
    const double gamma = 1.4;
    const double rho   = 2.0;
    Eigen::Vector3d m;
    m << 3.0, 5.0, 7.0;
    const double e = 163.0;
    double p, T, mu, lambda;

    p_T_mu_lambda(beta, gamma, rho, m, e, p, T, mu, lambda);

    BOOST_CHECK_CLOSE(p,       56.9,             close_enough);
    BOOST_CHECK_CLOSE(T,       39.83,            close_enough);
    BOOST_CHECK_CLOSE(mu,      11.6629085673383, close_enough);
    BOOST_CHECK_CLOSE(lambda,  -7.7752723782255, close_enough);
}
