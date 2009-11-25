#include <suzerain/config.h>
#include <suzerain/common.hpp>
#include <suzerain/timestepper.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( SMR91Method )

BOOST_AUTO_TEST_CASE( constants )
{
    using suzerain::timestepper::SMR91Method;

    {
        const float close_enough = std::numeric_limits<float>::epsilon();
        SMR91Method<float> m;
        for (int i = 0; i < m.substeps(); ++i) {
            const float res = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(res, close_enough);
        }
    }

    {
        const double close_enough = std::numeric_limits<double>::epsilon();
        SMR91Method<double> m;
        for (int i = 0; i < m.substeps(); ++i) {
            const double res = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(res, close_enough);
        }
    }

    {
        const long double close_enough
            = std::numeric_limits<long double>::epsilon();
        SMR91Method<long double> m;
        for (int i = 0; i < m.substeps(); ++i) {
            const long double res
                = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(res, close_enough);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
