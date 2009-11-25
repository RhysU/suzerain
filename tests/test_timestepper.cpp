#include <suzerain/config.h>
#include <suzerain/common.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timestepper.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( SMR91Method )

BOOST_AUTO_TEST_CASE( constants )
{
    using suzerain::timestepper::lowstorage::SMR91Method;

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

BOOST_AUTO_TEST_SUITE( MultiplicativeOperator )

BOOST_AUTO_TEST_CASE( applyOperator )
{
    const double close_enough = std::numeric_limits<double>::epsilon();

    using suzerain::RealState;
    using suzerain::timestepper::lowstorage::MultiplicativeOperator;

    RealState<double> a(1,1,1);
    a.data[0][0][0] = 1.0;

    MultiplicativeOperator<double> op(2.0);
    op.applyOperator(&a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 2.0, close_enough);
    op.applyOperator(&a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 4.0, close_enough);
    op.applyOperator(&a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 8.0, close_enough);
}

BOOST_AUTO_TEST_CASE( accumulateIdentityPlusScaledOperator )
{
    const double close_enough = std::numeric_limits<double>::epsilon();

    using suzerain::RealState;
    using suzerain::timestepper::lowstorage::MultiplicativeOperator;

    RealState<double> a(1,1,1), b(1,1,1), c(2,1,1);
    a.data[0][0][0] = 2.0;
    b.data[0][0][0] = 3.0;

    MultiplicativeOperator<double> op(5.0);
    op.accumulateIdentityPlusScaledOperator(7.0, &a, &b);
    BOOST_CHECK_CLOSE(b.data[0][0][0], 75.0, close_enough);
    op.accumulateIdentityPlusScaledOperator(0.0, &b, &a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 77.0, close_enough);

    // Ensure we catch an operation between two nonconforming states
    BOOST_CHECK_THROW(op.accumulateIdentityPlusScaledOperator(3.0, &b, &c),
                      std::logic_error);
}

BOOST_AUTO_TEST_CASE( invertIdentityPlusScaledOperator )
{
    const double close_enough = std::numeric_limits<double>::epsilon();

    using suzerain::RealState;
    using suzerain::timestepper::lowstorage::MultiplicativeOperator;

    RealState<double> a(1,1,1), b(2,1,1);
    a.data[0][0][0] = 2.0;

    MultiplicativeOperator<double> op(3.0);
    op.invertIdentityPlusScaledOperator(5.0, &a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 1.0/8.0, close_enough);
}

BOOST_AUTO_TEST_SUITE_END()
