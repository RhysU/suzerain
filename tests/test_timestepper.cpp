#include <suzerain/config.h>
#include <suzerain/common.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timestepper.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( SMR91Method )

BOOST_AUTO_TEST_CASE( name )
{
    const suzerain::timestepper::lowstorage::SMR91Method<float> m;
    BOOST_CHECK(m.name());
}

BOOST_AUTO_TEST_CASE( constants )
{
    using suzerain::timestepper::lowstorage::SMR91Method;

    {
        const float close_enough = std::numeric_limits<float>::epsilon();
        const SMR91Method<float> m;
        for (int i = 0; i < m.substeps(); ++i) {
            const float res = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(res, close_enough);
        }
    }

    {
        const double close_enough = std::numeric_limits<double>::epsilon();
        const SMR91Method<double> m;
        for (int i = 0; i < m.substeps(); ++i) {
            const double res = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(res, close_enough);
        }
    }

    {
        const long double close_enough
            = std::numeric_limits<long double>::epsilon();
        const SMR91Method<long double> m;
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
    op.applyOperator(a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 2.0, close_enough);
    op.applyOperator(a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 4.0, close_enough);
    op.applyOperator(a);
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
    op.accumulateIdentityPlusScaledOperator(7.0, a, b);
    BOOST_CHECK_CLOSE(b.data[0][0][0], 75.0, close_enough);
    op.accumulateIdentityPlusScaledOperator(0.0, b, a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 77.0, close_enough);

    // Ensure we catch an operation between two nonconforming states
    BOOST_CHECK_THROW(op.accumulateIdentityPlusScaledOperator(3.0, b, c),
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
    op.invertIdentityPlusScaledOperator(5.0, a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 1.0/8.0, close_enough);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( substep )

// Purely explicit Riccati equation nonlinear operator
// is the right hand side of (d/dt) y = y^2 + b y - a^2 -a b
template< typename FPT >
class RiccatiExplicitOperator
    : public suzerain::timestepper::INonlinearOperator<FPT>
{
private:
    const FPT a;
    const FPT b;

public:
    RiccatiExplicitOperator(const FPT a, const FPT b) : a(a), b(b) { };

    virtual void applyOperator(suzerain::IState<FPT> &state) const
                               throw(std::exception)
    {
        suzerain::RealState<FPT> &realstate
            = dynamic_cast<suzerain::RealState<FPT>&>(state);
        typedef typename suzerain::RealState<FPT>::index index;
        for (index k = 0; k < realstate.data.shape()[2]; ++k)
            for (index j = 0; j < realstate.data.shape()[1]; ++j)
                for (index i = 0; i < realstate.data.shape()[0]; ++i) {
                    FPT &y = realstate.data[i][j][k];
                    y = y*y + b*y - a*a - a*b;
                }
    };
};

BOOST_AUTO_TEST_CASE( substep_explicit )
{
    // See test_timestepper.sage for manufactured answers

    using suzerain::RealState;
    using suzerain::timestepper::lowstorage::SMR91Method;
    using suzerain::timestepper::lowstorage::substep;

    const double close_enough = std::numeric_limits<double>::epsilon()*100;
    const SMR91Method<double> m;
    const RiccatiExplicitOperator<double> riccati_op(2, 3);
    RealState<double> a(2,1,1), b(2,1,1);

    {
        a.data[0][0][0] =  5.0;
        a.data[1][0][0] =  7.0;
        b.data[0][0][0] = 11.0;
        b.data[1][0][0] = 13.0;

        substep(m, riccati_op, 17.0, a, b, 0);

        BOOST_CHECK_CLOSE(a.data[0][0][0],  30.0, close_enough);
        BOOST_CHECK_CLOSE(a.data[1][0][0],  60.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[0][0][0], 277.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[1][0][0], 551.0, close_enough);
    }

    {
        a.data[0][0][0] =  5.0;
        a.data[1][0][0] =  7.0;
        b.data[0][0][0] = 11.0;
        b.data[1][0][0] = 13.0;

        substep(m, riccati_op, 17.0, a, b, 1);

        BOOST_CHECK_CLOSE(a.data[0][0][0],    30.0,      close_enough);
        BOOST_CHECK_CLOSE(a.data[1][0][0],    60.0,      close_enough);
        BOOST_CHECK_CLOSE(b.data[0][0][0],  9871.0/60.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[1][0][0], 22163.0/60.0, close_enough);
    }

    {
        a.data[0][0][0] =  5.0;
        a.data[1][0][0] =  7.0;
        b.data[0][0][0] = 11.0;
        b.data[1][0][0] = 13.0;

        substep(m, riccati_op, 17.0, a, b, 2);

        BOOST_CHECK_CLOSE(a.data[0][0][0],  30.0,       close_enough);
        BOOST_CHECK_CLOSE(a.data[1][0][0],  60.0,       close_enough);
        BOOST_CHECK_CLOSE(b.data[0][0][0], 3715.0/12.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[1][0][0], 8159.0/12.0, close_enough);
    }

    // Requesting an out-of-bounds substep_index should balk
    BOOST_CHECK_THROW(substep(m, riccati_op, 17.0, a, b, 3),
                      std::invalid_argument);
}

// Nonlinear portion of a hybrid implicit/explicit Riccati operator is the
// right hand side of (d/dt) y = y^2 + b y - a^2 -a b minus the b y portion.
template< typename FPT >
class RiccatiNonlinearOperator
    : public suzerain::timestepper::INonlinearOperator<FPT>
{
private:
    const FPT a;
    const FPT b;

public:
    RiccatiNonlinearOperator(const FPT a, const FPT b) : a(a), b(b) {};

    virtual void applyOperator(suzerain::IState<FPT> &state) const
                               throw(std::exception)
    {
        suzerain::RealState<FPT> &realstate
            = dynamic_cast<suzerain::RealState<FPT>&>(state);
        typedef typename suzerain::RealState<FPT>::index index;
        for (index k = 0; k < realstate.data.shape()[2]; ++k)
            for (index j = 0; j < realstate.data.shape()[1]; ++j)
                for (index i = 0; i < realstate.data.shape()[0]; ++i) {
                    FPT &y = realstate.data[i][j][k];
                    y = y*y - a*a - a*b;
                }
    };
};

template< typename FPT >
class RiccatiLinearOperator
    : public suzerain::timestepper::lowstorage::MultiplicativeOperator<FPT>
{
public:
    RiccatiLinearOperator(const FPT a, const FPT b)
        : suzerain::timestepper::lowstorage::MultiplicativeOperator<FPT>(b) {};
};

BOOST_AUTO_TEST_CASE( substep_hybrid )
{
    // See test_timestepper.sage for manufactured answers

    using suzerain::RealState;
    using suzerain::timestepper::lowstorage::SMR91Method;
    using suzerain::timestepper::lowstorage::substep;

    const double close_enough = std::numeric_limits<double>::epsilon()*500;
    const SMR91Method<double> m;
    const RiccatiNonlinearOperator<double> nonlinear_op(2, 3);
    const RiccatiLinearOperator<double>    linear_op(2, 3);
    RealState<double> a(2,1,1), b(2,1,1);

    {
        a.data[0][0][0] =  5.0;
        a.data[1][0][0] =  7.0;
        b.data[0][0][0] = 11.0;
        b.data[1][0][0] = 13.0;

        substep(m, linear_op, nonlinear_op, 17.0, a, b, 0);

        BOOST_CHECK_CLOSE( a.data[0][0][0],            15.0, close_enough);
        BOOST_CHECK_CLOSE( a.data[1][0][0],            39.0, close_enough);
        BOOST_CHECK_CLOSE( b.data[0][0][0], -34885.0/1727.0, close_enough);
        BOOST_CHECK_CLOSE( b.data[1][0][0], -74951.0/1727.0, close_enough);
    }

    {
        a.data[0][0][0] =  5.0;
        a.data[1][0][0] =  7.0;
        b.data[0][0][0] = 11.0;
        b.data[1][0][0] = 13.0;

        substep(m, linear_op, nonlinear_op, 17.0, a, b, 1);

        BOOST_CHECK_CLOSE(a.data[0][0][0],            15.0, close_enough);
        BOOST_CHECK_CLOSE(a.data[1][0][0],            39.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[0][0][0], -   61.0/  15.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[1][0][0], -23263.0/1155.0, close_enough);
    }

    {
        a.data[0][0][0] =  5.0;
        a.data[1][0][0] =  7.0;
        b.data[0][0][0] = 11.0;
        b.data[1][0][0] = 13.0;

        substep(m, linear_op, nonlinear_op, 17.0, a, b, 2);

        BOOST_CHECK_CLOSE(a.data[0][0][0],       15.0, close_enough);
        BOOST_CHECK_CLOSE(a.data[1][0][0],       39.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[0][0][0], -193.0/9.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[1][0][0], -566.0/9.0, close_enough);
    }

    // Requesting an out-of-bounds substep_index should balk
    BOOST_CHECK_THROW(substep(m, linear_op, nonlinear_op, 17.0, a, b, 3),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
