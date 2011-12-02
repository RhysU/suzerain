#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timestepper.hpp>
#include <suzerain/richardson.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

// Shorthand
using suzerain::InterleavedState;
using suzerain::ContiguousState;
using suzerain::timestepper::INonlinearOperator;
using suzerain::timestepper::lowstorage::ILinearOperator;
using suzerain::timestepper::lowstorage::MultiplicativeOperator;
using suzerain::timestepper::lowstorage::SMR91Method;
using suzerain::timestepper::lowstorage::Yang11Method;
using suzerain::timestepper::lowstorage::LowStorageTimeController;
typedef MultiplicativeOperator<ContiguousState<3,double> >
    MultiplicativeOperatorD3;

// Explicit template instantiation to hopefully speed compilation
template class InterleavedState<3,double>;
template class ContiguousState<3,double>;

// Helper method for providing 3D size information
static boost::array<std::size_t,3> size3(
    std::size_t x, std::size_t y, std::size_t z)
{
    boost::array<std::size_t,3> a = {{ x, y, z }};
    return a;
}

// Shorthand for a double NaN used in may places below
static const double double_NaN = std::numeric_limits<double>::quiet_NaN();

// Nonlinear portion of a hybrid implicit/explicit Riccati operator is the
// right hand side of (d/dt) y = y^2 + b y - a^2 -a b minus the b y portion.
class RiccatiNonlinearOperator
    : public INonlinearOperator<ContiguousState<3,double> >
{
private:
    const double a;
    const double b;
    const double delta_t;

public:
    RiccatiNonlinearOperator(
            const double a,
            const double b,
            const double delta_t = std::numeric_limits<double>::infinity())
        : a(a), b(b), delta_t(delta_t) {};

    virtual std::vector<double> applyOperator(
            const double time,
            ContiguousState<3,double>& state,
            const double evmaxmag_real,
            const double evmaxmag_imag,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(time);
        SUZERAIN_UNUSED(evmaxmag_real);
        SUZERAIN_UNUSED(evmaxmag_imag);
        SUZERAIN_UNUSED(substep_index);

        ContiguousState<3,double>& s
            = dynamic_cast<ContiguousState<3,double>&>(state);

        typedef ContiguousState<3,double>::index index;
        for (index i = 0; i < (index) s.shape()[0]; ++i) {
            for (index k = 0; k < (index) s.shape()[2]; ++k) {
                for (index j = 0; j < (index) s.shape()[1]; ++j) {
                    double &y = s[i][j][k];
                    y = y*y - a*a - a*b;
                }
            }
        }

        return std::vector<double>(1, delta_t);
    };

};

class RiccatiLinearOperator
    : public MultiplicativeOperator<ContiguousState<3,double> >
{
public:
    RiccatiLinearOperator(
            const double a,
            const double b,
            const double delta_t = std::numeric_limits<double>::infinity())
        : MultiplicativeOperator<ContiguousState<3,double> >(b, delta_t)
    {
        SUZERAIN_UNUSED(a);
    }
};

// Functor returning the solution (d/dt) y = y^2 + b y - a^2 -a b
// where y(t) = a + (-(2*a+b)^(-1) + c*exp(-(2*a+b)*t))^(-1).
struct RiccatiSolution
{
    const double a, b, c;

    RiccatiSolution(double a, double b, double c) : a(a), b(b), c(c) {}

    double operator()(double t) const
    {
        return a + double(1)/(double(-1)/(2*a+b) + c*std::exp(-(2*a+b)*t));
    }
};

// Functor returning the solution (d/dt) y = a*y where y(t) = y0 * exp(a*t)
struct ExponentialSolution
{
    const double a, y0;

    ExponentialSolution(double a, double y0) : a(a), y0(y0) {}

    double operator()(double t) const { return y0*std::exp(a*t); }
};

// Purely explicit, time-dependent operator for (d/dt) y = cos(t);
class CosineExplicitOperator
    : public INonlinearOperator<ContiguousState<3,double> >
{
private:
    const double delta_t;

public:
    CosineExplicitOperator(
            const double delta_t = std::numeric_limits<double>::quiet_NaN())
        : delta_t(delta_t) { };

    virtual std::vector<double> applyOperator(
            const double time,
            ContiguousState<3,double> & state,
            const double evmaxmag_real,
            const double evmaxmag_imag,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(evmaxmag_real);
        SUZERAIN_UNUSED(evmaxmag_imag);
        SUZERAIN_UNUSED(substep_index);

        for (std::size_t i = 0; i < state.shape()[0]; ++i) {
            for (std::size_t k = 0; k < state.shape()[2]; ++k) {
                for (std::size_t j = 0; j < state.shape()[1]; ++j) {
                    state[i][j][k] = std::cos(time);
                }
            }
        }

        return std::vector<double>(1, delta_t);
    }
};

// Functor returning the solution (d/dt) y = cos(t) where y(t) = y0 + sin(t)
struct CosineSolution
{
    const double y0;

    CosineSolution(double y0) : y0(y0) {}

    double operator()(double t) const { return y0 + std::sin(t); }
};


BOOST_AUTO_TEST_SUITE( timestep_stability )

typedef boost::mpl::list< float ,double ,long double > constants_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( convective, T, constants_test_types )
{
    using suzerain::timestepper::convective_stability_criterion;

    const T pi            = boost::math::constants::pi<T>();
    const T u_x           = - 3;
    const T delta_x       =   5;
    const T u_y           = - 7;
    const T delta_y       =   9;
    const T u_z           =  11;
    const T delta_z       =  13;
    const T evmaxmag_imag =  17;

    { // Ensure default a behavior correct
        const T with_a = convective_stability_criterion(
                u_x, T(1)/delta_x, u_y, T(1)/delta_y, u_z, T(1)/delta_z,
                evmaxmag_imag, T(0));
        const T without_a = convective_stability_criterion(
                u_x, T(1)/delta_x, u_y, T(1)/delta_y, u_z, T(1)/delta_z,
                evmaxmag_imag);
        BOOST_CHECK_EQUAL(with_a, without_a);
    }

    { // Ensure <= restriction correct as documented
        const T a = 19;
        const T delta_t = convective_stability_criterion(
                u_x, T(1)/delta_x, u_y, T(1)/delta_y, u_z, T(1)/delta_z,
                evmaxmag_imag, a);
        const T lhs = pi*(  (std::abs(u_x)+a)/delta_x
                          + (std::abs(u_y)+a)/delta_y
                          + (std::abs(u_z)+a)/delta_z)*delta_t;
        BOOST_CHECK_LE(
                lhs, evmaxmag_imag + 10*std::numeric_limits<T>::epsilon());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( diffusive, T, constants_test_types )
{
    using suzerain::timestepper::diffusive_stability_criterion;

    const T pi            = boost::math::constants::pi<T>();
    const T delta_x       =  5;
    const T delta_y       =  9;
    const T delta_z       = 13;
    const T Re            = 17;
    const T Pr            = 19;
    const T gamma         = 23;
    const T evmaxmag_real = 29;
    const T nu            = 37;

    { // Ensure default a behavior correct
        const T with_nu0 = diffusive_stability_criterion(
                T(1)/delta_x, T(1)/delta_y, T(1)/delta_z, Re, Pr, gamma,
                evmaxmag_real, nu, T(0));
        const T without_nu0 = diffusive_stability_criterion(
                T(1)/delta_x, T(1)/delta_y, T(1)/delta_z, Re, Pr, gamma,
                evmaxmag_real, nu);
        BOOST_CHECK_EQUAL(with_nu0, without_nu0);
    }

    { // Ensure <= restriction correct as documented
        const T nu0 = 31;
        BOOST_REQUIRE_LE(nu0, nu); // Legit restriction?
        const T delta_t = diffusive_stability_criterion(
                T(1)/delta_x, T(1)/delta_y, T(1)/delta_z, Re, Pr, gamma,
                evmaxmag_real, nu, nu0);
        using std::abs;
        const T coeff = std::max(gamma*(abs(nu-nu0))/(Re*Pr), abs(nu-nu0)/Re);
        const T lhs   = coeff*pi*pi*(   T(1)/(delta_x*delta_x)
                                      + T(1)/(delta_y*delta_y)
                                      + T(1)/(delta_z*delta_z))*delta_t;
        BOOST_CHECK_LE(
                lhs, evmaxmag_real + 30*std::numeric_limits<T>::epsilon());
    }

    { // Ensure we handle nu < nu0 correctly
        const T nu0 = 11;
        BOOST_REQUIRE_LE(nu0, nu); // Legit restriction?
        const T delta_t = diffusive_stability_criterion(
                T(1)/delta_x, T(1)/delta_y, T(1)/delta_z, Re, Pr, gamma,
                evmaxmag_real, nu, nu0);
        using std::abs;
        const T coeff = std::max(gamma*(abs(nu-nu0))/(Re*Pr), abs(nu-nu0)/Re);
        const T lhs   = coeff*pi*pi*(   T(1)/(delta_x*delta_x)
                                      + T(1)/(delta_y*delta_y)
                                      + T(1)/(delta_z*delta_z))*delta_t;
        BOOST_CHECK_LE(
                lhs, evmaxmag_real + 30*std::numeric_limits<T>::epsilon());
    }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( SMR91Method_sanity )

BOOST_AUTO_TEST_CASE( name )
{
    const SMR91Method<float> m;
    BOOST_CHECK(m.name());
}

typedef boost::mpl::list< float ,double ,long double > constants_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( SMR91_constants, T, constants_test_types )
{
    { // Real-valued constants
        const T close_enough = std::numeric_limits<T>::epsilon();
        const SMR91Method<T> m;

        // Step size consistency
        for (std::size_t i = 0; i < m.substeps(); ++i) {
            const T res = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(res, close_enough);
        }

        // Time offset consistency
        for (std::size_t i = 0; i < m.substeps(); ++i) {
            T res = m.eta(i);
            for (std::size_t j = i; j-- > 0 ;) {
                res -= m.alpha(j);
                res -= m.beta(j);
            }
            BOOST_CHECK_SMALL(res, close_enough);
        }

        // Consistency of iota with eta
        BOOST_CHECK_EQUAL(m.iota(0), T(1));
        for (std::size_t i = 0; i < m.substeps()-1; ++i) {
            BOOST_CHECK_CLOSE(  (m.eta(i+1) - m.eta(i))
                               /(m.eta(i+1) - m.eta(0)),
                                 m.iota(i), /* Noisy */ 135*close_enough );
        }
        BOOST_CHECK_CLOSE(m.iota(m.substeps()-1),
                          1-m.eta(m.substeps()-1),
                          /* Noisy */ 135*close_enough);
    }

    { // Complex-valued constants
        const T close_enough = std::numeric_limits<T>::epsilon();
        typedef typename std::complex<T> complex;
        const SMR91Method<complex> m;

        // Step size consistency
        for (std::size_t i = 0; i < m.substeps(); ++i) {
            const complex res
                = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(std::abs(res), close_enough);
        }

        // Time offset consistency
        for (std::size_t i = 0; i < m.substeps(); ++i) {
            T res = m.eta(i);
            for (std::size_t j = i; j-- > 0 ;) {
                res -= m.alpha(j);
                res -= m.beta(j);
            }
            BOOST_CHECK_SMALL(res, close_enough);
        }

        // Consistency of iota with eta
        BOOST_CHECK_EQUAL(m.iota(0), T(1));
        for (std::size_t i = 0; i < m.substeps()-1; ++i) {
            BOOST_CHECK_CLOSE(  (m.eta(i+1) - m.eta(i))
                               /(m.eta(i+1) - m.eta(0)),
                                 m.iota(i),
                                 /* Noisy */ 135*close_enough );
        }
        BOOST_CHECK_CLOSE(m.iota(m.substeps()-1),
                          1-m.eta(m.substeps()-1),
                          /* Noisy */ 100*close_enough);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( Yang11_constants, T, constants_test_types )
{
    { // Real-valued constants
        const T close_enough = std::numeric_limits<T>::epsilon();
        const Yang11Method<T> m;

        // Step size consistency
        for (std::size_t i = 0; i < m.substeps(); ++i) {
            const T res = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(res, close_enough);
        }

        // Time offset consistency
        for (std::size_t i = 0; i < m.substeps(); ++i) {
            T res = m.eta(i);
            for (std::size_t j = i; j-- > 0 ;) {
                res -= m.alpha(j);
                res -= m.beta(j);
            }
            BOOST_CHECK_SMALL(res, close_enough);
        }

        // Consistency of iota with eta
        BOOST_CHECK_EQUAL(m.iota(0), T(1));
        for (std::size_t i = 0; i < m.substeps()-1; ++i) {
            BOOST_CHECK_CLOSE(  (m.eta(i+1) - m.eta(i))
                               /(m.eta(i+1) - m.eta(0)),
                                 m.iota(i),
                                 /* Noisy */ 135*close_enough );
        }
        BOOST_CHECK_CLOSE(m.iota(m.substeps()-1),
                          1-m.eta(m.substeps()-1),
                          /* Noisy */ 135*close_enough);
    }

    { // Complex-valued constants
        const T close_enough = std::numeric_limits<T>::epsilon();
        typedef typename std::complex<T> complex;
        const Yang11Method<complex> m;

        // Step size consistency
        for (std::size_t i = 0; i < m.substeps(); ++i) {
            const complex res
                = m.alpha(i) + m.beta(i) - m.gamma(i) - m.zeta(i);
            BOOST_CHECK_SMALL(std::abs(res), close_enough);
        }

        // Time offset consistency
        for (std::size_t i = 0; i < m.substeps(); ++i) {
            T res = m.eta(i);
            for (std::size_t j = i; j-- > 0 ;) {
                res -= m.alpha(j);
                res -= m.beta(j);
            }
            BOOST_CHECK_SMALL(res, close_enough);
        }

        // Consistency of iota with eta
        BOOST_CHECK_EQUAL(m.iota(0), T(1));
        for (std::size_t i = 0; i < m.substeps()-1; ++i) {
            BOOST_CHECK_CLOSE(  (m.eta(i+1) - m.eta(i))
                               /(m.eta(i+1) - m.eta(0)),
                                 m.iota(i),
                                 /* Noisy */ 135*close_enough );
        }
        BOOST_CHECK_CLOSE(m.iota(m.substeps()-1),
                          1-m.eta(m.substeps()-1),
                          /* Noisy */ 135*close_enough);
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( MultiplicativeOperator_sanity )

BOOST_AUTO_TEST_CASE( applyOperator )
{
    const double close_enough = std::numeric_limits<double>::epsilon();

    ContiguousState<3,double> a(size3(1,1,1));
    a[0][0][0] = 1.0;

    MultiplicativeOperatorD3 op(2.0);
    op.applyOperator(double_NaN, a, double_NaN, double_NaN);
    BOOST_CHECK_CLOSE(a[0][0][0], 2.0, close_enough);
    op.applyOperator(double_NaN, a, double_NaN, double_NaN);
    BOOST_CHECK_CLOSE(a[0][0][0], 4.0, close_enough);
    op.applyOperator(double_NaN, a, double_NaN, double_NaN);
    BOOST_CHECK_CLOSE(a[0][0][0], 8.0, close_enough);

    // Ensure we can instantiate
    MultiplicativeOperator<ContiguousState<3,double> > unused(2.0);
}

BOOST_AUTO_TEST_CASE( accumulateMassPlusScaledOperator )
{
    const double close_enough = std::numeric_limits<double>::epsilon();

    ContiguousState<3,double> a(size3(1,1,1)), b(size3(1,1,1));
    a[0][0][0] = 2.0;
    b[0][0][0] = 3.0;

    MultiplicativeOperatorD3 op(5.0);
    op.accumulateMassPlusScaledOperator(7.0, a, 1.0, b);
    BOOST_CHECK_CLOSE(b[0][0][0], 75.0, close_enough);
    op.accumulateMassPlusScaledOperator(0.0, b, 1.0, a);
    BOOST_CHECK_CLOSE(a[0][0][0], 77.0, close_enough);

    // Ensure we catch an operation between two nonconforming states
    ContiguousState<3,double> c(size3(2,1,1));
    BOOST_CHECK_THROW(op.accumulateMassPlusScaledOperator(3.0, b, 1.0, c),
                      std::logic_error);
}

BOOST_AUTO_TEST_CASE( invertMassPlusScaledOperator )
{
    const double close_enough = std::numeric_limits<double>::epsilon();

    ContiguousState<3,double> a(size3(1,1,1));
    a[0][0][0] = 2.0;

    MultiplicativeOperatorD3 op(3.0);
    op.invertMassPlusScaledOperator(5.0, a);
    BOOST_CHECK_CLOSE(a[0][0][0], 1.0/8.0, close_enough);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( substep_suite )

// Purely explicit Riccati equation nonlinear operator
// is the right hand side of (d/dt) y = y^2 + b y - a^2 -a b
class RiccatiExplicitOperator
    : public INonlinearOperator<ContiguousState<3,double> >
{
private:
    const double a, b, delta_t;

public:
    RiccatiExplicitOperator(
            const double a,
            const double b,
            const double delta_t = std::numeric_limits<double>::quiet_NaN())
        : a(a), b(b), delta_t(delta_t) { };

    virtual std::vector<double> applyOperator(
            const double time,
            ContiguousState<3,double> & state,
            const double evmaxmag_real,
            const double evmaxmag_imag,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(time);
        SUZERAIN_UNUSED(evmaxmag_real);
        SUZERAIN_UNUSED(evmaxmag_imag);
        SUZERAIN_UNUSED(substep_index);

        for (std::size_t i = 0; i < state.shape()[0]; ++i) {
            for (std::size_t k = 0; k < state.shape()[2]; ++k) {
                for (std::size_t j = 0; j < state.shape()[1]; ++j) {
                    double &y = state[i][j][k];
                    y = y*y + b*y - a*a - a*b;
                }
            }
        }

        return std::vector<double>(1, delta_t);
    }
};

BOOST_AUTO_TEST_CASE( substep_explicit_time_independent )
{
    using suzerain::timestepper::lowstorage::substep;
    // See test_timestepper.sage for manufactured answers

    const double delta_t = 17.0;
    const double close_enough = std::numeric_limits<double>::epsilon()*100;
    const SMR91Method<double> m;
    const MultiplicativeOperatorD3 trivial_linear_op(0);
    const RiccatiExplicitOperator riccati_op(2, 3);
    ContiguousState<3,double> a(size3(2,1,1)), b(size3(2,1,1));

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, trivial_linear_op, 1.0,
                riccati_op, double_NaN, a, b, delta_t, 0);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],  30.0, close_enough);
        BOOST_CHECK_CLOSE(a[1][0][0],  60.0, close_enough);
        BOOST_CHECK_CLOSE(b[0][0][0], 277.0, close_enough);
        BOOST_CHECK_CLOSE(b[1][0][0], 551.0, close_enough);
    }

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, trivial_linear_op, 1.0,
                riccati_op, double_NaN, a, b, delta_t, 1);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],    30.0,      close_enough);
        BOOST_CHECK_CLOSE(a[1][0][0],    60.0,      close_enough);
        BOOST_CHECK_CLOSE(b[0][0][0],  9871.0/60.0, close_enough);
        BOOST_CHECK_CLOSE(b[1][0][0], 22163.0/60.0, close_enough);
    }

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, trivial_linear_op, 1.0,
                riccati_op, double_NaN, a, b, delta_t, 2);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],  30.0,       close_enough);
        BOOST_CHECK_CLOSE(a[1][0][0],  60.0,       close_enough);
        BOOST_CHECK_CLOSE(b[0][0][0], 3715.0/12.0, close_enough);
        BOOST_CHECK_CLOSE(b[1][0][0], 8159.0/12.0, close_enough);
    }

    // Requesting an out-of-bounds substep_index should balk
    BOOST_CHECK_THROW(substep(m, trivial_linear_op, 1.0, riccati_op,
                              double_NaN, a, b, delta_t, 3),
            std::invalid_argument);
}

BOOST_AUTO_TEST_CASE( substep_hybrid_time_independent )
{
    using suzerain::timestepper::lowstorage::substep;
    // See test_timestepper.sage for manufactured answers

    const double delta_t = 17.0;
    const double close_enough = std::numeric_limits<double>::epsilon()*500;
    const SMR91Method<double> m;
    const RiccatiNonlinearOperator nonlinear_op(2, 3);
    const RiccatiLinearOperator    linear_op(2, 3);
    ContiguousState<3,double> a(size3(2,1,1)), b(size3(2,1,1));

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, linear_op, 1.0,
                nonlinear_op, double_NaN, a, b, delta_t, 0);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE( a[0][0][0],            15.0, close_enough);
        BOOST_CHECK_CLOSE( a[1][0][0],            39.0, close_enough);
        BOOST_CHECK_CLOSE( b[0][0][0], -34885.0/1727.0, close_enough);
        BOOST_CHECK_CLOSE( b[1][0][0], -74951.0/1727.0, close_enough);
    }

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, linear_op, 1.0,
                nonlinear_op, double_NaN, a, b, delta_t, 1);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],            15.0, close_enough);
        BOOST_CHECK_CLOSE(a[1][0][0],            39.0, close_enough);
        BOOST_CHECK_CLOSE(b[0][0][0], -   61.0/  15.0, close_enough);
        BOOST_CHECK_CLOSE(b[1][0][0], -23263.0/1155.0, close_enough);
    }

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, linear_op, 1.0,
                nonlinear_op, double_NaN, a, b, delta_t, 2);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],       15.0, close_enough);
        BOOST_CHECK_CLOSE(a[1][0][0],       39.0, close_enough);
        BOOST_CHECK_CLOSE(b[0][0][0], -193.0/9.0, close_enough);
        BOOST_CHECK_CLOSE(b[1][0][0], -566.0/9.0, close_enough);
    }

    // Requesting an out-of-bounds substep_index should balk
    BOOST_CHECK_THROW(substep(m, linear_op, 1.0, nonlinear_op,
                              double_NaN, a, b, delta_t, 3),
            std::invalid_argument);
}

BOOST_AUTO_TEST_CASE( substep_explicit_time_dependent )
{
    using suzerain::timestepper::lowstorage::substep;
    // See test_timestepper.sage for manufactured answers

    const double delta_t = 17.0;
    const double close_enough = std::numeric_limits<double>::epsilon()*500;
    const double pi = boost::math::constants::pi<double>();
    const double time = pi / 3.0;
    const SMR91Method<double> m;
    const MultiplicativeOperatorD3 trivial_linear_op(0);
    const CosineExplicitOperator cosine_op;
    ContiguousState<3,double> a(size3(2,1,1)), b(size3(2,1,1));

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, trivial_linear_op, 1.0,
                cosine_op, time, a, b, delta_t, 0);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],   1.0/ 2.0, close_enough);
        BOOST_CHECK_CLOSE(a[1][0][0],   1.0/ 2.0, close_enough);
        BOOST_CHECK_CLOSE(b[0][0][0], 143.0/15.0, close_enough);
        BOOST_CHECK_CLOSE(b[1][0][0], 173.0/15.0, close_enough);
    }

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, trivial_linear_op, 1.0,
                cosine_op, time, a, b, delta_t, 1);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],
                std::cos(1.0/3.0*pi + 136.0/15.0), close_enough);

        BOOST_CHECK_CLOSE(a[1][0][0],
                std::cos(1.0/3.0*pi + 136.0/15.0), close_enough);

        BOOST_CHECK_CLOSE(b[0][0][0],
                85.0/12.0*std::cos(1.0/3.0*pi + 136.0/15.0) - 2879.0/60.0,
                close_enough);

        BOOST_CHECK_CLOSE(b[1][0][0],
                85.0/12.0*std::cos(1.0/3.0*pi + 136.0/15.0) - 3337.0/60.0,
                close_enough);
    }

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, trivial_linear_op, 1.0,
                cosine_op, time, a, b, delta_t, 2);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],
                std::cos(1.0/3.0*pi + 34.0/3.0),
                close_enough);
        BOOST_CHECK_CLOSE(a[1][0][0],
                std::cos(1.0/3.0*pi + 34.0/3.0),
                close_enough);
        BOOST_CHECK_CLOSE(b[0][0][0],
                51.0/4.0*std::cos(1.0/3.0*pi + 34.0/3.0) - 875.0/12.0,
                close_enough);
        BOOST_CHECK_CLOSE(b[1][0][0],
                51.0/4.0*std::cos(1.0/3.0*pi + 34.0/3.0) - 1021.0/12.0,
                close_enough);
    }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( noninterleaved_storage )

// Run the explicit timestepper against (d/dt) y = a*y
// where it is expected to be third order.
BOOST_AUTO_TEST_CASE( step_explicit_time_independent )
{
    // Fix test problem parameters
    const ExponentialSolution soln(2.0, 1.0);
    const double t_initial = 0.140, t_final = 0.145; // Asymptotic regime

    // Fix method, operators, and storage space
    const SMR91Method<double> m;
    const MultiplicativeOperatorD3 trivial_linear_op(0);
    ContiguousState<3,double> a(size3(1,1,1)), b(size3(1,1,1));

    // Coarse grid calculation using explicitly provided time step
    const std::size_t coarse_nsteps = 16;
    const double delta_t_coarse = (t_final - t_initial)/coarse_nsteps;
    a[0][0][0] = soln(t_initial);
    {
        const MultiplicativeOperatorD3 nonlinear_op(soln.a);
        for (std::size_t i = 0; i < coarse_nsteps; ++i) {
            const double delta_t_used = suzerain::timestepper::lowstorage::step(
                    m, trivial_linear_op, 1.0, nonlinear_op,
                    double_NaN, a, b, delta_t_coarse);
            BOOST_CHECK_EQUAL(delta_t_used, delta_t_coarse);
        }
    }
    const double coarse_final = a[0][0][0];
    const double coarse_error = fabs(coarse_final - soln(t_final));
    BOOST_CHECK_SMALL(coarse_error, 1.0e-12); // Tolerance found using Octave

    // Finer grid calculation using "dynamically" computed time step
    const std::size_t finer_nsteps = 2*coarse_nsteps;
    const double delta_t_finer = (t_final - t_initial)/finer_nsteps;
    a[0][0][0] = soln(t_initial);
    {
        const MultiplicativeOperatorD3 nonlinear_op(soln.a, delta_t_finer);
        for (std::size_t i = 0; i < finer_nsteps; ++i) {
            const double delta_t_used = suzerain::timestepper::lowstorage::step(
                    m, trivial_linear_op, 1.0, nonlinear_op, double_NaN, a, b);
            BOOST_CHECK_EQUAL(delta_t_used, delta_t_finer);
        }
    }
    const double finer_final = a[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-13); // Tolerance found using Octave

    // SMR91 is third order against the exponential problem
    const double expected_order = 2.95; // allows for floating point losses
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    BOOST_CHECK(!(boost::math::isnan)(observed_order));
    BOOST_CHECK_GE(observed_order, expected_order);

    // Richardson extrapolation should show h^3 term elimination gives a better
    // result than h^1, h^2, h^4, etc...
    {
        gsl_matrix * data = gsl_matrix_alloc(1,2);
        gsl_vector * k = gsl_vector_alloc(1);
        gsl_matrix * normtable = gsl_matrix_alloc(data->size2, data->size2);
        gsl_vector * exact = gsl_vector_alloc(1);
        gsl_vector_set(exact, 0, soln(t_final));

        double richardson_h_error[4];
        for (std::size_t i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
#pragma warning(push,disable:810 2259)
            gsl_vector_set(k,0,i+1);
#pragma warning(pop)
            suzerain_richardson_extrapolation(
                    data, finer_nsteps/coarse_nsteps,
                    k, normtable, exact);
            richardson_h_error[i] = gsl_matrix_get(
                    normtable, normtable->size2-1, normtable->size2-1);
        }

        gsl_matrix_free(normtable);
        gsl_vector_free(exact);
        gsl_vector_free(k);
        gsl_matrix_free(data);

        BOOST_CHECK_LT(richardson_h_error[2], richardson_h_error[0]);
        BOOST_CHECK_LT(richardson_h_error[2], richardson_h_error[1]);
        BOOST_CHECK_LT(richardson_h_error[2], richardson_h_error[3]);
    }
}

// Run the explicit timestepper against (d/dt) y = cos(t)
// where it is expected to be third order.
BOOST_AUTO_TEST_CASE( step_explicit_time_dependent )
{
    // Fix test problem parameters
    const CosineSolution soln(0.0);
    const double t_initial = 0.000, t_final = 0.0125; // Asymptotic regime

    // Fix method, operators, and storage space
    const SMR91Method<double> m;
    const MultiplicativeOperatorD3 trivial_linear_op(0);
    ContiguousState<3,double> a(size3(1,1,1)), b(size3(1,1,1));

    // Coarse grid calculation using explicitly provided time step
    const std::size_t coarse_nsteps = 16;
    const double delta_t_coarse = (t_final - t_initial)/coarse_nsteps;
    a[0][0][0] = soln(t_initial);
    {
        const CosineExplicitOperator nonlinear_op(t_final - t_initial);
        double t = t_initial;
        for (std::size_t i = 0; i < coarse_nsteps; ++i) {
            const double delta_t_used = suzerain::timestepper::lowstorage::step(
                    m, trivial_linear_op, 1.0, nonlinear_op,
                    t, a, b, delta_t_coarse);
            BOOST_CHECK_EQUAL(delta_t_used, delta_t_coarse);
            t += delta_t_used;
        }
        BOOST_CHECK_GE(t, t_final); // Some slop allowed
    }
    const double coarse_final = a[0][0][0];
    const double coarse_error = fabs(coarse_final - soln(t_final));
    BOOST_CHECK_SMALL(coarse_error, 1.0e-12);

    // Finer grid calculation using explicitly provided time step
    const std::size_t finer_nsteps = 2*coarse_nsteps;
    const double delta_t_finer = (t_final - t_initial)/finer_nsteps;
    a[0][0][0] = soln(t_initial);
    {
        const CosineExplicitOperator nonlinear_op(t_final - t_initial);
        double t = t_initial;
        for (std::size_t i = 0; i < finer_nsteps; ++i) {
            const double delta_t_used = suzerain::timestepper::lowstorage::step(
                    m, trivial_linear_op, 1.0, nonlinear_op,
                    t, a, b, delta_t_finer);
            BOOST_CHECK_EQUAL(delta_t_used, delta_t_finer);
            t += delta_t_used;
        }
        BOOST_CHECK_GE(t, t_final); // Some slop allowed
    }
    const double finer_final = a[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-13);

    // SMR91 is third order against the cosine problem
    const double expected_order = 2.95; // allows for floating point losses
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    BOOST_CHECK(!(boost::math::isnan)(observed_order));
    BOOST_CHECK_GE(observed_order, expected_order);

    // Richardson extrapolation should show h^3 term elimination gives a better
    // result than h^1, h^2, h^4, etc...
    {
        gsl_matrix * data = gsl_matrix_alloc(1,2);
        gsl_vector * k = gsl_vector_alloc(1);
        gsl_matrix * normtable = gsl_matrix_alloc(data->size2, data->size2);
        gsl_vector * exact = gsl_vector_alloc(1);
        gsl_vector_set(exact, 0, soln(t_final));

        double richardson_h_error[4];
        for (std::size_t i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
#pragma warning(push,disable:810 2259)
            gsl_vector_set(k,0,i+1);
#pragma warning(pop)
            suzerain_richardson_extrapolation(
                    data, finer_nsteps/coarse_nsteps,
                    k, normtable, exact);
            richardson_h_error[i] = gsl_matrix_get(
                    normtable, normtable->size2-1, normtable->size2-1);
        }

        gsl_matrix_free(normtable);
        gsl_vector_free(exact);
        gsl_vector_free(k);
        gsl_matrix_free(data);

        BOOST_CHECK_LT(richardson_h_error[2], richardson_h_error[0]);
        BOOST_CHECK_LT(richardson_h_error[2], richardson_h_error[1]);
        BOOST_CHECK_LT(richardson_h_error[2], richardson_h_error[3]);
    }
}

// Run the hybrid implicit/explicit timestepper against a Riccati problem
// where it is expected to be second order.
BOOST_AUTO_TEST_CASE( step_hybrid )
{
    // Fix test problem parameters
    const RiccatiSolution soln(2.0, 2.0, -50.0);
    const double t_initial = 0.140, t_final = 0.145; // Asymptotic regime

    // Fix method, operators, and storage space
    const SMR91Method<double> m;
    const RiccatiNonlinearOperator nonlinear_op(soln.a, soln.b);
    const RiccatiLinearOperator    linear_op(soln.a, soln.b);
    ContiguousState<3,double> a(size3(1,1,1)), b(size3(1,1,1));

    // Coarse grid calculation
    const std::size_t coarse_nsteps = 16;
    a[0][0][0] = soln(t_initial);
    for (std::size_t i = 0; i < coarse_nsteps; ++i) {
        suzerain::timestepper::lowstorage::step(
                m, linear_op, 1.0, nonlinear_op, double_NaN, a, b,
                (t_final - t_initial)/coarse_nsteps);
    }
    const double coarse_final = a[0][0][0];
    const double coarse_error = fabs(coarse_final - soln(t_final));
    BOOST_CHECK_SMALL(coarse_error, 1.0e-10); // Tolerance found using Octave

    // Finer grid calculation
    const std::size_t finer_nsteps = 2*coarse_nsteps;
    a[0][0][0] = soln(t_initial);
    for (std::size_t i = 0; i < finer_nsteps; ++i) {
        suzerain::timestepper::lowstorage::step(
                m, linear_op, 1.0, nonlinear_op, double_NaN, a, b,
                (t_final - t_initial)/finer_nsteps);
    }
    const double finer_final = a[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-11); // Tolerance found using Octave

    // SMR91 is second order against the Riccati problem
    const double expected_order = 1.98; // allows for floating point losses
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    BOOST_CHECK(!(boost::math::isnan)(observed_order));
    BOOST_CHECK_GE(observed_order, expected_order);

    // Richardson extrapolation should show h^2 term elimination gives a better
    // result than h^1, h^3, etc...
    {
        gsl_matrix * data = gsl_matrix_alloc(1,2);
        gsl_vector * k = gsl_vector_alloc(1);
        gsl_matrix * normtable = gsl_matrix_alloc(data->size2, data->size2);
        gsl_vector * exact = gsl_vector_alloc(1);
        gsl_vector_set(exact, 0, soln(t_final));

        double richardson_h_error[4];
        for (std::size_t i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
#pragma warning(push,disable:810 2259)
            gsl_vector_set(k,0,i+1);
#pragma warning(pop)
            suzerain_richardson_extrapolation(
                    data, finer_nsteps/coarse_nsteps,
                    k, normtable, exact);
            richardson_h_error[i] = gsl_matrix_get(
                    normtable, normtable->size2-1, normtable->size2-1);
        }

        gsl_matrix_free(normtable);
        gsl_vector_free(exact);
        gsl_vector_free(k);
        gsl_matrix_free(data);

        BOOST_CHECK_LT(richardson_h_error[1], richardson_h_error[0]);
        BOOST_CHECK_LT(richardson_h_error[1], richardson_h_error[2]);
        BOOST_CHECK_LT(richardson_h_error[1], richardson_h_error[3]);
    }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( mixed_storage )

// FIXME: Implement explicit/hybrid tests against Interleaved/Noninterleaved
// See Redmine ticket #1194

BOOST_AUTO_TEST_SUITE_END()


// Tests for control logic of LowStorageTimeController in test_timecontroller.
// Presumably getting LowStorageTimeController to type check is the big deal.
// Explicitly instantiate it to ensure the template looks okay.
// FIXME Instantiate LowStorageTimeController on mixed state
template class LowStorageTimeController<
        InterleavedState<3,double>, InterleavedState<3,double>,
        suzerain::timestepper::DeltaTReducer
    >;
template class LowStorageTimeController<
        InterleavedState<3,double>, InterleavedState<3,double>,
        void // Default Reducer behavior
    >;
template class LowStorageTimeController<
        ContiguousState<3,double>, ContiguousState<3,double>,
        suzerain::timestepper::DeltaTReducer
    >;
template class LowStorageTimeController<
        ContiguousState<3,double>, ContiguousState<3,double>,
        void // Default reducer behavior
    >;

BOOST_AUTO_TEST_SUITE( low_storage_controller_suite )

BOOST_AUTO_TEST_CASE ( make_controller )
{
    const SMR91Method<double> m;
    const MultiplicativeOperatorD3 trivial_linear_op(0);
    const RiccatiNonlinearOperator riccati_op(2, 3);
    ContiguousState<3,double> a(size3(2,1,1)), b(size3(2,1,1));

    // Compilation and instantiation is half the battle.  Go Joe!
    boost::scoped_ptr<suzerain::timestepper::TimeController<double> > p(
        make_LowStorageTimeController(m, trivial_linear_op,
                                      1.0, riccati_op, a, b));

    BOOST_REQUIRE(p);
}

BOOST_AUTO_TEST_SUITE_END()
