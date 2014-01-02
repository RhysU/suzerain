/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2009-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

#include <suzerain/lowstorage.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <suzerain/common.hpp>
#include <suzerain/richardson.h>
#include <suzerain/state.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

// Shorthand
namespace mpl = boost::mpl;
using suzerain::contiguous_state;
using suzerain::interleaved_state;
using suzerain::multi_array::ref;
using suzerain::lowstorage::controller;
using suzerain::lowstorage::linear_operator;
using suzerain::lowstorage::method;
using suzerain::lowstorage::method_interface;
using suzerain::lowstorage::multiplicative_operator;
using suzerain::lowstorage::nonlinear_operator;
using suzerain::lowstorage::smr91;
using suzerain::lowstorage::yang11;

// Explicit template instantiation to hopefully speed compilation
template class interleaved_state<3,double>;
template class contiguous_state<3,double>;

// State types to be tested against themselves
typedef mpl::list<
        contiguous_state<3,double>,
        interleaved_state<3,double>
    > state_types;

// Pairs of state types to be tested for interoperability
typedef mpl::list<
        mpl::vector<contiguous_state <3,double>, contiguous_state <3,double> >,
        mpl::vector<interleaved_state<3,double>, interleaved_state<3,double> >,
        mpl::vector<contiguous_state <3,double>, interleaved_state<3,double> >,
        mpl::vector<interleaved_state<3,double>, contiguous_state <3,double> >
    > state_type_pairs;

// Helper method for providing 3D size information
static suzerain::array<std::size_t,3> size3(
    std::size_t x, std::size_t y, std::size_t z)
{
    suzerain::array<std::size_t,3> a = {{ x, y, z }};
    return a;
}

// Shorthand for a double NaN used in may places below
static const double double_NaN = std::numeric_limits<double>::quiet_NaN();

// Nonlinear portion of a hybrid implicit/explicit Riccati operator is the
// right hand side of (d/dt) y = y^2 + b y - a^2 -a b minus the b y portion.
class riccati_nonlinear_operator
    : public nonlinear_operator<ref<double,3> >
{
private:
    const double a;
    const double b;
    const double delta_t;

public:
    riccati_nonlinear_operator(
            const double a,
            const double b,
            const double delta_t = std::numeric_limits<double>::infinity())
        : a(a), b(b), delta_t(delta_t) {};

    virtual std::vector<double> apply_operator(
            const double time,
            ref<double,3>& state,
            const method_interface<double>& method,
            const std::size_t substep_index) const
    {
        SUZERAIN_UNUSED(time);
        SUZERAIN_UNUSED(method);
        SUZERAIN_UNUSED(substep_index);

        typedef contiguous_state<3,double>::index index;
        for (index i = 0; i < (index) state.shape()[0]; ++i) {
            for (index k = 0; k < (index) state.shape()[2]; ++k) {
                for (index j = 0; j < (index) state.shape()[1]; ++j) {
                    double &y = state[i][j][k];
                    y = y*y - a*a - a*b;
                }
            }
        }

        return std::vector<double>(1, delta_t);
    };

};

template< typename StateA, typename StateB = StateA >
class riccati_linear_operator
    : public multiplicative_operator<StateA,StateB>
{
public:
    riccati_linear_operator(
            const double a,
            const double b,
            const double delta_t = std::numeric_limits<double>::infinity())
        : multiplicative_operator<StateA,StateB>(b, delta_t)
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
class cosine_explicit_operator
    : public nonlinear_operator<ref<double,3> >
{
private:
    const double delta_t;

public:
    cosine_explicit_operator(const double delta_t = double_NaN)
        : delta_t(delta_t) { };

    virtual std::vector<double> apply_operator(
            const double time,
            ref<double,3>& state,
            const method_interface<double>& method,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(method);
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

typedef mpl::list< float ,double ,long double > constants_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( convective, T, constants_test_types )
{
    using suzerain::lowstorage::convective_stability_criterion;

    const T pi            = boost::math::constants::pi<T>();
    const T u_x           = - 3;
    const T delta_x       =   5;
    const T u_y           = - 7;
    const T delta_y       =   9;
    const T u_z           =  11;
    const T delta_z       =  13;
    const T evmaxmag_imag =  17;

    { // Ensure <= restriction correct as documented
        const T a = 19;
        const T delta_t = convective_stability_criterion(
                u_x, a, pi/delta_x, u_y, a, pi/delta_y, u_z, a, pi/delta_z,
                evmaxmag_imag);
        const T lhs = pi*(  (std::abs(u_x)+a)/delta_x
                          + (std::abs(u_y)+a)/delta_y
                          + (std::abs(u_z)+a)/delta_z)*delta_t;
        BOOST_CHECK_LE(
                lhs, evmaxmag_imag + 10*std::numeric_limits<T>::epsilon());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( diffusive, T, constants_test_types )
{
    using suzerain::lowstorage::diffusive_stability_criterion;

    const T pi            = boost::math::constants::pi<T>();
    const T delta_x       =  5;
    const T delta_y       =  9;
    const T delta_z       = 13;
    const T Re            = 17;
    const T Pr            = 19;
    const T gamma         = 23;
    const T evmaxmag_real = 29;
    const T nu            = 37;

    { // Ensure <= restriction correct as documented
        const T nu0 = 31;
        BOOST_REQUIRE_LE(nu0, nu); // Legit restriction?
        const T delta_t = diffusive_stability_criterion(
                pi*pi/(delta_x*delta_x),
                pi*pi/(delta_y*delta_y),
                pi*pi/(delta_z*delta_z), Re, Pr, gamma,
                evmaxmag_real, nu, nu0, T(0), T(0));
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
                pi*pi/(delta_x*delta_x),
                pi*pi/(delta_y*delta_y),
                pi*pi/(delta_z*delta_z), Re, Pr, gamma,
                evmaxmag_real, nu, nu0, T(0), T(0));
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


BOOST_AUTO_TEST_SUITE( SMR91_sanity )

BOOST_AUTO_TEST_CASE( name )
{
    const method<smr91,float> m;
    BOOST_CHECK(m.name());
}

typedef mpl::list< float ,double ,long double > constants_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( SMR91_constants, T, constants_test_types )
{
    { // Real-valued constants
        const T close_enough = std::numeric_limits<T>::epsilon();
        const method<smr91,T> m;

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
        const method<smr91,complex> m;

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
        const method<yang11,T> m;

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
        const method<yang11,complex> m;

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

BOOST_AUTO_TEST_SUITE( multiplicative_operator_sanity )

BOOST_AUTO_TEST_CASE( apply_operator )
{
    const method<smr91,double> m;
    typedef multiplicative_operator<contiguous_state<3,double> > op_type;
    const double close_enough = std::numeric_limits<double>::epsilon();

    contiguous_state<3,double> a(size3(1,1,1));
    a[0][0][0] = 1.0;

    op_type op(2.0);
    op.apply_operator(double_NaN, a, m, 0);
    BOOST_CHECK_CLOSE(a[0][0][0], 2.0, close_enough);
    op.apply_operator(double_NaN, a, m, 1);
    BOOST_CHECK_CLOSE(a[0][0][0], 4.0, close_enough);
    op.apply_operator(double_NaN, a, m, 2);
    BOOST_CHECK_CLOSE(a[0][0][0], 8.0, close_enough);

    // Ensure we can instantiate
    op_type unused(2.0);
}

BOOST_AUTO_TEST_CASE( accumulate_mass_plus_scaled_operator )
{
    typedef multiplicative_operator<contiguous_state<3,double> > op_type;
    const double close_enough = std::numeric_limits<double>::epsilon();

    contiguous_state<3,double> a(size3(1,1,1)), b(size3(1,1,1));
    a[0][0][0] = 2.0;
    b[0][0][0] = 3.0;

    op_type op(5.0);
    op.accumulate_mass_plus_scaled_operator(7.0, a, 1.0, b);
    BOOST_CHECK_CLOSE(b[0][0][0], 75.0, close_enough);
    op.accumulate_mass_plus_scaled_operator(0.0, b, 1.0, a);
    BOOST_CHECK_CLOSE(a[0][0][0], 77.0, close_enough);

    // Ensure we catch an operation between two nonconforming states
    contiguous_state<3,double> c(size3(2,1,1));
    BOOST_CHECK_THROW(op.accumulate_mass_plus_scaled_operator(3.0, b, 1.0, c),
                      std::logic_error);
}

BOOST_AUTO_TEST_CASE( invert_mass_plus_scaled_operator )
{
    const method<smr91,double> m;
    typedef multiplicative_operator<contiguous_state<3,double> > op_type;
    const double close_enough = std::numeric_limits<double>::epsilon();

    contiguous_state<3,double> a(size3(1,1,1));
    a[0][0][0] = 2.0;

    op_type op(3.0);
    op.invert_mass_plus_scaled_operator(5.0, a, m);
    BOOST_CHECK_CLOSE(a[0][0][0], 1.0/8.0, close_enough);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( substep_suite )

// Purely explicit Riccati equation nonlinear operator
// is the right hand side of (d/dt) y = y^2 + b y - a^2 -a b
class riccati_explicit_operator
    : public nonlinear_operator<ref<double,3> >
{
private:
    const double a, b, delta_t;

public:
    riccati_explicit_operator(const double a,
                            const double b,
                            const double delta_t = double_NaN)
        : a(a), b(b), delta_t(delta_t) { };

    virtual std::vector<double> apply_operator(
            const double time,
            ref<double,3>& state,
            const method_interface<double>& method,
            const std::size_t substep_index = 0) const
    {
        SUZERAIN_UNUSED(time);
        SUZERAIN_UNUSED(method);
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

BOOST_AUTO_TEST_CASE_TEMPLATE( substep_explicit_time_independent,
                               State, state_types )
{
    using suzerain::lowstorage::substep;
    // See test_lowstorage.sage for manufactured answers

    const double delta_t = 17.0;
    const double close_enough = std::numeric_limits<double>::epsilon()*100;
    const method<smr91,double> m;
    const multiplicative_operator<State> trivial_linop(0);
    const riccati_explicit_operator riccati_op(2, 3);
    State a(size3(2,1,1)), b(size3(2,1,1));

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, trivial_linop, 1.0,
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

        const double delta_t_used = substep(m, trivial_linop, 1.0,
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

        const double delta_t_used = substep(m, trivial_linop, 1.0,
                riccati_op, double_NaN, a, b, delta_t, 2);
        BOOST_CHECK_EQUAL(delta_t, delta_t_used);

        BOOST_CHECK_CLOSE(a[0][0][0],  30.0,       close_enough);
        BOOST_CHECK_CLOSE(a[1][0][0],  60.0,       close_enough);
        BOOST_CHECK_CLOSE(b[0][0][0], 3715.0/12.0, close_enough);
        BOOST_CHECK_CLOSE(b[1][0][0], 8159.0/12.0, close_enough);
    }

    // Requesting an out-of-bounds substep_index should balk
    BOOST_CHECK_THROW(substep(m, trivial_linop, 1.0, riccati_op,
                              double_NaN, a, b, delta_t, 3),
            std::invalid_argument);
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( substep_hybrid_time_independent,
                                State, state_types )
{
    using suzerain::lowstorage::substep;
    // See test_lowstorage.sage for manufactured answers

    const double delta_t = 17.0;
    const double close_enough = std::numeric_limits<double>::epsilon()*500;
    const method<smr91,double> m;
    const riccati_nonlinear_operator nonlinear_op(2, 3);
    const riccati_linear_operator<State> linear_op(2, 3);
    State a(size3(2,1,1)), b(size3(2,1,1));

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

BOOST_AUTO_TEST_CASE_TEMPLATE( substep_explicit_time_dependent,
                               State, state_types )
{
    using suzerain::lowstorage::substep;
    // See test_lowstorage.sage for manufactured answers

    const double delta_t = 17.0;
    const double close_enough = std::numeric_limits<double>::epsilon()*500;
    const double pi = boost::math::constants::pi<double>();
    const double time = pi / 3.0;
    const method<smr91,double> m;
    const multiplicative_operator<State> trivial_linop(0);
    const cosine_explicit_operator cosine_op;
    State a(size3(2,1,1)), b(size3(2,1,1));

    {
        a[0][0][0] =  5.0;
        a[1][0][0] =  7.0;
        b[0][0][0] = 11.0;
        b[1][0][0] = 13.0;

        const double delta_t_used = substep(m, trivial_linop, 1.0,
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

        const double delta_t_used = substep(m, trivial_linop, 1.0,
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

        const double delta_t_used = substep(m, trivial_linop, 1.0,
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


BOOST_AUTO_TEST_SUITE( step_suite )

// Run the explicit timestepper against (d/dt) y = a*y
// where it is expected to be third order.
BOOST_AUTO_TEST_CASE_TEMPLATE( step_explicit_time_independent,
                               StatePair, state_type_pairs )
{
    typedef typename mpl::at<StatePair,mpl::int_<0> >::type   state_a_type;
    typedef typename mpl::at<StatePair,mpl::int_<1> >::type   state_b_type;
    typedef multiplicative_operator<state_a_type,state_b_type> mult_op_type;

    // Fix test problem parameters
    const ExponentialSolution soln(2.0, 1.0);
    const double t_initial = 0.140, t_final = 0.145; // Asymptotic regime

    // Fix method, operators, and storage space
    const method<smr91,double> m;
    const mult_op_type trivial_linop(0);
    state_a_type a(size3(1,1,1));
    state_b_type b(size3(1,1,1));

    // Coarse grid calculation using explicitly provided time step
    const std::size_t coarse_nsteps = 16;
    const double delta_t_coarse = (t_final - t_initial)/coarse_nsteps;
    a[0][0][0] = soln(t_initial);
    {
        const mult_op_type nonlinear_op(soln.a);
        for (std::size_t i = 0; i < coarse_nsteps; ++i) {
            const double delta_t_used = suzerain::lowstorage::step(
                    m, trivial_linop, 1.0, nonlinear_op,
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
        const mult_op_type nonlinear_op(soln.a, delta_t_finer);
        for (std::size_t i = 0; i < finer_nsteps; ++i) {
            const double delta_t_used = suzerain::lowstorage::step(
                    m, trivial_linop, 1.0, nonlinear_op, double_NaN, a, b);
            BOOST_CHECK_EQUAL(delta_t_used, delta_t_finer);
        }
    }
    const double finer_final = a[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-13); // Tolerance found using Octave

    // smr91 is third order against the exponential problem
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
BOOST_AUTO_TEST_CASE_TEMPLATE( step_explicit_time_dependent,
                               StatePair, state_type_pairs )
{
    typedef typename mpl::at<StatePair,mpl::int_<0> >::type   state_a_type;
    typedef typename mpl::at<StatePair,mpl::int_<1> >::type   state_b_type;
    typedef multiplicative_operator<state_a_type,state_b_type> mult_op_type;

    // Fix test problem parameters
    const CosineSolution soln(0.0);
    const double t_initial = 0.000, t_final = 0.0125; // Asymptotic regime

    // Fix method, operators, and storage space
    const method<smr91,double> m;
    const mult_op_type trivial_linop(0);
    state_a_type a(size3(1,1,1));
    state_b_type b(size3(1,1,1));

    // Coarse grid calculation using explicitly provided time step
    const std::size_t coarse_nsteps = 16;
    const double delta_t_coarse = (t_final - t_initial)/coarse_nsteps;
    a[0][0][0] = soln(t_initial);
    {
        const cosine_explicit_operator nonlinear_op(t_final - t_initial);
        double t = t_initial;
        for (std::size_t i = 0; i < coarse_nsteps; ++i) {
            const double delta_t_used = suzerain::lowstorage::step(
                    m, trivial_linop, 1.0, nonlinear_op,
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
        const cosine_explicit_operator nonlinear_op(t_final - t_initial);
        double t = t_initial;
        for (std::size_t i = 0; i < finer_nsteps; ++i) {
            const double delta_t_used = suzerain::lowstorage::step(
                    m, trivial_linop, 1.0, nonlinear_op,
                    t, a, b, delta_t_finer);
            BOOST_CHECK_EQUAL(delta_t_used, delta_t_finer);
            t += delta_t_used;
        }
        BOOST_CHECK_GE(t, t_final); // Some slop allowed
    }
    const double finer_final = a[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-13);

    // smr91 is third order against the cosine problem
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
BOOST_AUTO_TEST_CASE_TEMPLATE( step_hybrid, StatePair, state_type_pairs )
{
    typedef typename mpl::at<StatePair,mpl::int_<0> >::type state_a_type;
    typedef typename mpl::at<StatePair,mpl::int_<1> >::type state_b_type;

    // Fix test problem parameters
    const RiccatiSolution soln(2.0, 2.0, -50.0);
    const double t_initial = 0.140, t_final = 0.145; // Asymptotic regime

    // Fix method, operators, and storage space
    const method<smr91,double> m;
    const riccati_nonlinear_operator nonlinear_op(soln.a, soln.b);
    const riccati_linear_operator<
                state_a_type, state_b_type
            > linear_op(soln.a, soln.b);
    state_a_type a(size3(1,1,1));
    state_b_type b(size3(1,1,1));

    // Coarse grid calculation
    const std::size_t coarse_nsteps = 16;
    a[0][0][0] = soln(t_initial);
    for (std::size_t i = 0; i < coarse_nsteps; ++i) {
        suzerain::lowstorage::step(
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
        suzerain::lowstorage::step(
                m, linear_op, 1.0, nonlinear_op, double_NaN, a, b,
                (t_final - t_initial)/finer_nsteps);
    }
    const double finer_final = a[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-11); // Tolerance found using Octave

    // smr91 is second order against the Riccati problem
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


// Tests for control logic of timecontroller in test_timecontroller.
// Presumably getting lowstorage::controller to type check is the big deal.
// Explicitly instantiate it to ensure the template looks okay.

// controller for interleaved_state
template class controller<
        interleaved_state<3,double>, interleaved_state<3,double>,
        suzerain::lowstorage::delta_t_reducer
    >;
template class controller<
        interleaved_state<3,double>, interleaved_state<3,double>,
        void // Default Reducer behavior
    >;

// controller for contiguous_state
template class controller<
        contiguous_state<3,double>, contiguous_state<3,double>,
        suzerain::lowstorage::delta_t_reducer
    >;
template class controller<
        contiguous_state<3,double>, contiguous_state<3,double>,
        void // Default reducer behavior
    >;

// controller for {Interleaved,Contiguous}State
template class controller<
        interleaved_state<3,double>, contiguous_state<3,double>,
        suzerain::lowstorage::delta_t_reducer
    >;
template class controller<
        interleaved_state<3,double>, contiguous_state<3,double>,
        void // Default reducer behavior
    >;

// controller for {Contiguous,Interleaved}State
template class controller<
        contiguous_state<3,double>, interleaved_state<3,double>,
        suzerain::lowstorage::delta_t_reducer
    >;
template class controller<
        contiguous_state<3,double>, interleaved_state<3,double>,
        void // Default reducer behavior
    >;

BOOST_AUTO_TEST_SUITE( controller_suite )

BOOST_AUTO_TEST_CASE_TEMPLATE ( invoke_make_controller, StatePair, state_type_pairs )
{
    typedef typename mpl::at<StatePair,mpl::int_<0> >::type state_a_type;
    typedef typename mpl::at<StatePair,mpl::int_<1> >::type state_b_type;

    const method<smr91,double> m;
    const multiplicative_operator<state_a_type, state_b_type> trivial_linop(0);
    const riccati_nonlinear_operator riccati_op(2, 3);
    state_a_type a(size3(2,1,1));
    state_b_type b(size3(2,1,1));

    // Compilation and instantiation is half the battle.  Go Joe!
    suzerain::scoped_ptr<suzerain::timecontroller<double> > p(
        make_controller(m, trivial_linop, 1., riccati_op, a, b));

    BOOST_REQUIRE(p);
}

BOOST_AUTO_TEST_SUITE_END()
