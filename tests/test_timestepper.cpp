#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timestepper.hpp>
#include <suzerain/richardson.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

// Purely explicit Riccati equation nonlinear operator
// is the right hand side of (d/dt) y = y^2 + b y - a^2 -a b
template< typename FPT >
class RiccatiExplicitOperator
    : public suzerain::timestepper::INonlinearOperator<FPT>
{
private:
    const FPT a;
    const FPT b;
    const FPT delta_t;

public:
    RiccatiExplicitOperator(
            const FPT a,
            const FPT b,
            const FPT delta_t = std::numeric_limits<FPT>::quiet_NaN())
        : a(a), b(b), delta_t(delta_t) { };

    virtual FPT applyOperator(suzerain::IState<FPT> &state,
                              const bool delta_t_requested = false) const
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
        return delta_t;
    };
};

// Nonlinear portion of a hybrid implicit/explicit Riccati operator is the
// right hand side of (d/dt) y = y^2 + b y - a^2 -a b minus the b y portion.
template< typename FPT >
class RiccatiNonlinearOperator
    : public suzerain::timestepper::INonlinearOperator<FPT>
{
private:
    const FPT a;
    const FPT b;
    const FPT delta_t;

public:
    RiccatiNonlinearOperator(
            const FPT a,
            const FPT b,
            const FPT delta_t = std::numeric_limits<FPT>::quiet_NaN())
        : a(a), b(b), delta_t(delta_t) {};

    virtual FPT applyOperator(suzerain::IState<FPT> &state,
                              const bool delta_t_requested = false) const
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
        return delta_t;
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

// Functor returning the solution (d/dt) y = y^2 + b y - a^2 -a b
// where y(t) = a + (-(2*a+b)^(-1) + c*exp(-(2*a+b)*t))^(-1).
template< typename FPT >
struct RiccatiSolution
{
    const FPT a, b, c;

    RiccatiSolution(FPT a, FPT b, FPT c) : a(a), b(b), c(c) {}

    FPT operator()(FPT t) const
    {
        return a + FPT(1)/(FPT(-1)/(2*a+b) + c*exp(-(2*a+b)*t));
    }
};

// Functor returning the solution (d/dt) y = a*y where y(t) = y0 * exp(a*t)
template< typename FPT >
struct ExponentialSolution
{
    const FPT a, y0;

    ExponentialSolution(FPT a, FPT y0) : a(a), y0(y0) {}

    FPT operator()(FPT t) const { return y0*exp(a*t); }
};


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

BOOST_AUTO_TEST_CASE( accumulateMassPlusScaledOperator )
{
    const double close_enough = std::numeric_limits<double>::epsilon();

    using suzerain::RealState;
    using suzerain::timestepper::lowstorage::MultiplicativeOperator;

    RealState<double> a(1,1,1), b(1,1,1), c(2,1,1);
    a.data[0][0][0] = 2.0;
    b.data[0][0][0] = 3.0;

    MultiplicativeOperator<double> op(5.0);
    op.accumulateMassPlusScaledOperator(7.0, a, b);
    BOOST_CHECK_CLOSE(b.data[0][0][0], 75.0, close_enough);
    op.accumulateMassPlusScaledOperator(0.0, b, a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 77.0, close_enough);

    // Ensure we catch an operation between two nonconforming states
    BOOST_CHECK_THROW(op.accumulateMassPlusScaledOperator(3.0, b, c),
                      std::logic_error);
}

BOOST_AUTO_TEST_CASE( invertMassPlusScaledOperator )
{
    const double close_enough = std::numeric_limits<double>::epsilon();

    using suzerain::RealState;
    using suzerain::timestepper::lowstorage::MultiplicativeOperator;

    RealState<double> a(1,1,1), b(2,1,1);
    a.data[0][0][0] = 2.0;

    MultiplicativeOperator<double> op(3.0);
    op.invertMassPlusScaledOperator(5.0, a);
    BOOST_CHECK_CLOSE(a.data[0][0][0], 1.0/8.0, close_enough);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( substep )

BOOST_AUTO_TEST_CASE( substep_explicit )
{
    // See test_timestepper.sage for manufactured answers

    using suzerain::RealState;
    using suzerain::timestepper::lowstorage::MultiplicativeOperator;
    using suzerain::timestepper::lowstorage::SMR91Method;
    using suzerain::timestepper::lowstorage::substep;

    const double close_enough = std::numeric_limits<double>::epsilon()*100;
    const SMR91Method<double> m;
    const MultiplicativeOperator<double>  trivial_linear_op(0);
    const RiccatiExplicitOperator<double> riccati_op(2, 3);
    RealState<double> a(2,1,1), b(2,1,1);

    {
        a.data[0][0][0] =  5.0;
        a.data[1][0][0] =  7.0;
        b.data[0][0][0] = 11.0;
        b.data[1][0][0] = 13.0;

        substep(m, trivial_linear_op, riccati_op, 17.0, a, b, 0);

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

        substep(m, trivial_linear_op, riccati_op, 17.0, a, b, 1);

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

        substep(m, trivial_linear_op, riccati_op, 17.0, a, b, 2);

        BOOST_CHECK_CLOSE(a.data[0][0][0],  30.0,       close_enough);
        BOOST_CHECK_CLOSE(a.data[1][0][0],  60.0,       close_enough);
        BOOST_CHECK_CLOSE(b.data[0][0][0], 3715.0/12.0, close_enough);
        BOOST_CHECK_CLOSE(b.data[1][0][0], 8159.0/12.0, close_enough);
    }

    // Requesting an out-of-bounds substep_index should balk
    BOOST_CHECK_THROW(substep(m, trivial_linear_op, riccati_op, 17.0, a, b, 3),
                      std::invalid_argument);
}

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


BOOST_AUTO_TEST_SUITE( step_delta_t_provided )

// Run the explicit timestepper against (d/dt) y = a*y
// where it is expected to be third order.
BOOST_AUTO_TEST_CASE( step_explicit )
{
    // Fix test problem parameters
    const ExponentialSolution<double> soln(2.0, 1.0);
    const double t_initial = 0.140, t_final = 0.145; // Asymptotic regime

    // Fix method, operators, and storage space
    using namespace suzerain::timestepper;
    const lowstorage::SMR91Method<double> m;
    const lowstorage::MultiplicativeOperator<double> trivial_linear_op(0);
    const lowstorage::MultiplicativeOperator<double> nonlinear_op(soln.a);
    suzerain::RealState<double> a(1,1,1), b(1,1,1);

    // Coarse grid calculation
    const std::size_t coarse_nsteps = 16;
    a.data[0][0][0] = soln(t_initial);
    for (std::size_t i = 0; i < coarse_nsteps; ++i) {
        suzerain::timestepper::lowstorage::step(
                m, trivial_linear_op, nonlinear_op,
                (t_final - t_initial)/coarse_nsteps, a, b);
    }
    const double coarse_final = a.data[0][0][0];
    const double coarse_error = fabs(coarse_final - soln(t_final));
    BOOST_CHECK_SMALL(coarse_error, 1.0e-12); // Tolerance found using Octave

    // Finer grid calculation
    const std::size_t finer_nsteps = 2*coarse_nsteps;
    a.data[0][0][0] = soln(t_initial);
    for (std::size_t i = 0; i < finer_nsteps; ++i) {
        suzerain::timestepper::lowstorage::step(
                m, trivial_linear_op, nonlinear_op,
                (t_final - t_initial)/finer_nsteps, a, b);
    }
    const double finer_final = a.data[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-13); // Tolerance found using Octave

    // SMR91 is third order against the exponential problem
    const double expected_order = 2.95; // allows for floating point losses
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    BOOST_CHECK(!boost::math::isnan(observed_order));
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
        for (int i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
            gsl_vector_set(k,0,i+1);
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
    const RiccatiSolution<double> soln(2.0, 2.0, -50.0);
    const double t_initial = 0.140, t_final = 0.145; // Asymptotic regime

    // Fix method, operators, and storage space
    const suzerain::timestepper::lowstorage::SMR91Method<double> m;
    const RiccatiNonlinearOperator<double> nonlinear_op(soln.a, soln.b);
    const RiccatiLinearOperator<double>    linear_op(soln.a, soln.b);
    suzerain::RealState<double> a(1,1,1), b(1,1,1);

    // Coarse grid calculation
    const std::size_t coarse_nsteps = 16;
    a.data[0][0][0] = soln(t_initial);
    for (std::size_t i = 0; i < coarse_nsteps; ++i) {
        suzerain::timestepper::lowstorage::step(
                m, linear_op, nonlinear_op,
                (t_final - t_initial)/coarse_nsteps, a, b);
    }
    const double coarse_final = a.data[0][0][0];
    const double coarse_error = fabs(coarse_final - soln(t_final));
    BOOST_CHECK_SMALL(coarse_error, 1.0e-10); // Tolerance found using Octave

    // Finer grid calculation
    const std::size_t finer_nsteps = 2*coarse_nsteps;
    a.data[0][0][0] = soln(t_initial);
    for (std::size_t i = 0; i < finer_nsteps; ++i) {
        suzerain::timestepper::lowstorage::step(
                m, linear_op, nonlinear_op,
                (t_final - t_initial)/finer_nsteps, a, b);
    }
    const double finer_final = a.data[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-11); // Tolerance found using Octave

    // SMR91 is second order against the Riccati problem
    const double expected_order = 1.98; // allows for floating point losses
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    BOOST_CHECK(!boost::math::isnan(observed_order));
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
        for (int i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
            gsl_vector_set(k,0,i+1);
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


BOOST_AUTO_TEST_SUITE( step_delta_t_computed )

// Run the timestepper purely explicitly against (d/dt) y = a*y
// where it is expected to be third order.
BOOST_AUTO_TEST_CASE( step_explicit )
{
    // Fix test problem parameters
    const ExponentialSolution<double> soln(2.0, 1.0);
    const double t_initial = 0.140, t_final = 0.145; // Asymptotic regime

    // Fix method, operators, and storage space
    using namespace suzerain::timestepper;
    const lowstorage::SMR91Method<double> m;
    const lowstorage::MultiplicativeOperator<double> trivial_linear_op(0);
    suzerain::RealState<double> a(1,1,1), b(1,1,1);

    // Coarse grid calculation
    const std::size_t coarse_nsteps = 16;
    a.data[0][0][0] = soln(t_initial);
    {
        // Nonlinear operator provides timestep size information
        const lowstorage::MultiplicativeOperator<double> nonlinear_op(
                soln.a, (t_final - t_initial)/coarse_nsteps);
        for (std::size_t i = 0; i < coarse_nsteps; ++i) {
            suzerain::timestepper::lowstorage::step(
                    m, trivial_linear_op, nonlinear_op, a, b);
        }
    }
    const double coarse_final = a.data[0][0][0];
    const double coarse_error = fabs(coarse_final - soln(t_final));
    BOOST_CHECK_SMALL(coarse_error, 1.0e-12); // Tolerance found using Octave

    // Finer grid calculation
    const std::size_t finer_nsteps = 2*coarse_nsteps;
    a.data[0][0][0] = soln(t_initial);
    {
        // Nonlinear operator provides timestep size information
        const lowstorage::MultiplicativeOperator<double> nonlinear_op(
                soln.a, (t_final - t_initial)/finer_nsteps);
        for (std::size_t i = 0; i < finer_nsteps; ++i) {
            suzerain::timestepper::lowstorage::step(
                    m, trivial_linear_op, nonlinear_op, a, b);
        }
    }
    const double finer_final = a.data[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-13); // Tolerance found using Octave

    // SMR91 is third order against the exponential problem
    const double expected_order = 2.95; // allows for floating point losses
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    BOOST_CHECK(!boost::math::isnan(observed_order));
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
        for (int i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
            gsl_vector_set(k,0,i+1);
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
    const RiccatiSolution<double> soln(2.0, 2.0, -50.0);
    const double t_initial = 0.140, t_final = 0.145; // Asymptotic regime

    // Fix method, operators, and storage space
    const suzerain::timestepper::lowstorage::SMR91Method<double> m;
    const RiccatiLinearOperator<double>    linear_op(soln.a, soln.b);
    suzerain::RealState<double> a(1,1,1), b(1,1,1);

    // Coarse grid calculation
    const std::size_t coarse_nsteps = 16;
    a.data[0][0][0] = soln(t_initial);
    {
        // Nonlinear operator provides timestep size information
        const RiccatiNonlinearOperator<double> nonlinear_op(
                soln.a, soln.b, (t_final - t_initial)/coarse_nsteps);
        for (std::size_t i = 0; i < coarse_nsteps; ++i) {
            suzerain::timestepper::lowstorage::step(
                    m, linear_op, nonlinear_op, a, b);
        }
    }
    const double coarse_final = a.data[0][0][0];
    const double coarse_error = fabs(coarse_final - soln(t_final));
    BOOST_CHECK_SMALL(coarse_error, 1.0e-10); // Tolerance found using Octave

    // Finer grid calculation
    const std::size_t finer_nsteps = 2*coarse_nsteps;
    a.data[0][0][0] = soln(t_initial);
    {
        // Nonlinear operator provides timestep size information
        const RiccatiNonlinearOperator<double> nonlinear_op(
                soln.a, soln.b, (t_final - t_initial)/finer_nsteps);
        for (std::size_t i = 0; i < finer_nsteps; ++i) {
            suzerain::timestepper::lowstorage::step(
                    m, linear_op, nonlinear_op, a, b);
        }
    }
    const double finer_final = a.data[0][0][0];
    const double finer_error = fabs(finer_final - soln(t_final));
    BOOST_CHECK_SMALL(finer_error, 1.0e-11); // Tolerance found using Octave

    // SMR91 is second order against the Riccati problem
    const double expected_order = 1.98; // allows for floating point losses
    const double observed_order
        = log(coarse_error/finer_error)/log(finer_nsteps/coarse_nsteps);
    BOOST_CHECK(!boost::math::isnan(observed_order));
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
        for (int i = 0;
             i < sizeof(richardson_h_error)/sizeof(richardson_h_error[0]);
             ++i) {
            gsl_matrix_set(data, 0, 0, coarse_final);
            gsl_matrix_set(data, 0, 1, finer_final);
            gsl_vector_set(k,0,i+1);
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
