#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#include <fftw3.h>
#include <boost/test/included/unit_test.hpp>
#include "test_tools.hpp"

template<int NR = 8, int NC = (NR/2) + 1>
class data {
public:
    typedef double                 Real;
    typedef std::complex<Real>     Complex;
    typedef std::valarray<Real>    VectorReal;
    typedef std::valarray<Complex> VectorComplex;

    const Real L;           /* domain size in physical space */
    VectorReal x;           /* discrete grid in physical space */
    VectorReal r;           /* r_i = f(x_i) */
    VectorReal k;           /* k_i = wavenumber i */
    VectorComplex c;        /* c_i = \hat{f}_{k_i} */

    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan_r2c;
    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan_c2r;

    data() throw(std::logic_error)
        :   L(2*M_PI),
            x(NR),
            r(Real(0), NR),
            k(NC),
            c(Complex(0), NC),
            plan_r2c(fftw_plan_dft_r2c_1d(
                       NR,
                       &r[0],
                       reinterpret_cast<fftw_complex*>(&c[0]),
                       FFTW_ESTIMATE),
                    std::ptr_fun(fftw_destroy_plan)),
            plan_c2r(fftw_plan_dft_c2r_1d(
                       NR,
                       reinterpret_cast<fftw_complex*>(&c[0]),
                       &r[0],
                       FFTW_ESTIMATE),
                    std::ptr_fun(fftw_destroy_plan))
    {
        if (!plan_r2c.get()) throw std::logic_error("invalid r2c plan");
        if (!plan_c2r.get()) throw std::logic_error("invalid c2r plan");

        x[0] = 0.0;
        for (std::size_t i = 1; i < NR; ++i) x[i] = x[i-1] + L/NR;
        for (std::size_t i = 0; i < NC; ++i) k[i] = 2*M_PI*i/L;
    }

    void to_wave_space()
    {
        fftw_execute(plan_r2c.get());
        c /= NR;
    }

    void differentiate()
    {
        for (int i = 0; i < NC; ++i) c[i] *= Complex(0,k[i]);
    }

    void to_physical_space()
    {
        fftw_execute(plan_c2r.get());
    }
};

BOOST_AUTO_TEST_CASE( forward_and_backward )
{
    typedef data<> Data;
    Data d;

    d.r = sin(2.0*d.x);
    const Data::VectorReal expected(d.r);

    d.to_wave_space();
    d.to_physical_space();

    check_close_collections(
            &d.r[0], &d.r[d.r.size()],
            &expected[0], &expected[expected.size()],
            std::numeric_limits<Data::Real>::epsilon() * 1.0e+3);
}

BOOST_AUTO_TEST_CASE( poor_mans_dealiasing )
{
    typedef data<8> Data8;
    typedef data<16> Data16;
    Data8 d8;
    Data16 d16;

    d8.r = sin(2.0*d8.x);
    const Data8::VectorReal expected(d8.r);
    d8.to_wave_space();

    for (std::size_t i = 0; i < d8.c.size(); ++i) {
        d16.c[i] = d8.c[i];
    }
    d8.to_physical_space();
    d16.to_physical_space();

    const Data8::Real close_enough =
        std::numeric_limits<Data8::Real>::epsilon() * 1.0e+3;

    check_close_collections(
            &d8.r[0], &d8.r[d8.r.size()],
            &expected[0], &expected[expected.size()],
            close_enough);

    // Make sure we recovered the same data on the 16 mode grid
    Data8::VectorReal error(d16.r[std::slice(0,d8.r.size(),2)]);
    error = abs(error - d8.r);
    BOOST_CHECK_SMALL(error.sum(), close_enough);
}

BOOST_AUTO_TEST_CASE( differentiate )
{
    typedef data<32> Data;
    Data d;

    d.r = sin(2.0*d.x);
    d.to_wave_space();
    d.differentiate();
    d.to_physical_space();

    Data::VectorReal expected = 2.0*cos(2.0*d.x);

    BOOST_CHECK_SMALL(
            (abs(expected-d.r)).sum(),
            std::numeric_limits<Data::Real>::epsilon() * 1.0e+4);
}
