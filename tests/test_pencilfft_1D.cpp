//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#include <suzerain/pencilfft.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <fftw3.h>

#include <suzerain/common.hpp>

#include "test_tools.hpp"

#pragma warning(disable:1418)

using namespace suzerain;

template<class MultiArray>
void debug_dump(const std::string &prefix, const MultiArray &x)
{
    BOOST_STATIC_ASSERT(MultiArray::dimensionality == 1);
    using namespace std;

    cout << prefix;
    copy(x.begin(), x.end(),
         ostream_iterator<typename MultiArray::element>(std::cout, " "));
    cout << endl;
}

void debug_dump(const std::string &prefix,
               const fftw_complex * begin,
               const fftw_complex * end) {
    using namespace std;

    cout << prefix;
    while (begin != end) {
        cout << "("
             << (*begin)[0]
             << ","
             << (*begin)[1]
             << ")";
        begin++;
        if (begin != end) cout << ", ";
    }
    cout << endl;
}

// Helper function that kicks the tires of a 1D c2c transform
template<class ComplexMultiArray1, class ComplexMultiArray2>
void symmetry_1D_complex_forward(ComplexMultiArray1 &in,
                                 ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    typedef typename pencilfft::internal::transform_traits<
        typename ComplexMultiArray1::element>::real_type real_type;
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const real_type close_enough
        = std::numeric_limits<real_type>::epsilon()*10*NR*NR*NR;
    const periodic_function<real_type,int> pf(NR, (NR+1)/2);

    // Load a real-valued function into the input array and transform it
    fill_with_NaN(in);
    fill_with_NaN(out);
    for (int i = 0; i < NR; ++i) {
        suzerain::complex::assign_complex(in[i], pf.physical(i), 0.0);
    }
    pencilfft::forward_c2c(0, in, out);

    // Real input should exhibit conjugate symmetry in wave space...
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) { // ...up to grid modes
        real_type a_real, a_imag, b_real, b_imag;
        suzerain::complex::assign_components(a_real, a_imag, out[i]);
        suzerain::complex::assign_components(b_real, b_imag, out[NC-i]);
        BOOST_REQUIRE_CLOSE(a_real,  b_real, close_enough);
        BOOST_REQUIRE_CLOSE(a_imag, -b_imag, close_enough);

        // We should also see the expected frequency content
        const std::complex<real_type> expected = pf.wave(i);
        BOOST_REQUIRE_CLOSE(a_real, expected.real(), sqrt(close_enough));
        BOOST_REQUIRE_CLOSE(a_imag, expected.imag(), sqrt(close_enough));
    }
    // Ensure we see the expected zero mode magnitude
    {
        real_type z_real, z_imag;
        suzerain::complex::assign_components(z_real, z_imag, out[0]);
        const std::complex<real_type> expected = pf.wave(0);
        BOOST_REQUIRE_CLOSE(z_real, expected.real(), close_enough);
        BOOST_REQUIRE_SMALL(z_imag, close_enough);
    }

    // Load an imaginary-valued function into the input array and transform it
    fill_with_NaN(in);
    fill_with_NaN(out);
    for (int i = 0; i < NR; ++i) {
        suzerain::complex::assign_complex(in[i], 0.0, pf.physical(i));
    }
    pencilfft::forward_c2c(0, in, out);

    // Imaginary input should exhibit a similar symmetry in wave space...
    // Re(X_k) = - Re(X_{N-k}), Im(X_k) = Im(X_{N-k})
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) { // ...up to grid modes
        real_type a_real, a_imag, b_real, b_imag;
        suzerain::complex::assign_components(a_real, a_imag, out[i]);
        suzerain::complex::assign_components(b_real, b_imag, out[NC-i]);
        BOOST_REQUIRE_CLOSE(a_real, -b_real, close_enough);
        BOOST_REQUIRE_CLOSE(a_imag,  b_imag, close_enough);

        // We should also see the expected frequency content
        std::complex<real_type> expected
            = std::complex<real_type>(0,1) * pf.wave(i);
        BOOST_REQUIRE_CLOSE(a_real, expected.real(), sqrt(close_enough));
        BOOST_REQUIRE_CLOSE(a_imag, expected.imag(), sqrt(close_enough));
    }
    // Ensure we see the expected zero mode magnitude
    {
        real_type z_real, z_imag;
        suzerain::complex::assign_components(z_real, z_imag, out[0]);
        std::complex<real_type> expected = pf.wave(0);
        BOOST_REQUIRE_SMALL(z_real, close_enough);
        BOOST_REQUIRE_CLOSE(z_imag, expected.real(), close_enough);
    }
}

// Helper function that kicks the tires of a 1D r2c transform
template<class RealMultiArray, class ComplexMultiArray>
void check_1D_forward_r2c(RealMultiArray &in, ComplexMultiArray &out)
{
    BOOST_STATIC_ASSERT(RealMultiArray::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray::dimensionality == 1);
    typedef typename pencilfft::internal::transform_traits<
        typename ComplexMultiArray::element>::real_type real_type;
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const real_type close_enough
        = std::numeric_limits<real_type>::epsilon()*10*NR*NR*NR;
    const periodic_function<real_type,int> pf(NR, (NR+1)/2);

    // Load a real-valued function into the input array and transform it
    fill_with_NaN(out);
    for (int i = 0; i < NR; ++i) {
        in[i] = pf.physical(i);
    }
    pencilfft::forward_r2c(0, in, out);

    // We should see the expected frequency content up to grid modes
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) {
        real_type z_real, z_imag;
        suzerain::complex::assign_components(z_real, z_imag, out[i]);
        const std::complex<real_type> expected = pf.wave(i);
        BOOST_REQUIRE_CLOSE(z_real, expected.real(), sqrt(close_enough));
        BOOST_REQUIRE_CLOSE(z_imag, expected.imag(), sqrt(close_enough));
    }
    // Ensure we see the expected zero mode magnitude
    {
        real_type z_real, z_imag;
        suzerain::complex::assign_components(z_real, z_imag, out[0]);
        const std::complex<real_type> expected = pf.wave(0);
        BOOST_REQUIRE_CLOSE(z_real, expected.real(), close_enough);
        BOOST_REQUIRE_SMALL(z_imag, close_enough);
    }
}

// Compare our results with FFTW's directly computed result
template<class ComplexMultiArray1, class ComplexMultiArray2>
void compare_1D_complex_forward(ComplexMultiArray1 &in,
                                ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    typedef typename pencilfft::internal::transform_traits<
        typename ComplexMultiArray1::element>::real_type real_type;
    const int NC = in.shape()[0];
    const int NR = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*5e2*NC*NC;
    const periodic_function<real_type,int> pf1(NR, (NR+1)/2, 3.0/M_PI);
    const periodic_function<real_type,int> pf2(NR, (NR+1)/2, 5.0/M_PI);

    // Plan before loading in the data since planning overwrites in
    suzerain::shared_array<fftw_complex> buffer(
        static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex)*NR)),
        std::ptr_fun(fftw_free));
    BOOST_REQUIRE(buffer.get() != NULL);
    typename ComplexMultiArray1::element * const in_data
        = in.origin() + std::inner_product(
            in.index_bases(), in.index_bases()+1, in.strides(), 0);
    suzerain::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan(
            fftw_plan_many_dft(1,
                               &NR,
                               1,
                               reinterpret_cast<fftw_complex*>(in_data),
                               NULL,
                               in.strides()[0],
                               NR,
                               buffer.get(),
                               NULL,
                               1,
                               NR,
                               FFTW_FORWARD,
                               FFTW_PRESERVE_INPUT),
            std::ptr_fun(fftw_destroy_plan));
    BOOST_REQUIRE(plan.get() != NULL);

    // Load a complex-valued function into the input array
    fill_with_NaN(in);
    fill_with_NaN(out);
    for (int i = 0; i < NR; ++i) {
        suzerain::complex::assign_complex(
                in[i], pf1.physical(i), pf2.physical(i));
    }

    // Transform input FFTW directly and also our wrapper
    fftw_execute(plan.get());  // Important to be first for in == out
    pencilfft::forward_c2c(0, in, out);

    // Ensure we got the same result, up to scaling differences
    for (int i = 0; i < NR; ++i) {
        double out_real, out_imag;
        suzerain::complex::assign_components(
                out_real, out_imag, out[i]);
        // FFTW with a stride gives a different result than with stride 1
        // BOOST_REQUIRE_EQUAL would be nice, but it fails here
        if (fabs(buffer[i][0]) < close_enough) {
            BOOST_REQUIRE_SMALL(out_real*NC, close_enough);
        } else {
            BOOST_REQUIRE_CLOSE(out_real*NC, buffer[i][0], close_enough);
        }
        if (fabs(buffer[i][1]) < close_enough) {
            BOOST_REQUIRE_SMALL(out_imag*NC, close_enough);
        } else {
            BOOST_REQUIRE_CLOSE(out_imag*NC, buffer[i][1], close_enough);
        }
    }
}

// Helper function that kicks the tires of a 1D c2c transform
template<class ComplexMultiArray1, class ComplexMultiArray2>
void symmetry_1D_complex_backward(ComplexMultiArray1 &in, ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    typedef typename pencilfft::internal::transform_traits<
        typename ComplexMultiArray2::element>::real_type real_type;
    const int NC = in.shape()[0];
    const int NR = out.shape()[0];
    const real_type close_enough
        = std::numeric_limits<real_type>::epsilon()*10*NC*NC;

    // Load a conjugate-symmetric function into the input array...
    // ...with known frequency content and constant offset 17
    fill_with_NaN(out);
    suzerain::multi_array::fill(in, 0);
    suzerain::complex::assign_complex(in[0], 17.0, 0.0);
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) {
        const real_type val = i/2.0;
        suzerain::complex::assign_complex(in[i],    val, -val);
        suzerain::complex::assign_complex(in[NC-i], val,  val);
    }
    // TODO Test behavior of highest even half mode

    pencilfft::backward_c2c(0, in, out);

    // Output should be a real-valued function in physical space...
    for (int i = 0; i < NR; ++i) {
        real_type a_real, a_imag;
        suzerain::complex::assign_components(a_real, a_imag, out[i]);
        BOOST_REQUIRE_SMALL(a_imag, close_enough);

        // ...with known physical space values
        const real_type xi = i*2*M_PI/NR;
        real_type expected_real = 17;
        for (int j = 1; j < (std::min(NC,NR)+1)/2; ++j) {
            expected_real += j*sin(j*xi) + j*cos(j*xi);
        }
        BOOST_REQUIRE_CLOSE(a_real, expected_real, sqrt(close_enough));
    }

    // Load a real-symmetric function into the input array...
    // ...with known frequency content and constant offset 17
    fill_with_NaN(out);
    suzerain::multi_array::fill(in, 0);
    suzerain::complex::assign_complex(in[0], 0.0, 17.0);
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) {
        const real_type val = i/2.0;
        suzerain::complex::assign_complex(in[i],     val, val);
        suzerain::complex::assign_complex(in[NC-i], -val, val);
    }
    // TODO Test behavior of highest even half mode

    pencilfft::backward_c2c(0, in, out);

    // Output should be an imaginary-valued function in physical space...
    for (int i = 0; i < NR; ++i) {
        real_type a_real, a_imag;
        suzerain::complex::assign_components(a_real, a_imag, out[i]);
        BOOST_REQUIRE_SMALL(a_real, close_enough);

        // ...with known physical space values
        const real_type xi = i*2*M_PI/NR;
        real_type expected_imag = 17;
        for (int j = 1; j < (std::min(NC,NR)+1)/2; ++j) {
            expected_imag += j*sin(j*xi) + j*cos(j*xi);
        }
        BOOST_REQUIRE_CLOSE(a_imag, expected_imag, sqrt(close_enough));
    }
}

// Helper function that kicks the tires of a 1D c2r transform
template<class ComplexMultiArray, class RealMultiArray>
void check_1D_backward_c2r(ComplexMultiArray &in, RealMultiArray &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray::dimensionality == 1);
    BOOST_STATIC_ASSERT(RealMultiArray::dimensionality == 1);
    typedef typename pencilfft::internal::transform_traits<
        typename ComplexMultiArray::element>::real_type real_type;
    const int NC = in.shape()[0];
    const int NR = out.shape()[0];
    const real_type close_enough
        = std::numeric_limits<real_type>::epsilon()*10*NR*NR;

    // Load a function into the input array...
    // ...with known frequency content and constant offset 17
    fill_with_NaN(out);
    suzerain::multi_array::fill(in,0);
    suzerain::complex::assign_complex(in[0], 17.0, 0.0);
    for (int i = 1; i < (std::min(NC,NR)+1)/2; ++i) {
        const real_type val = i/2.0;
        suzerain::complex::assign_complex(in[i], val, -val);
    }
    // TODO Test behavior of highest even half mode

    pencilfft::backward_c2r(0, in, out);

    // Output should be the expected function
    for (int i = 0; i < NR; ++i) {
        real_type a_real = out[i];

        const real_type xi = i*2*M_PI/NR;
        real_type expected_real = 17;
        for (int j = 1; j < (std::min(NC,NR)+1)/2; ++j) {
            expected_real += j*sin(j*xi) + j*cos(j*xi);
        }
        BOOST_REQUIRE_CLOSE(a_real, expected_real, sqrt(close_enough));
    }
}

// Compare our results with FFTW's directly computed result
template<class ComplexMultiArray1, class ComplexMultiArray2>
void compare_1D_complex_backward(ComplexMultiArray1 &in,
                                 ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    const int NC = in.shape()[0];
    const int NR = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*1e2*NC*NC;

    // Plan before loading in the data since planning overwrites in
    suzerain::shared_array<fftw_complex> buffer(
        static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex)*NC)),
        std::ptr_fun(fftw_free));
    BOOST_REQUIRE(buffer.get() != NULL);
    typename ComplexMultiArray1::element * const in_data
        = in.origin() + std::inner_product(
            in.index_bases(), in.index_bases()+1, in.strides(), 0);
    suzerain::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan(
            fftw_plan_many_dft(1,
                               &NC,
                               1,
                               reinterpret_cast<fftw_complex*>(in_data),
                               NULL,
                               in.strides()[0],
                               NC,
                               buffer.get(),
                               NULL,
                               1,
                               NC,
                               FFTW_BACKWARD,
                               FFTW_PRESERVE_INPUT),
            std::ptr_fun(fftw_destroy_plan));
    BOOST_REQUIRE(plan.get() != NULL);

    // Load a complex-valued function into the input array
    fill_with_NaN(in);
    fill_with_NaN(out);
    for (int i = 0; i < NC; ++i) {
        suzerain::complex::assign_complex(in[i], i, i);
    }

    // Transform input FFTW directly and also our wrapper
    // Important for FFTW to be first to handle in === out

    // Ramp up amplitudes for FFTW's backwards transform
    for (int i = 0; i < NC; ++i) {
        suzerain::complex::assign_complex_scaled(in[i], in[i], NC);
    }
    fftw_execute(plan.get());
    // Ramp down amplitudes for our backwards transform
    for (int i = 0; i < NC; ++i) {
        suzerain::complex::assign_complex_scaled(in[i], in[i], 1.0/NC);
    }
    pencilfft::backward_c2c(0, in, out);

    // Ensure we got exactly the same result after normalization
    for (int i = 0; i < NR; ++i) {
        double out_real, out_imag;
        suzerain::complex::assign_components(
                out_real, out_imag, out[i]);
        // FFTW with a stride gives a different result than with stride 1
        // BOOST_REQUIRE_EQUAL would be nice, but it fails here
        BOOST_REQUIRE_CLOSE(out_real, buffer[i][0] / NR, close_enough);
        BOOST_REQUIRE_CLOSE(out_imag, buffer[i][1] / NR, close_enough);
    }
}

template<class ComplexMultiArray1, class ComplexMultiArray2>
void differentiate_on_forward_1D_c2c(ComplexMultiArray1 &in,
                                     ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    typedef typename pencilfft::internal::transform_traits<
        typename ComplexMultiArray1::element>::real_type real_type;
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*NR*NR*NR;

    const double length[2] = { 2.0 * M_PI, 10.0 };
    for (std::size_t l = 0; l < sizeof(length)/sizeof(length[0]); ++l) {
        for (int derivative = 0; derivative < 8; ++derivative) {
            // Load a complex function into the input array
            fill_with_NaN(in);
            fill_with_NaN(out);
            const periodic_function<real_type,int> pf(
                    NR, (std::min(NR,NC)+1)/2, M_PI/3.0, length[l]);
            for (int i = 0; i < NR; ++i) {
                const double val =  pf.physical(i, 0);
                suzerain::complex::assign_complex(in[i], val, -val);
            }

            // Forward transform and differentiate
            pencilfft::forward_c2c(0, in, out, length[l], derivative);

            // Backwards transform without differentiating
            pencilfft::backward_c2c(0, out, in, length[l], 0);

            // Ensure we see what we expect
            for (int i = 0; i < NR; ++i) {
                const double expected_val = pf.physical(i, derivative);
                double real, imag;
                suzerain::complex::assign_components(real, imag, in[i]);
                if (fabs(expected_val) < close_enough) {
                    BOOST_REQUIRE_SMALL(real, close_enough);
                    BOOST_REQUIRE_SMALL(imag, close_enough);
                } else {
                    BOOST_REQUIRE_CLOSE(
                            expected_val, real, sqrt(close_enough));
                    BOOST_REQUIRE_CLOSE(
                           -expected_val, imag, sqrt(close_enough));
                }
            }
        }
    }
}

template<class ComplexMultiArray1, class ComplexMultiArray2>
void differentiate_on_backward_1D_c2c(ComplexMultiArray1 &in,
                                      ComplexMultiArray2 &out)
{
    BOOST_STATIC_ASSERT(ComplexMultiArray1::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray2::dimensionality == 1);
    typedef typename pencilfft::internal::transform_traits<
        typename ComplexMultiArray1::element>::real_type real_type;
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*10*NR*NR*NR;

    const double length[2] = { 2.0 * M_PI, 10.0 };
    for (std::size_t l = 0; l < sizeof(length)/sizeof(length[0]); ++l) {
        for (int derivative = 0; derivative < 8; ++derivative) {
            // Load a complex function into the input array
            fill_with_NaN(in);
            fill_with_NaN(out);
            const periodic_function<real_type,int> pf(
                    NR, (std::min(NR,NC)+1)/2, M_PI/3.0, length[l]);
            for (int i = 0; i < NR; ++i) {
                const double val = pf.physical(i, 0);
                suzerain::complex::assign_complex(in[i], val, -val);
            }

            // Forward transform without differentiating
            pencilfft::forward_c2c(0, in, out, length[l], 0);

            // Backwards transform and differentiate
            pencilfft::backward_c2c(0, out, in, length[l], derivative);

            // Ensure we see what we expect as the derivative
            for (int i = 0; i < NR; ++i) {
                const double expected_val = pf.physical(i, derivative);
                double real, imag;
                suzerain::complex::assign_components(real, imag, in[i]);
                if (fabs(expected_val) < close_enough) {
                    BOOST_REQUIRE_SMALL(real, close_enough);
                    BOOST_REQUIRE_SMALL(imag, close_enough);
                } else {
                    BOOST_REQUIRE_CLOSE(
                             expected_val, real, sqrt(close_enough));
                    BOOST_REQUIRE_CLOSE(
                            -expected_val, imag, sqrt(close_enough));
                }
            }
        }
    }
}

template<class RealMultiArray, class ComplexMultiArray>
void differentiate_on_forward_1D_r2c(RealMultiArray &in,
                                     ComplexMultiArray &out)
{
    BOOST_STATIC_ASSERT(RealMultiArray::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray::dimensionality == 1);
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*100*NR*NR*NR;
    typedef typename RealMultiArray::element real_type;

    const double length[2] = { 2.0 * M_PI, 10.0 };
    for (std::size_t l = 0; l < sizeof(length)/sizeof(length[0]); ++l) {
        for (int derivative = 0; derivative < 8; ++derivative) {
            // Load a real-valued function into the input array
            fill_with_NaN(in);
            fill_with_NaN(out);
            const periodic_function<real_type,int> pf(
                    NR, (std::min(NR,NC)+1)/2, M_PI/3.0, length[l]);
            for (int i = 0; i < NR; ++i) {
                in[i] = pf.physical(i, 0);
            }

            // Forward transform and differentiate
            pencilfft::forward_r2c(0, in, out, length[l], derivative);

            // Backwards transform without differentiating
            pencilfft::backward_c2r(0, out, in, length[l], 0);

            // Ensure we see what we expect
            for (int i = 0; i < NR; ++i) {
                const double expected_val = pf.physical(i, derivative);
                const double observed_val = in[i];
                if (fabs(expected_val) < close_enough) {
                    BOOST_REQUIRE_SMALL(observed_val, close_enough);
                } else {
                    BOOST_REQUIRE_CLOSE(
                            expected_val, observed_val, sqrt(close_enough));
                }
            }
        }
    }
}

template<class RealMultiArray, class ComplexMultiArray>
void differentiate_on_backward_1D_c2r(RealMultiArray &in,
                                      ComplexMultiArray &out)
{
    BOOST_STATIC_ASSERT(RealMultiArray::dimensionality == 1);
    BOOST_STATIC_ASSERT(ComplexMultiArray::dimensionality == 1);
    typedef typename RealMultiArray::element real_type;
    const int NR = in.shape()[0];
    const int NC = out.shape()[0];
    const double close_enough
        = std::numeric_limits<double>::epsilon()*100*NR*NR*NR;

    const double length[2] = { 2.0 * M_PI, 10.0 };
    for (std::size_t l = 0; l < sizeof(length)/sizeof(length[0]); ++l) {
        for (int derivative = 0; derivative < 8; ++derivative) {
            // Load a real-valued function into the input array
            fill_with_NaN(in);
            fill_with_NaN(out);
            const periodic_function<real_type,int> pf(
                    NR, (std::min(NR,NC)+1)/2, M_PI/3.0, length[l]);
            for (int i = 0; i < NR; ++i) {
                in[i] = pf.physical(i, 0);
            }

            // Forward transform without differentiating
            pencilfft::forward_r2c(0, in, out, length[l], 0);

            // Backward transform and differentiate
            pencilfft::backward_c2r(0, out, in, length[l], derivative);

            // Ensure we see what we expect
            for (int i = 0; i < NR; ++i) {
                const double expected_val = pf.physical(i, derivative);
                const double observed_val = in[i];
                if (fabs(expected_val) < close_enough) {
                    BOOST_REQUIRE_SMALL(observed_val, close_enough);
                } else {
                    BOOST_REQUIRE_CLOSE(
                            expected_val, observed_val, sqrt(close_enough));
                }
            }
        }
    }
}


/* powers of 2; simple composites of form 2*prime; prime numbers; 1 */
#define TRANSFORM_1D_SIZE_SEQ \
        (2)(4)(8)(16)(32)(64) \
        (6)(10)(14)(22)(26)(34) \
        (3)(5)(7)(9)(11)(13)(17)(19)(23)(29) \
        (1)

BOOST_AUTO_TEST_SUITE( test_1d_out_of_place );
#define TEST_1D_OUT_OF_PLACE(r, data, elem) \
        BOOST_AUTO_TEST_CASE( BOOST_PP_CAT(test_1d_out_of_place_,elem) ) \
        { test_1d_out_of_place(elem); }
void test_1d_out_of_place(const int N)
{
    // C2C: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,1> array_type;
        array_type in(boost::extents[N]), out(boost::extents[N]);

        // No dealiasing in effect: size in == size out
        symmetry_1D_complex_forward(in, out);
        compare_1D_complex_forward(in, out);
        symmetry_1D_complex_backward(in, out);
        compare_1D_complex_backward(in, out);
        differentiate_on_forward_1D_c2c(in, out);
        differentiate_on_backward_1D_c2c(in, out);
    }

    // C2C: Test multi_array_ref using fftw_complex
    {
        suzerain::scoped_array<fftw_complex> in_data(new fftw_complex[N]);
        suzerain::scoped_array<fftw_complex> out_data(new fftw_complex[N]);
        typedef boost::multi_array_ref<fftw_complex, 1> array_ref_type;
        array_ref_type in(in_data.get(), boost::extents[N]);
        array_ref_type out(out_data.get(), boost::extents[N]);

        // No dealiasing in effect: size in == size out
        symmetry_1D_complex_forward(in, out);
        compare_1D_complex_forward(in, out);
        symmetry_1D_complex_backward(in, out);
        compare_1D_complex_backward(in, out);
        differentiate_on_forward_1D_c2c(in, out);
        differentiate_on_backward_1D_c2c(in, out);
    }

    // R2C: Test multi_array using std::complex
    {
        boost::multi_array<double,1>               in(boost::extents[N]);
        boost::multi_array<std::complex<double>,1> out(boost::extents[N/2+1]);

        // No dealiasing in effect
        check_1D_forward_r2c(in, out);
    }

    // R2C: Test multi_array using fftw_complex
    {
        boost::multi_array<double,1>      in(boost::extents[N]);
        suzerain::scoped_array<fftw_complex> out_data(new fftw_complex[N/2+1]);
        boost::multi_array_ref<fftw_complex, 1> out(
                out_data.get(), boost::extents[N/2+1]);

        // No dealiasing in effect
        check_1D_forward_r2c(in, out);
    }

    // C2R: Test multi_array using std::complex
    {
        boost::multi_array<std::complex<double>,1> in(boost::extents[N/2+1]);
        boost::multi_array<double,1>               out(boost::extents[N]);

        // No dealiasing in effect
        check_1D_backward_c2r(in, out);
    }

    // C2R: Test multi_array using fftw_complex
    {
        suzerain::scoped_array<fftw_complex> in_data(new fftw_complex[N/2+1]);
        boost::multi_array_ref<fftw_complex, 1> in(
                in_data.get(), boost::extents[N/2+1]);
        boost::multi_array<double,1>      out(boost::extents[N]);

        // No dealiasing in effect
        check_1D_backward_c2r(in, out);
    }
}
BOOST_PP_SEQ_FOR_EACH(TEST_1D_OUT_OF_PLACE,_,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE( test_1d_out_of_place_one_reversed );
#define TEST_1D_OUT_OF_PLACE_ONE_REVERSED(r, data, elem) \
        BOOST_AUTO_TEST_CASE( \
          BOOST_PP_CAT(test_1d_out_of_place_one_reversed_,elem) ) \
        { test_1d_out_of_place_one_reversed(elem); }
void test_1d_out_of_place_one_reversed(const int N)
{
    const int ndim = 1;

    // C2C: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,ndim> array_type;
        typedef boost::general_storage_order<ndim> storage;
        array_type::size_type ordering[ndim] = { 0 };
        const bool ascending[ndim] = { false };

        array_type in(boost::extents[N], storage(ordering, ascending));
        array_type out(boost::extents[N]);

        symmetry_1D_complex_forward(in, out);
        compare_1D_complex_forward(in, out);
        symmetry_1D_complex_backward(in, out);
        compare_1D_complex_backward(in, out);
        differentiate_on_forward_1D_c2c(in, out);
        differentiate_on_backward_1D_c2c(in, out);
    }

    // R2C: Test multi_array using std::complex when real storage reversed
    {
        typedef boost::multi_array<std::complex<double>,ndim> complex_array_type;
        typedef boost::multi_array<double,ndim>               real_array_type;

        typedef boost::general_storage_order<ndim> storage;
        real_array_type::size_type ordering[ndim] = { 0 };
        const bool ascending[ndim] = { false };

        real_array_type in(boost::extents[N], storage(ordering, ascending));
        complex_array_type out(boost::extents[N/2+1]);

        // No dealiasing in effect
        check_1D_forward_r2c(in, out);
    }

    // R2C: Test multi_array using std::complex when complex storage reversed
    {
        typedef boost::multi_array<std::complex<double>,ndim> complex_array_type;
        typedef boost::multi_array<double,ndim>               real_array_type;

        typedef boost::general_storage_order<ndim> storage;
        complex_array_type::size_type ordering[ndim] = { 0 };
        const bool ascending[ndim] = { false };

        real_array_type in(boost::extents[N]);
        complex_array_type out(boost::extents[N/2+1],
                               storage(ordering, ascending));

        // No dealiasing in effect
        check_1D_forward_r2c(in, out);
    }

    // C2R: Test multi_array using std::complex when real storage reversed
    {
        typedef boost::multi_array<std::complex<double>,ndim> complex_array_type;
        typedef boost::multi_array<double,ndim>               real_array_type;

        typedef boost::general_storage_order<ndim> storage;
        real_array_type::size_type ordering[ndim] = { 0 };
        const bool ascending[ndim] = { false };

        complex_array_type in(boost::extents[N/2+1]);
        real_array_type out(boost::extents[N], storage(ordering, ascending));

        // No dealiasing in effect
        check_1D_backward_c2r(in, out);
    }

    // C2R: Test multi_array using std::complex when complex storage reversed
    {
        typedef boost::multi_array<std::complex<double>,ndim> complex_array_type;
        typedef boost::multi_array<double,ndim>               real_array_type;

        typedef boost::general_storage_order<ndim> storage;
        complex_array_type::size_type ordering[ndim] = { 0 };
        const bool ascending[ndim] = { false };

        complex_array_type in(boost::extents[N/2+1],
                               storage(ordering, ascending));
        real_array_type out(boost::extents[N]);

        // No dealiasing in effect
        check_1D_backward_c2r(in, out);
    }
}
BOOST_PP_SEQ_FOR_EACH(TEST_1D_OUT_OF_PLACE_ONE_REVERSED,\
    _,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE( test_1d_out_of_place_two_reversed );
#define TEST_1D_OUT_OF_PLACE_TWO_REVERSED(r, data, elem) \
        BOOST_AUTO_TEST_CASE( \
          BOOST_PP_CAT(test_1d_out_of_place_two_reversed_,elem) ) \
        { test_1d_out_of_place_two_reversed(elem); }
void test_1d_out_of_place_two_reversed(const int N)
{
    const int ndim = 1;

    // C2C: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,ndim> array_type;
        typedef boost::general_storage_order<ndim> storage;
        array_type::size_type ordering[ndim] = { 0 };
        const bool ascending[ndim] = { false };

        array_type in(boost::extents[N], storage(ordering, ascending));
        array_type out(boost::extents[N], storage(ordering, ascending));

        symmetry_1D_complex_forward(in, out);
        compare_1D_complex_forward(in, out);
        symmetry_1D_complex_backward(in, out);
        compare_1D_complex_backward(in, out);
        differentiate_on_forward_1D_c2c(in, out);
        differentiate_on_backward_1D_c2c(in, out);
    }

    // R2C: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,ndim> complex_array_type;
        typedef boost::multi_array<double,ndim>               real_array_type;

        typedef boost::general_storage_order<ndim> storage;
        real_array_type::size_type ordering[ndim] = { 0 };
        const bool ascending[ndim] = { false };

        real_array_type in(boost::extents[N],
                           storage(ordering, ascending));
        complex_array_type out(boost::extents[N/2+1],
                               storage(ordering, ascending));

        // No dealiasing in effect
        check_1D_forward_r2c(in, out);
    }

    // C2R: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,ndim> complex_array_type;
        typedef boost::multi_array<double,ndim>               real_array_type;

        typedef boost::general_storage_order<ndim> storage;
        real_array_type::size_type ordering[ndim] = { 0 };
        const bool ascending[ndim] = { false };

        complex_array_type in(boost::extents[N/2+1],
                              storage(ordering, ascending));
        real_array_type out(boost::extents[N],
                           storage(ordering, ascending));

        // No dealiasing in effect
        check_1D_backward_c2r(in, out);
    }
}
BOOST_PP_SEQ_FOR_EACH(TEST_1D_OUT_OF_PLACE_TWO_REVERSED,\
    _,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE( test_1d_in_place );
#define TEST_1D_IN_PLACE(r, data, elem) \
        BOOST_AUTO_TEST_CASE( BOOST_PP_CAT(test_1d_in_place_,elem) ) \
        { test_1d_in_place(elem); }
void test_1d_in_place(const int N)
{
    // C2C: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,1> array_type;
        array_type both(boost::extents[N]);

        // No dealiasing in effect for in place transform: NR == NC
        symmetry_1D_complex_forward(both, both);
        compare_1D_complex_forward(both, both);
        symmetry_1D_complex_backward(both, both);
        compare_1D_complex_forward(both, both);
        differentiate_on_forward_1D_c2c(both, both);
        differentiate_on_backward_1D_c2c(both, both);
    }

    // R2C: Test multi_array_ref using std::complex
    // Test performed in place; API requires both real and complex views
    {
        typedef std::complex<double> complex;
        suzerain::scoped_array<complex> raw(new complex[N/2+1]);

        typedef boost::multi_array_ref<complex::value_type,1> real_array_type;
        real_array_type in(reinterpret_cast<complex::value_type *>(raw.get()),
                           boost::extents[N]);

        typedef boost::multi_array_ref<complex,1> complex_array_type;
        complex_array_type out(raw.get(), boost::extents[N/2+1]);

        // No dealiasing in effect
        check_1D_forward_r2c(in, out);
    }

    // C2R: Test multi_array_ref using std::complex
    // Test performed in place; API requires both real and complex views
    {
        typedef std::complex<double> complex;
        suzerain::scoped_array<complex> raw(new complex[N/2+1]);

        typedef boost::multi_array_ref<complex,1> complex_array_type;
        complex_array_type in(raw.get(), boost::extents[N/2+1]);

        typedef boost::multi_array_ref<complex::value_type,1> real_array_type;
        real_array_type out(reinterpret_cast<complex::value_type *>(raw.get()),
                            boost::extents[N]);

        // No dealiasing in effect
        check_1D_backward_c2r(in, out);
    }
}
BOOST_PP_SEQ_FOR_EACH(TEST_1D_IN_PLACE,_,TRANSFORM_1D_SIZE_SEQ);
BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE( test_1d_out_of_place_dealiased );
#define TEST_1D_OUT_OF_PLACE_DEALIASED(r, product) \
        BOOST_AUTO_TEST_CASE( \
          BOOST_PP_CAT(test_1d_out_of_place_dealiased_, \
            BOOST_PP_CAT(BOOST_PP_SEQ_ELEM(0,product), \
              BOOST_PP_CAT(_, \
                BOOST_PP_CAT(BOOST_PP_SEQ_ELEM(1,product), \
                  BOOST_PP_CAT(_,BOOST_PP_SEQ_ELEM(2,product))))))) \
        {\
            BOOST_PP_CAT(test_1d_out_of_place_dealiased_, \
                         BOOST_PP_SEQ_ELEM(0,product))( \
                BOOST_PP_SEQ_ELEM(1,product), \
                BOOST_PP_SEQ_ELEM(2,product) ); \
        }
#define TEST_1D_OUT_OF_PLACE_DEALIASED_LESS_ONLY(r,product) \
    BOOST_PP_IIF(BOOST_PP_LESS(BOOST_PP_SEQ_ELEM(1,product),\
                BOOST_PP_SEQ_ELEM(2,product)), \
                TEST_1D_OUT_OF_PLACE_DEALIASED(r,product), \
                BOOST_PP_EMPTY());
#define TEST_1D_OUT_OF_PLACE_DEALIASED_GREATER_ONLY(r,product) \
    BOOST_PP_IIF(BOOST_PP_GREATER(BOOST_PP_SEQ_ELEM(1,product),\
                BOOST_PP_SEQ_ELEM(2,product)), \
                TEST_1D_OUT_OF_PLACE_DEALIASED(r,product), \
                BOOST_PP_EMPTY());
void test_1d_out_of_place_dealiased_forward(const int NR, const int NC)
{
    // C2C: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,1> array_type;
        array_type in(boost::extents[NR]), out(boost::extents[NC]);
        symmetry_1D_complex_forward(in, out); // Dealiasing in effect
        differentiate_on_forward_1D_c2c(in, out);
    }

    // R2C and C2R: Test multi_array using std::complex
    {
        typedef boost::multi_array<double,1>               real_array_type;
        typedef boost::multi_array<std::complex<double>,1> complex_array_type;
        real_array_type in(boost::extents[NR]);
        complex_array_type out(boost::extents[NC]);
        differentiate_on_forward_1D_r2c(in, out);
    }
}
void test_1d_out_of_place_dealiased_backward(const int NC, const int NR)
{
    // C2C: Test multi_array using std::complex
    {
        typedef boost::multi_array<std::complex<double>,1> array_type;
        array_type in(boost::extents[NC]), out(boost::extents[NR]);
        symmetry_1D_complex_backward(in, out); // Dealiasing in effect
        differentiate_on_backward_1D_c2c(out, in); // NB reversed
    }

    // C2R and C2R: Test multi_array using std::complex
    {
        typedef boost::multi_array<double,1>               real_array_type;
        typedef boost::multi_array<std::complex<double>,1> complex_array_type;
        real_array_type in(boost::extents[NR]);
        complex_array_type out(boost::extents[NC]);
        differentiate_on_backward_1D_c2r(in, out);
    }
}
BOOST_PP_SEQ_FOR_EACH_PRODUCT(\
        TEST_1D_OUT_OF_PLACE_DEALIASED_LESS_ONLY, \
        ((forward))(TRANSFORM_1D_SIZE_SEQ)(TRANSFORM_1D_SIZE_SEQ) );
BOOST_PP_SEQ_FOR_EACH_PRODUCT(\
        TEST_1D_OUT_OF_PLACE_DEALIASED_GREATER_ONLY, \
        ((forward))(TRANSFORM_1D_SIZE_SEQ)(TRANSFORM_1D_SIZE_SEQ) );
BOOST_PP_SEQ_FOR_EACH_PRODUCT(\
        TEST_1D_OUT_OF_PLACE_DEALIASED_LESS_ONLY, \
        ((backward))(TRANSFORM_1D_SIZE_SEQ)(TRANSFORM_1D_SIZE_SEQ) );
BOOST_PP_SEQ_FOR_EACH_PRODUCT(\
        TEST_1D_OUT_OF_PLACE_DEALIASED_GREATER_ONLY, \
        ((backward))(TRANSFORM_1D_SIZE_SEQ)(TRANSFORM_1D_SIZE_SEQ) );
BOOST_AUTO_TEST_SUITE_END();
