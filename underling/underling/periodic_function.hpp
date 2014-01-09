//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$

#ifndef UNDERLING_PERIODIC_FUNCTION_HPP
#define UNDERLING_PERIODIC_FUNCTION_HPP

#include <cassert>
#include <complex>
#include <cmath>

/** @file
 * Provides the test-ready periodic_function template.
 */

namespace underling {

/**
 * Provides a periodic real-valued function useful for testing FFT behavior.
 * Wave space coefficients are stored "in-order".  Derivatives of the function
 * can be obtained.
 * */
template<typename FPT = double, typename Integer = int>
class periodic_function {
public:

    /**
     * Produce a periodic, univariate, real-valued signal with known frequency
     * content on domain of supplied length.
     *
     * @param N Number of points in the physical domain
     * @param max_mode_exclusive Exclusive upper bound on the signal's
     *        frequency content expressed as a wavenumber index.  Setting
     *        <tt>max_mode_exclusive <= 0</tt> indicates that the maximum
     *        supportable frequency given \c N should be used.
     * @param shift Phase shift in the signal
     * @param length Domain size over which the signal is periodic
     * @param constant The constant content of the signal and an amplitude
     *                 factor used to scale all modes.
     */
     explicit periodic_function(const Integer N,
                                const Integer max_mode_exclusive = 0,
                                const FPT shift = M_PI/3,
                                const FPT length = 2*M_PI,
                                const FPT constant = 17)
        : N(N),
          max_mode_exclusive(
                    max_mode_exclusive > 0 ? max_mode_exclusive : N/2+1),
          shift(shift),
          length(length),
          constant(constant)
        {
            assert(N > 0);
            assert(max_mode_exclusive <= (N/2+1));
        }

    /**
     * Retrieve the real-valued signal amplitude at the given location.
     *
     * @param x The desired grid location.
     *          Should be within <tt>[0,length)</tt>.
     * @param derivative Desired derivative of the signal with zero
     *        indicating the signal itself.
     *
     * @param Returns the requested value.
     */
    FPT physical_evaluate(const FPT x, const Integer derivative = 0) const;

    /**
     * Retrieve the real-valued signal amplitude at the given gridpoint.
     *
     * @param i Zero-indexed value of the desired grid location.
     *        Must be within <tt>[0,N)</tt>
     * @param derivative Desired derivative of the signal with zero
     *        indicating the signal itself.
     *
     * @param Returns the requested value.
     */
    FPT physical(const Integer i, const Integer derivative = 0) const;

    /**
     * Retrieve the requested complex-valued signal modes in wave space.
     *
     * @param i Zero-indexed mode number, where zero is the constant
     *        mode.
     * @param derivative Desired derivative of the signal with zero
     *        indicating the signal itself.
     *
     * @param Returns the requested value.
     */
    typename std::complex<FPT> wave(const Integer i,
                                    const Integer derivative = 0) const;

    const Integer N;
    const Integer max_mode_exclusive;
    const FPT shift;
    const FPT length;
    const FPT constant;
};

template<typename FPT, typename Integer>
FPT periodic_function<FPT,Integer>::physical_evaluate(
        const FPT x,
        const Integer derivative) const
{
    assert(0 <= (derivative % 4) && (derivative % 4) <= 3);

    FPT retval = (max_mode_exclusive > 0 && derivative == 0 )
        ? constant : 0;
    for (Integer j = 1; j < max_mode_exclusive; ++j) {
        switch (derivative % 4) {
            case 0:
                retval +=   j * pow(j*(2.0*M_PI/length), derivative)
                          * constant
                          * sin(j*(2.0*M_PI/length)*x + shift);
                break;
            case 1:
                retval +=   j * pow(j*(2.0*M_PI/length), derivative)
                          * constant
                          * cos(j*(2.0*M_PI/length)*x + shift);
                break;
            case 2:
                retval -=   j * pow(j*(2.0*M_PI/length), derivative)
                          * constant
                          * sin(j*(2.0*M_PI/length)*x + shift);
                break;
            case 3:
                retval -=   j * pow(j*(2.0*M_PI/length), derivative)
                          * constant
                          * cos(j*(2.0*M_PI/length)*x + shift);
                break;
        }
    }
    return retval;
}

template<typename FPT, typename Integer>
FPT periodic_function<FPT,Integer>::physical(
        const Integer i,
        const Integer derivative) const
{
    assert(0 <= i);
    assert(i <  N);

    const FPT xi = i*length/N;
    return physical_evaluate(xi, derivative);
}

template<typename FPT, typename Integer>
typename std::complex<FPT> periodic_function<FPT,Integer>::wave(
        const Integer i,
        const Integer derivative) const
{
    assert(0 <= i && i < N);

    typedef typename std::complex<FPT> complex_type;
    complex_type retval;
    if (i == 0) {
        if (i < max_mode_exclusive) {
            retval = complex_type(1, 0);
        }
    } else if (i < (N/2)) {
        if (i < max_mode_exclusive) {
            retval = complex_type(   (i)*sin(shift)/2, -  (i)*cos(shift)/2 );
        }
    } else if (i == (N/2) &&  (N%2)) {
        if (i < max_mode_exclusive) {
            retval = complex_type(   (i)*sin(shift)/2, -  (i)*cos(shift)/2 );
        }
    } else if (i == (N/2) && !(N%2)) { // Highest half mode on even grid
        if (i < max_mode_exclusive) {
            retval = complex_type(   (i)*sin(shift), 0    );
        }
    } else {
        if ((N-i) < max_mode_exclusive) {
            retval = complex_type( (N-i)*sin(shift)/2, +(N-i)*cos(shift)/2 );
        }
    }
    retval *= constant;

    for (int j = 0; j < derivative; ++j) {
        retval *= complex_type(0, 2*M_PI/length);
    }

    return retval;
}

} // namespace underling

#endif // UNDERLING_PERIODIC_FUNCTION_HPP
