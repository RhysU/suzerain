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

#ifndef SUZERAIN_SUPPORT_FFTW_DEFINITION_HPP
#define SUZERAIN_SUPPORT_FFTW_DEFINITION_HPP

/** @file
 * Declarations to configure FFTW operation from the CLI.
 */

#include <suzerain/common.hpp>
#include <suzerain/fftw.hpp>
#include <suzerain/support/definition_base.hpp>

namespace suzerain {

namespace support {

/** Holds FFTW-usage parameters, e.g. the planning rigor. */
class definition_fftw
    : public virtual definition_base
{
public:

    /** Default constructor */
    explicit definition_fftw(
            const fftw::rigor rigor_fft = fftw::estimate,
            const fftw::rigor rigor_mpi = fftw::estimate);

    /** @copydoc support::definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

    /**
     * The FFTW rigor flag intended for FFT planning.
     * @see rigor() for more details.
     */
    fftw::rigor rigor_fft;

    /**
     * The FFTW rigor flag intended for FFTW MPI communication planning.
     * @see rigor() for more details.
     */
    fftw::rigor rigor_mpi;

    /**
     * The file location used to accumulate FFTW wisdom.
     */
    std::string plan_wisdom;

    /**
     * The planning time limit for use with <tt>fftw_set_timelimit</tt>.
     * Will be a concrete time limit if one was specified or FFTW_NO_TIMELIMIT.
     */
    double plan_timelimit;

private:

    /** Normalizes and stores a value in <tt>this->rigor_fft_</tt>. */
    void normalize_rigor_fft(std::string input);

    /** Normalizes and stores a value in <tt>this->rigor_mpi_</tt>. */
    void normalize_rigor_mpi(std::string input);
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_FFTW_DEFINITION_HPP
