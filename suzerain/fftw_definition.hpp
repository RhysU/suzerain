//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// fftw_definition.hpp: A definition_base to configure FFTW from the CLI
// $Id$

#ifndef SUZERAIN_FFTW_DEFINITION_HPP
#define SUZERAIN_FFTW_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/definition_base.hpp>
#include <suzerain/fftw.hpp>

/** @file
 * A \ref definition_base to configure FFTW operation from the CLI.
 */

namespace suzerain {

/** Holds FFTW-usage parameters, e.g. the planning rigor. */
class FFTWDefinition : public definition_base
{
public:

    /** Default constructor */
    explicit FFTWDefinition(const fftw::rigor rigor_fft = fftw::estimate,
                            const fftw::rigor rigor_mpi = fftw::estimate);

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

} // namespace suzerain

#endif // SUZERAIN_FFTW_DEFINITION_HPP
