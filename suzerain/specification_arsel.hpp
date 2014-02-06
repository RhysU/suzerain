//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
// Copyright (C) 2014 The PECOS Development Team
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

#ifndef SUZERAIN_SPECIFICATION_ARSEL_HPP
#define SUZERAIN_SPECIFICATION_ARSEL_HPP

/** @file
 * Provides \ref specification_arsel.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Encapsulates settings for invoking \ref arsel.
 */
class specification_arsel
{

public:

    /** Construct an instance with the recommended default values. */
    specification_arsel(const bool         absrho    = false,
                        const std::string& criterion = "CIC",
                        const std::size_t  maxorder  = 512,
                        const std::size_t  minorder  = 0,
                        const real_t       wlenT0    = 7);

    /** Use absolute autocorrelation when integrating for T0. */
    bool absrho;

    /**
     * The model order will be selected using the specified criterion.
     *
     * Criteria are specified using the following abbreviations:
     * \li \c AIC  - Akaike information criterion
     * \li \c AICC - asymptotically-corrected Akaike information criterion
     * \li \c BIC  - consistent criterion BIC
     * \li \c CIC  - combined information criterion
     * \li \c FIC  - finite information criterion
     * \li \c FSIC - finite sample information criterion
     * \li \c GIC  - generalized information criterion
     * \li \c MCC  - minimally consistent criterion
     *
     * \throws std::invalid_argument on unknown \c criterion.
     * \returns The newly set criterion abbreviation.
     */
    const std::string& criterion(const std::string& abbrev);

    /** Retrieve the currently-specified criterion. */
    const std::string& criterion() const;

    /** Consider only models of at most order AR(<code>maxorder</code>). */
    std::size_t maxorder;

    /** Consider only models of at least order AR(<code>minorder</code>). */
    std::size_t minorder;

    /** Subtract the mean from the signal iif true.  Always.  */
    static const bool submean = true;

    /** Integrate for T0 until WLEN times the input length. */
    real_t wlenT0;

    /**
     * Perform best model selection per the selected criterion.
     * Meant for use by \ref suzerain::arsel.
     */
    std::vector<real_t>::difference_type
    best_model(std::size_t          N,
               std::size_t          minorder,
               std::vector<real_t>& params,
               std::vector<real_t>& sigma2e,
               std::vector<real_t>& gain,
               std::vector<real_t>& autocor) const;

private:

    /** Tracks current criterion for reporting purposes. */
    std::string abbrev;

    /**
     * Encapsulates the best model function.
     * Void star ugliness reduces header coupling.
     */
    void * bmf;

};

} // namespace suzerain

#endif // SUZERAIN_SPECIFICATION_ARSEL_HPP
