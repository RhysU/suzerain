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

/** @file
 * @copydoc specification_arsel.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/specification_arsel.hpp>

#include <suzerain/ar.hpp>

namespace suzerain {

specification_arsel::specification_arsel(
        const bool         absrho,
        const std::string& criterion,
        const std::size_t  maxorder,
        const std::size_t  minorder,
        const double       wlenT0)
    : absrho(absrho)
    , maxorder(maxorder)
    , minorder(minorder)
    , wlenT0(wlenT0)
    , bmf(NULL)
{
    this->criterion(criterion);
}

// Once-and-for-all typedef around the void star crimes we commit
typedef ar::best_model_function<
            ar::Burg, std::size_t, std::size_t, std::vector<double>
        > best_model_function;
// Then make sure the coppers will never find us, never.  We hope.
enum {
    the_crime = 1/(sizeof(best_model_function::type) == sizeof(void *))
};

const std::string&
specification_arsel::criterion(const std::string& abbrev)
{
    const best_model_function::type best_model
            = best_model_function::lookup(abbrev, submean);
    if (!best_model) {
        throw std::invalid_argument(
                "Unknown model selection criterion abbreviation");
    }
    bmf = reinterpret_cast<void *>(best_model);
    return this->abbrev = abbrev;
}

const std::string&
specification_arsel::criterion() const
{
    return this->abbrev;
}

std::vector<double>::difference_type
specification_arsel::best_model(
        std::size_t          N,
        std::vector<double>& params,
        std::vector<double>& sigma2e,
        std::vector<double>& gain,
        std::vector<double>& autocor) const
{
    return reinterpret_cast<best_model_function::type>(bmf)(
            N, minorder, params, sigma2e, gain, autocor);
}

} // namespace suzerain
