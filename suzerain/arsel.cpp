//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
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
 * @copydoc arsel.hpp
 */

#include <suzerain/arsel.hpp>

#include <suzerain/ar.hpp>
#include <suzerain/specification_arsel.hpp>

namespace suzerain {

running_statistics<double,1>
arsel(const std::size_t Ny,
      const std::vector<double>& t,
      const double * data,
      const std::size_t ld,
      const specification_arsel& spec,
      std::vector<double>& eff_N,
      std::vector<double>& eff_var,
      std::vector<double>& mu,
      std::vector<double>& mu_sigma,
      std::vector<double>& p,
      std::vector<double>& T)
{
    const std::size_t Nt = t.size();
    if (Ny == 0) throw std::invalid_argument("arsel received Ny == 0");
    if (Nt == 0) throw std::invalid_argument("arsel received Nt == 0");
    if (ld < Ny) throw std::invalid_argument("arsel received ld < Ny");

    // Prepare output vectors to gather burg_method()
    eff_N   .resize(Ny);
    eff_var .resize(Ny);
    mu      .resize(Ny);
    mu_sigma.resize(Ny);
    p       .resize(Ny);
    T       .resize(Ny);

    // Compute sampling rate statistics
    running_statistics<double,1> dt;
    for (std::size_t i = 1; i < t.size(); ++i) {
        const double diff = t[i] - t[i-1];
        dt(&diff);
    }

    // Prepare vectors to capture burg_method() output
    std::vector<double> params, sigma2e, gain, autocor;
    params .reserve(spec.maxorder*(spec.maxorder + 1)/2);
    sigma2e.reserve(spec.maxorder + 1);
    gain   .reserve(spec.maxorder + 1);
    autocor.reserve(spec.maxorder + 1);

    // Prepare repeatedly-used working storage for burg_method().
    std::vector<double> f, b, Ak, ac;

    for (std::size_t j = 0; j < Ny; ++j) {

        std::size_t maxorder = spec.maxorder;
        params .clear();
        sigma2e.clear();
        gain   .clear();
        autocor.clear();
        ar::burg_method(ar::strided_adaptor<const double*>(&data[j+ 0*ld], ld),
                        ar::strided_adaptor<const double*>(&data[j+Nt*ld], ld),
                        mu[j],
                        maxorder,
                        back_inserter(params),
                        back_inserter(sigma2e),
                        back_inserter(gain),
                        back_inserter(autocor),
                        spec.submean,
                        true /* output hierarchy? */,
                        f, b, Ak, ac);

        spec.best_model(Nt, params, sigma2e, gain, autocor);

        const double T0 = ar::decorrelation_time(
                static_cast<std::size_t>(spec.wlenT0*Nt),
                ar::autocorrelation(params.begin(), params.end(),
                                    gain[0], autocor.begin()),
                spec.absrho);

        eff_var[j]  = (Nt*gain[0]*sigma2e[0]) / (Nt - T0);
        eff_N[j]    = Nt / T0;
        mu_sigma[j] = std::sqrt(eff_var[j] / eff_N[j]);
        p[j]        = params.size();
        T[j]        = T0 * dt.avg(0);

    }

    return dt;
}

} // namespace suzerain
