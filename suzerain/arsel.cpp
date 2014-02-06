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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/arsel.hpp>
#include <suzerain/ar.hpp>
#include <suzerain/running_statistics.hpp>

namespace suzerain {

real_t
arsel(const std::vector<real_t> t,
      const ArrayXXr& data,
      const specification_arsel& spec,
      std::vector<real_t>& eff_N,
      std::vector<real_t>& eff_var,
      std::vector<real_t>& mu,
      std::vector<real_t>& mu_sigma,
      std::vector<real_t>& p,
      std::vector<real_t>& T)
{
    size_t Ny, Nt;
    SUZERAIN_ENSURE((Ny = data.rows()) > 0);
    SUZERAIN_ENSURE((Nt = data.cols()) > 0);
    SUZERAIN_ENSURE(t.size() == Nt);

    // Prepare output vectors to gather burg_method()
    eff_N   .resize(Ny);
    eff_var .resize(    Ny);
    mu      .resize(    Ny);
    mu_sigma.resize(    Ny);
    p       .resize(    Ny);
    T       .resize(    Ny);

    // Compute sampling rate statistics
    running_statistics<real_t,1> dt;
    for (size_t i = 1; i < t.size(); ++i) {
        const real_t diff = t[i] - t[i-1];
        dt(&diff);
    }

    // Prepare vectors to capture burg_method() output
    std::vector<real_t> params, sigma2e, gain, autocor;
    params .reserve(spec.maxorder*(spec.maxorder + 1)/2);
    sigma2e.reserve(spec.maxorder + 1);
    gain   .reserve(spec.maxorder + 1);
    autocor.reserve(spec.maxorder + 1);

    // Prepare repeatedly-used working storage for burg_method().
    std::vector<real_t> f, b, Ak, ac;

    for (size_t j = 0; j < Ny; ++j) {

        std::size_t maxorder = spec.maxorder;
        params .clear();
        sigma2e.clear();
        gain   .clear();
        autocor.clear();
        ar::strided_adaptor<const real_t*> signal_begin(&data.coeff(j, 0),Ny);
        ar::strided_adaptor<const real_t*> signal_end  (&data.coeff(j,Nt),Ny);
        ar::burg_method(signal_begin,
                        signal_end,
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

        const real_t T0 = ar::decorrelation_time(
                static_cast<size_t>(spec.wlenT0*Nt),
                ar::autocorrelation(params.begin(), params.end(),
                                    gain[0], autocor.begin()),
                spec.absrho);

        eff_var[j]  = (Nt*gain[0]*sigma2e[0]) / (Nt - T0);
        eff_N[j]    = Nt / T0;
        mu_sigma[j] = std::sqrt(eff_var[j] / eff_N[j]);
        p[j]        = params.size();
        T[j]        = T0 * dt.avg(0);

    }

    return dt.avg(0);
}

} // namespace suzerain
