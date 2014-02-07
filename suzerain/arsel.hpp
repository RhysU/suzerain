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

#ifndef SUZERAIN_ARSEL_HPP
#define SUZERAIN_ARSEL_HPP

/** @file
 * Automatically fit autoregressive models to input signals.
 * Serves as a high-level wrapper of functionality within \ref ar.hpp.
 */

#include <cstddef>
#include <vector>

#include <suzerain/running_statistics.hpp>

namespace suzerain {

// Forward declarations
class specification_arsel;

/**
 * Automatically fit autoregressive models to input signals. Use ar::burg_method
 * and ar::best_model to fit an autoregressive process for signals contained in
 * the rows of matrix data.
 *
 * \param[in   Ny       Number of independent signals to process.
 * \param[in]  t        Times of input data samples also informing number of
 *                      samples.  Though \e only equispaced data is
 *                      currently supported, this information will be
 *                      used to compute a mean sample rate and therefore
 *                      inform output \c T.
 * \param[in]  data     Column-major <code>Ny</code> by <code>t.size()</code>
 *                      data to process.
 * \param[in   ld       Leading dimension between columns of \c data.
 * \param[in]  spec     Specification governing algorithmic choices.
 * \param[out] eff_N    Number of effectively independent samples in each row.
 * \param[out] eff_var  Estimated effective signal variance for each row.
 * \param[out] mu       Sample mean for each row.
 * \param[out] mu_sigma Sampling error in each row.
 *                      (that is, the standard deviation of the sample mean).
 * \param[out] p        Selected autoregressive model order for each row.
 * \param[out] T        Decorrelation time for each row
 *                      computed from decorrelatio separation \f$T_0\f$
 *                      and the mean sampling rate in \c t.
 *
 * \returns Statistical information about input sample times \c t.
 * \see \ref ar for more details.
 */
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
      std::vector<double>& T);

} // end namespace suzerain

#endif /* SUZERAIN_ARSEL_HPP */
