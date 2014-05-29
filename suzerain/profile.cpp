//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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
 * @copydoc profile.hpp
 */

#include <suzerain/profile.hpp>

#include <suzerain/samples.hpp>
#include <suzerain/summary.hpp>

namespace suzerain {

profile::profile()
{
}

profile::profile(profile::storage_type::Index Ny)
    : storage(storage_type::Zero(Ny, storage_type::ColsAtCompileTime))
{
}

profile& profile::operator=(const samples &q)
{
    // Resize our storage and defensively NaN in case we miss something
    this->storage.setConstant(q.storage.rows(),
            storage_type::ColsAtCompileTime,
            std::numeric_limits<storage_type::Scalar>::quiet_NaN());

    // Copy profiles of interest from q
    this->rho()          = q.rho();
    this->rho_u()        = q.rho_u();
    this->rho_E()        = q.rho_E();
    this->a()            = q.a();
    this->H0()           = q.H0();
    this->ke()           = q.ke();
    this->T()            = q.T();
    this->mu()           = q.mu();
    this->u()            = q.u();

    return *this;
}

profile& profile::operator=(const summary &q)
{
    // Resize our storage and defensively NaN in case we miss something
    this->storage.setConstant(q.storage.rows(),
            storage_type::ColsAtCompileTime,
            std::numeric_limits<storage_type::Scalar>::quiet_NaN());

    // Copy profiles of interest from q
    this->rho()          = q.bar_rho();
    this->rho_u().col(0) = q.bar_rho_u();
    this->rho_u().col(1) = q.bar_rho_v();
    this->rho_u().col(2) = q.bar_rho_w();
    this->rho_E()        = q.bar_rho_E();
    this->a()            = q.bar_a();
    this->H0()           = q.bar_H0();
    this->ke()           = q.bar_ke();
    this->T()            = q.bar_T();
    this->mu()           = q.bar_mu();
    this->u().col(0)     = q.bar_u();
    this->u().col(1)     = q.bar_v();
    this->u().col(2)     = q.bar_w();

    return *this;
}

} // namespace suzerain
