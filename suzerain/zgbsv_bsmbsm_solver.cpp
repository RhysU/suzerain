//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc zgbsv_bsmbsm_solver.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/zgbsv_bsmbsm_solver.hpp>

namespace suzerain {

// TODO Implement

zgbsv_bsmbsm_solver::zgbsv_bsmbsm_solver(const suzerain_bsmbsm& bsmbsm)
    : suzerain_bsmbsm(bsmbsm) // Copy the given BSMBSM problem size
    , A(0, 0, 0)              // See Eigen "Changing the mapped array" docs
    , LU(0, 0, 0)             // Ditto
    , ipiv(bsmbsm.N)          // Required storage known
    , b(bsmbsm.N)             // Required storage known
    , x(b.data(), bsmbsm.N)   // Changeable by subclass, defaults to b
{
    // NOP
}

} // end namespace suzerain
