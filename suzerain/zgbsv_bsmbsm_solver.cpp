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

zgbsv_bsmbsm_solver::zgbsv_bsmbsm_solver(
        const zgbsv_specification& specification,
        const suzerain_bsmbsm&    bsmbsm)
    : zgbsv_specification(specification)
    , suzerain_bsmbsm(bsmbsm)
    , A(0, 0, 0)       // Degenerate, requires subclass override
    , LU(0, 0, 0)      // Degenerate, requires subclass override
    , ipiv(N)          // Required storage already known per bsmbsm
    , b(N)             // Required storage already known per bsmbsm
    , x(0, 0)          // Degenerate, requires subclass override
{
    // NOP
}

zgbsv_bsmbsm_solver_zgbsv::zgbsv_bsmbsm_solver_zgbsv(
        const zgbsv_specification& specification,
        const suzerain_bsmbsm&    bsmbsm)
    : zgbsv_bsmbsm_solver(specification, bsmbsm)
{
    // See Eigen "Changing the mapped array" for details on placement new calls
    if (method() != zgbsv_specification::zgbsv) throw std::invalid_argument(
            "Invalid method() in zgbsv_bsmbsm_solver_zgbsv");

    // Allocate storage for in-place LU which also stores A pre-factorization
    // Beginning A offset different from LU per LAPACK ZGBTRF requirements
    buf.resize(LD + KL, N);
    new (&A)  A_type (buf.data() + KL, LD,         buf.cols(), buf.colStride());
    new (&LU) LU_type(buf.data(),      buf.rows(), buf.cols()                 );

    // As solve is in-place, output x simply aliases right hand side b
    assert(in_place() == true);
    new (&x) x_type(b.data(), N);
}

} // end namespace suzerain
