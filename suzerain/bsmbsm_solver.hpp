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

#ifndef SUZERAIN_BSMBSM_SOLVER_HPP
#define SUZERAIN_BSMBSM_SOLVER_HPP

/** @file
 * Encapsulates solving BSMBSM problems per a \ref zgbsv_specification.
 * @see The definition of a BSMBSM problem in \ref bsmbsm.h.
 */

#include <suzerain/common.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/zgbsv_specification.hpp>

namespace suzerain {

// TODO Document

class bsmbsm_solver : public suzerain_bsmbsm, public zgbsv_specification
{

protected:

    bsmbsm_solver(const suzerain_bsmbsm&     bsmbsm,
                  const zgbsv_specification& specification);

public:

    ArrayXXc LU;

    ArrayXc b;

    typedef Map<ArrayXXc, Aligned, OuterStride<Dynamic> > A_type;

    A_type A;

    typedef Map<ArrayXc, Aligned> x_type;

    x_type x;

    ArrayXi ipiv;

};

class bsmbsm_solver_zgbsv : public bsmbsm_solver
{
public:

    bsmbsm_solver_zgbsv(const suzerain_bsmbsm&     bsmbsm,
                        const zgbsv_specification& specification);


};

class bsmbsm_solver_zgbsvx : public bsmbsm_solver
{
public:

    bsmbsm_solver_zgbsvx(const suzerain_bsmbsm&     bsmbsm,
                         const zgbsv_specification& specification);

    ArrayXr  r;
    ArrayXr  c;
    ArrayXc  work;
    ArrayXr  rwork;

private:

    ArrayXXc A_;
    ArrayXc  x_;

};

class bsmbsm_solver_zcgbsvx : public bsmbsm_solver
{
public:

    bsmbsm_solver_zcgbsvx(const suzerain_bsmbsm&     bsmbsm,
                          const zgbsv_specification& specification);

    ArrayXc  work;

private:

    ArrayXXc A_;
    ArrayXc  x_;

};

} // namespace suzerain

#endif // SUZERAIN_BSMBSM_SOLVER_HPP
