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
#include <suzerain/timers.h>
#include <suzerain/zgbsv_specification.hpp>

namespace suzerain {

// TODO Document

class bsmbsm_solver : public suzerain_bsmbsm
{

public:

    zgbsv_specification spec;

    ArrayXXc LU;

    ArrayXc Pb;

    typedef Map<ArrayXXc, Aligned, OuterStride<Dynamic> > PAPT_type;

    PAPT_type PAPT;

    typedef Map<ArrayXc, Aligned> Px_type;

    Px_type Px;

    ArrayXi ipiv;

    int supply_b(const complex_double *b, int incb = 1)
    {
        SUZERAIN_TIMER_SCOPED("suzerain_bsmbsm_zaPxpby");
        return suzerain_bsmbsm_zaPxpby('N', S, n, 1, b, incb, 0, Pb.data(), 1);
    }

    int solve(const char trans)
    {
        SUZERAIN_TIMER_SCOPED(spec.mname());
        return solve_internal(trans);
    }

    int demand_x(complex_double *x, int incx = 1) const
    {
        SUZERAIN_TIMER_SCOPED("suzerain_bsmbsm_zaPxpby");
        return suzerain_bsmbsm_zaPxpby('T', S, n, 1, Px.data(), 1, 0, x, incx);
    }

protected:

    bsmbsm_solver(const suzerain_bsmbsm&     bsmbsm,
                  const zgbsv_specification& spec);

    virtual int solve_hook(const char trans) = 0;

private:

    int solve_internal(const char trans);

};

class bsmbsm_solver_zgbsv : public bsmbsm_solver
{
public:

    bsmbsm_solver_zgbsv(const suzerain_bsmbsm&     bsmbsm,
                        const zgbsv_specification& spec);

protected:

    virtual int solve_hook(const char trans);

};

class bsmbsm_solver_zgbsvx : public bsmbsm_solver
{
public:

    bsmbsm_solver_zgbsvx(const suzerain_bsmbsm&     bsmbsm,
                         const zgbsv_specification& spec);

    ArrayXr r;
    ArrayXr c;
    ArrayXc work;
    ArrayXr rwork;

protected:

    virtual int solve_hook(const char trans);

private:

    ArrayXXc PAPT_;
    ArrayXc  Px_;

};

class bsmbsm_solver_zcgbsvx : public bsmbsm_solver
{
public:

    bsmbsm_solver_zcgbsvx(const suzerain_bsmbsm&     bsmbsm,
                          const zgbsv_specification& spec);

    ArrayXc work;

protected:

    virtual int solve_hook(const char trans);

private:

    ArrayXXc PAPT_;
    ArrayXc  Px_;

};

} // namespace suzerain

#endif // SUZERAIN_BSMBSM_SOLVER_HPP
