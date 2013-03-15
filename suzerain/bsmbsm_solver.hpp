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

    ArrayXXc PB;

    typedef Map<ArrayXXc, Aligned, OuterStride<Dynamic> > PAPT_type;

    PAPT_type PAPT;

    typedef Map<ArrayXXc, Aligned> PX_type;

    PX_type PX;

    ArrayXi ipiv;

    int supply_b(const complex_double *b, const int j, const int incb = 1)
    {
        SUZERAIN_TIMER_SCOPED("suzerain_bsmbsm_zaPxpby");
        return suzerain_bsmbsm_zaPxpby('N', S, n, 1, b,                incb,
                                                  0, PB.col(j).data(), 1);
    }

    int supply_B(const complex_double *B, const int ldB, const int incB = 1)
    {
        int info = 0, j = -1;
        while (!info && ++j < PB.cols())
            info = supply_b(B + j*ldB, j, incB);
        return info;
    }

    int supply_B(const complex_double *B)
    {
        return supply_B(B, N, 1);
    }

    char fact() const
    {
        return fact_;
    }

    char default_fact() const
    {
        return spec.equil() ? 'E' : 'N';
    }

    bool apprx() const
    {
        return apprx_;
    }

    bool apprx(const bool value)
    {
        const bool old = apprx_;
        if (spec.reuse())
            apprx_ = value;
        return old;
    }

    void supplied_PAPT()
    {
        if (spec.reuse()) {
            apprx_ = fact_ != default_fact();
        } else {
            fact_ = default_fact();
        }
    }

    int solve(const char trans, const int nrhs)
    {
        SUZERAIN_TIMER_SCOPED(spec.mname());
        return solve_internal(trans, nrhs);
    }

    int solve(const char trans)
    {
        return solve(trans, PB.cols());
    }

    int demand_x(complex_double *x, const int j, const int incx = 1) const
    {
        SUZERAIN_TIMER_SCOPED("suzerain_bsmbsm_zaPxpby");
        return suzerain_bsmbsm_zaPxpby('T', S, n, 1, PX.col(j).data(), 1,
                                                  0, x,                incx);
    }

    int demand_X(complex_double *X, const int ldX, const int incX = 1) const
    {
        int info = 0, j = -1;
        while (!info && ++j < PB.cols())
            info = demand_x(X + j*ldX, j, incX);
        return info;
    }

    int demand_X(complex_double *X) const { return demand_X(X, N, 1); }

protected:

    bsmbsm_solver(const suzerain_bsmbsm&     bsmbsm,
                  const zgbsv_specification& spec,
                  const int                  nrhs);

    virtual int solve_hook(const char trans, const int nrhs) = 0;

    char fact_;

    int  apprx_;

private:

    int solve_internal(const char trans, const int nrhs);

};

class bsmbsm_solver_zgbsv : public bsmbsm_solver
{
public:

    bsmbsm_solver_zgbsv(const suzerain_bsmbsm&     bsmbsm,
                        const zgbsv_specification& spec,
                        const int                  nrhs);

protected:

    virtual int solve_hook(const char trans, const int nrhs);

};

class bsmbsm_solver_zgbsvx : public bsmbsm_solver
{
public:

    bsmbsm_solver_zgbsvx(const suzerain_bsmbsm&     bsmbsm,
                         const zgbsv_specification& spec,
                         const int                  nrhs);

    ArrayXr r;
    ArrayXr c;
    ArrayXc work;
    ArrayXr rwork;

protected:

    virtual int solve_hook(const char trans, const int nrhs);

private:

    ArrayXXc PAPT_;
    ArrayXXc PX_;

};

class bsmbsm_solver_zcgbsvx : public bsmbsm_solver
{
public:

    bsmbsm_solver_zcgbsvx(const suzerain_bsmbsm&     bsmbsm,
                          const zgbsv_specification& spec,
                          const int                  nrhs);

    ArrayXc work;

protected:

    virtual int solve_hook(const char trans, const int nrhs);

private:

    ArrayXXc PAPT_;
    ArrayXXc PX_;

};

} // namespace suzerain

#endif // SUZERAIN_BSMBSM_SOLVER_HPP
