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
#include <suzerain/running_statistics.hpp>
#include <suzerain/timers.h>
#include <suzerain/zgbsv_specification.hpp>

namespace suzerain {

// TODO Document

class bsmbsm_solver : public suzerain_bsmbsm
{

public:

    zgbsv_specification spec;

    typedef Matrix<complex_double, Dynamic, Dynamic, ColMajor> LU_type;

    LU_type LU;

    typedef LU_type PB_type;

    PB_type PB;

    typedef Map<LU_type, Aligned, OuterStride<Dynamic> > PAPT_type;

    PAPT_type PAPT;

    typedef Map<PB_type, Aligned> PX_type;

    PX_type PX;

    typedef Matrix<int, Dynamic, 1> ipiv_type;

    ipiv_type ipiv;

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

    virtual bsmbsm_solver& supplied_PAPT();

    virtual bool apprx(const bool value);

    int solve(const char trans, const int nrhs)
    {
        SUZERAIN_TIMER_SCOPED(spec.mname());
        return solve_internal(trans, nrhs);
    }

    int solve(const char trans)
    {
        return solve(trans, PB.cols());
    }

    // TODO Add statistics reporting capabilities as virtual method

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

    int apprx_;

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
    /** Type of the contiguous storage housing r, c, and rwork. */
    typedef Matrix<double, Dynamic, 3, ColMajor> rcwork_type;

    /** Type of the contiguous storage housing ferr and berr. */
    typedef Matrix<double, 2, Dynamic, ColMajor> err_type;

public:

    bsmbsm_solver_zgbsvx(const suzerain_bsmbsm&     bsmbsm,
                         const zgbsv_specification& spec,
                         const int                  nrhs);

    char equed() const
    {
        return equed_;
    }

    double rcond() const
    {
        return rcond_;
    }

    double cond() const
    {
        return 1 / rcond_;
    }

    rcwork_type::ConstColXpr r()     const { return rcwork_.col(0); }
    rcwork_type::ConstColXpr c()     const { return rcwork_.col(1); }
    rcwork_type::ConstColXpr rwork() const { return rcwork_.col(2); }

    err_type::ConstRowXpr ferr() const { return err_.row(0); }
    err_type::ConstRowXpr berr() const { return err_.row(1); }

    // Track statistics on various quantities after each solve
    typedef running_statistics<double, 5> stats_type;
    static const char * const stats_names[stats_type::static_size];
    stats_type stats;

protected:

    virtual int solve_hook(const char trans, const int nrhs);

    char equed_;

    rcwork_type rcwork_;
    rcwork_type::ColXpr r_()     { return rcwork_.col(0); }
    rcwork_type::ColXpr c_()     { return rcwork_.col(1); }
    rcwork_type::ColXpr rwork_() { return rcwork_.col(2); }


    double rcond_;

    err_type err_;
    err_type::RowXpr ferr_() { return err_.row(0); }
    err_type::RowXpr berr_() { return err_.row(1); }

    Matrix<complex_double, Dynamic, 2, ColMajor> work_;

private:

    LU_type PAPT_;

    PB_type PX_;
};

class bsmbsm_solver_zcgbsvx : public bsmbsm_solver
{
    /** Type of the contiguous storage housing siter and diter. */
    typedef Matrix<int, 2, Dynamic, ColMajor> iter_type;

    /** Type of the contiguous storage housing tolsc and res. */
    typedef Matrix<double, 2, Dynamic, ColMajor> tolscres_type;

public:

    bsmbsm_solver_zcgbsvx(const suzerain_bsmbsm&     bsmbsm,
                          const zgbsv_specification& spec,
                          const int                  nrhs);

    virtual bsmbsm_solver& supplied_PAPT();

    iter_type::ConstRowXpr siter() const { return iter_.row(0); }
    iter_type::ConstRowXpr diter() const { return iter_.row(1); }

    tolscres_type::ConstRowXpr tolsc() const { return tolscres_.row(0); }
    tolscres_type::ConstRowXpr res()   const { return tolscres_.row(1); }

    // Track statistics on various quantities after each solve
    typedef running_statistics<double, 8> stats_type;
    static const char * const stats_names[stats_type::static_size];
    stats_type stats;

protected:

    virtual int solve_hook(const char trans, const int nrhs);

    double afrob_;

    iter_type iter_;
    iter_type::RowXpr siter_() { return iter_.row(0); }
    iter_type::RowXpr diter_() { return iter_.row(1); }

    tolscres_type tolscres_;
    tolscres_type::RowXpr tolsc_() { return tolscres_.row(0); }
    tolscres_type::RowXpr res_()   { return tolscres_.row(1); }

    Matrix<complex_double, Dynamic, 1> r_;

private:

    LU_type PAPT_;

    PB_type PX_;

};

} // namespace suzerain

#endif // SUZERAIN_BSMBSM_SOLVER_HPP
