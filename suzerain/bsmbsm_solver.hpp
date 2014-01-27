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

#ifndef SUZERAIN_BSMBSM_SOLVER_HPP
#define SUZERAIN_BSMBSM_SOLVER_HPP

/** @file
 * Encapsulates solving BSMBSM problems per a \ref zgbsv_specification.
 * @see The definition of a BSMBSM problem in \ref bsmbsm.h.
 */

#include <suzerain/common.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/running_statistics.hpp>
#include <suzerain/specification_zgbsv.hpp>
#include <suzerain/timers.h>

namespace suzerain {

// TODO Too much logic is inlined in the bsmbsm_solver base class
//      (Inlining done so that timers are invoked only when available)

/**
 * An abstract base class for solving BSMBSM problems per a \ref
 * zgbsv_specification.
 *
 * To use the interface, users
 * <ol>
 *   <li>obtain a subclass instance,</li>
 *   <li>invoke #supply_B or #supply_b to set the right hand side(s),</li>
 *   <li>supply an unfactored operator into #PAPT per suzerain_bsmbsm.h,</li>
 *   <li>invoke #supplied_PAPT,</li>
 *   <li>invoke #solve to factorize and perform backsubstitution,</li>
 *   <li>and finally invoke #demand_X or #demand_x to obtain the solution.</li>
 * </ol>
 * If <tt>spec.in_place() == true</tt>, storage #LU and #PAPT are coincident as
 * are #PB and #PX.  For maximum flexibility, it should be assumed that calling
 * #solve destroys the contents of #PB and #PAPT.
 *
 * Users can, and should, reuse an instance for repeatedly solving problems of
 * some given size.  Factorizations will be reused when possible.  When a new
 * operator is of interest, invoke #supplied_PAPT again to inform the instance
 * that a new operator is stored in #PAPT.  If <tt>spec.reuse()</tt>,
 * subclasses may attempt to reuse the previous factorization unless
 * apprx(const bool) is called after #supplied_PAPT.  For example,
 * <tt>supplied_PAPT.apprx(false)</tt>.
 *
 * @see The definition of a BSMBSM problem in \ref bsmbsm.h.
 */
class bsmbsm_solver : public suzerain_bsmbsm
{

protected:

    /**
     * Public construction verboten.
     * Either construct a subclass or use the static build() method.
     */
    bsmbsm_solver(const suzerain_bsmbsm&     bsmbsm,
                  const zgbsv_specification& spec,
                  const int                  nrhs);

public:

    /**
     * Construct an instance for the given problem size,
     * solver specification, and maximum number of right hand sides.
     *
     * @param bsmbsm Defines the BSMBSM problem.
     * @param spec   Defines the solver behavior.
     * @param nrhs   Maximum number of right hand sides per #solve invocation.
     *
     * @return a <tt>new</tt> subclass instance matching the request.
     *         Callers must subsequently <tt>delete</tt> the instance.
     *
     * @throw <tt>std::invalid_argument</tt> if
     *        <tt>spec.method()</tt> is unknown.
     */
    static bsmbsm_solver* build(const suzerain_bsmbsm&     bsmbsm,
                                const zgbsv_specification& spec,
                                const int                  nrhs);

    /** Virtual destructor as appropriate for an abstract base class. */
    virtual ~bsmbsm_solver() { /* NOP */ }

    /** What solve specification does this instance use? */
    zgbsv_specification spec;

    /** What is the storage type for the LU factorization? */
    typedef Matrix<complex_double, Dynamic, Dynamic, ColMajor> LU_type;

    /** Storage for the LU factorization in conjunction with #ipiv. */
    LU_type LU;

    /** What is the storage type for permuted right hand sides \f$P B\f$? */
    typedef LU_type PB_type;

    /** Storage for permuted right hand sides \f$P B\f$. */
    PB_type PB;

    /** How is the \f$P A P^{\mbox{T}}\f$ operator stored? */
    typedef Map<LU_type, Aligned, OuterStride<Dynamic> > PAPT_type;

    /** Storage for the \f$P A P^{\mbox{T}}\f$ operator. */
    PAPT_type PAPT;

    /** What is the storage type for permuted solutions \f$P X\f$? */
    typedef Map<PB_type, Aligned> PX_type;

    /** Storage for permuted solutions \f$P X\f$. */
    PX_type PX;

    /** What is the storage type for LU-related pivots? */
    typedef Matrix<int, Dynamic, 1> ipiv_type;

    /** Storage for the factorization pivots in conjunction with #LU */
    ipiv_type ipiv;

    /**
     * Supply a non-permuted right hand side to #PB.
     * This routine automatically accounts for the BSMBSM-related renumbering.
     *
     * @param b    Data containing one right hand side.
     * @param j    Index of the right hand side.  After #solve, column \c j of
     *             #PX will contain the corresponding permuted solution.
     * @param incb Stride between adjacent elements of \c b.
     *
     * @return Zero on success.  Nonzero otherwise.
     */
    int supply_b(const complex_double *b, const int j, const int incb = 1)
    {
        assert(0 <= j && j < PB.cols());
        SUZERAIN_TIMER_SCOPED("suzerain_bsmbsm_zaPxpby");
        return suzerain_bsmbsm_zaPxpby('N', S, n, 1, b,                incb,
                                                  0, PB.col(j).data(), 1);
    }

    /**
     * Supply <tt>PB.cols()</tt> non-permuted right hand sides to #PB.
     * This routine automatically accounts for the BSMBSM-related renumbering.
     *
     * @param B    Data containing <tt>PB.cols()</tt> right hand sides.
     * @param ldB  Leading dimension between adjacent right hand sides.
     * @param incB Stride between adjacent elements in each right hand side.
     *
     * @return Zero on success.  Nonzero otherwise.
     */
    int supply_B(const complex_double *B, const int ldB, const int incB = 1)
    {
        int info = 0, j = -1;
        while (!info && ++j < PB.cols())
            info = supply_b(B + j*ldB, j, incB);
        return info;
    }

    /**
     * Supply <tt>PB.cols()</tt> non-permuted, contiguous right hand sides to
     * #PB.  This routine automatically accounts for the BSMBSM-related
     * renumbering.
     *
     * @param B Data containing <tt>PB.cols()</tt> right hand sides.
     *
     * @return Zero on success.  Nonzero otherwise.
     */
    int supply_B(const complex_double *B)
    {
        return supply_B(B, N, 1);
    }

    /**
     * Is a factorized operator available in #LU and #ipiv?
     *
     * @return 'N' if an operator is not factorized.
     *         Otherwise, some non-'N' subclass-dependent value.
     */
    char fact() const
    {
        return fact_;
    }

    /**
     * Should factorizations be equilibrated by default per #spec?
     *
     * @return 'E' if <tt>spec.equil() == true</tt>.
     *         Otherwise, 'N'.
     */
    char default_fact() const
    {
        return spec.equil() ? 'E' : 'N';
    }

    /** Could the next #solve invocation use an approximate factorization? */
    bool apprx() const
    {
        return apprx_;
    }

    /** Inform the solver that a new operator has been specified with #PAPT. */
    virtual bsmbsm_solver& supplied_PAPT();

    /**
     * Inform the solver that using an approximate factorization for the next
     * call to #solve is either acceptable or not.  This might be called
     * immediately after #supplied_PAPT().
     *
     * @return True if a factorization could have been reused
     *         prior to this invocation.  False otherwise.
     */
    virtual bool apprx(const bool acceptable);

    /**
     * Solve \f$ LU P X =  P B \f$ or \f$ {LU}^{\mbox{T}} P X = P B \f$.
     * Factorization is performed if appropriate.  Any subclass-specific,
     * right-hand-side-related outputs are only defined on indices <tt>0, ...,
     * nrhs - 1</tt>.
     *
     * @param trans Should \f$A\f$ or \f$A^{\mbox{T}}\f$ be used?
     * @param nrhs  Number of right hand sides up to a maximum of
     *              <tt>PB.cols()</tt>.
     *
     * @return Zero on success.  Nonzero otherwise.
     */
    int solve(const char trans, const int nrhs)
    {
        if (nrhs == 0) return 0;              // Avoid timing degenerate calls
        SUZERAIN_TIMER_SCOPED(spec.mname());
        return solve_internal(trans, nrhs);
    }

    /**
     * Solve \f$ LU P X =  P B \f$ or \f$ {LU}^{\mbox{T}} P X = P B \f$.
     *
     * @param trans Should \f$A\f$ or \f$A^{\mbox{T}}\f$ be used?
     *
     * @return Zero on success.  Nonzero otherwise.
     */
    int solve(const char trans)
    {
        return solve(trans, PB.cols());
    }

    /**
     * Demand a non-permuted solution from #PX.  This routine automatically
     * accounts for the BSMBSM-related renumbering.
     *
     * @param x    Destination for requested one right hand side.
     * @param j    Index of the right hand side to request.
     *             Before #solve, column \c j of
     *             #PB contained the corresponding permuted solution.
     * @param incx Stride between adjacent elements of \c x.
     *
     * @return Zero on success.  Nonzero otherwise.
     */
    int demand_x(complex_double *x, const int j, const int incx = 1) const
    {
        assert(0 <= j && j < PX.cols());
        SUZERAIN_TIMER_SCOPED("suzerain_bsmbsm_zaPxpby");
        return suzerain_bsmbsm_zaPxpby('T', S, n, 1, PX.col(j).data(), 1,
                                                  0, x,                incx);
    }

    /**
     * Demand<tt>PX.cols()</tt> non-permuted solutions from #PX.  This routine
     * automatically accounts for the BSMBSM-related renumbering.
     *
     * @param X    Destination for requested <tt>PX.cols</tt> right hand sides.
     * @param ldX  Leading dimension between adjacent right hand sides.
     * @param incX Stride between adjacent elements in each right hand side.
     *
     * @return Zero on success.  Nonzero otherwise.
     */
    int demand_X(complex_double *X, const int ldX, const int incX = 1) const
    {
        int info = 0, j = -1;
        while (!info && ++j < PB.cols())
            info = demand_x(X + j*ldX, j, incX);
        return info;
    }

    /**
     * Demand <tt>PX.cols()</tt> non-permuted, contiguous solutions from #PX.
     * This routine automatically accounts for the BSMBSM-related renumbering.
     *
     * @param X Destination for requested <tt>PX.cols</tt> right hand sides.
     *
     * @return Zero on success.  Nonzero otherwise.
     */
    int demand_X(complex_double *X) const { return demand_X(X, N, 1); }

    /**
     * Prepare a human-readable summary of the solution procedure.
     * Possibly prepare additional output describing operational statistics.
     *
     * @return Each string is one self-contained message suitable for logging.
     */
    virtual std::vector<std::string> summarize_statistics() const;

protected:

    virtual int solve_hook(const char trans, const int nrhs) = 0;

    char fact_;

    int apprx_;

private:

    int solve_internal(const char trans, const int nrhs);

};

/**
 * Solve BSMBSM problems using LAPACK's ZGBSV, ZGBTRF, and ZGBTRS.
 */
class bsmbsm_solver_zgbsv : public bsmbsm_solver
{
public:

    /**
     * Construct an instance for the given problem size,
     * solver specification, and maximum number of right hand sides.
     *
     * @param bsmbsm Defines the BSMBSM problem.
     * @param spec   Defines the solver behavior.
     * @param nrhs   Maximum number of right hand sides per #solve invocation.
     *
     * @throw <tt>std::invalid_argument</tt> if
     *        <tt>spec.method() != zgbsv_specification::zgbsv</tt>.
     */
    bsmbsm_solver_zgbsv(const suzerain_bsmbsm&     bsmbsm,
                        const zgbsv_specification& spec,
                        const int                  nrhs);

protected:

    virtual int solve_hook(const char trans, const int nrhs);

};

/**
 * Solve BSMBSM problems using LAPACK's ZGBSVX.
 */
class bsmbsm_solver_zgbsvx : public bsmbsm_solver
{
    /** Type of the contiguous storage housing r, c, and rwork. */
    typedef Matrix<double, Dynamic, 3, ColMajor> rcwork_type;

    /** Type of the contiguous storage housing ferr and berr. */
    typedef Matrix<double, 2, Dynamic, ColMajor> err_type;

public:

    /**
     * Construct an instance for the given problem size,
     * solver specification, and maximum number of right hand sides.
     *
     * @param bsmbsm Defines the BSMBSM problem.
     * @param spec   Defines the solver behavior.
     * @param nrhs   Maximum number of right hand sides per #solve invocation.
     *
     * @throw <tt>std::invalid_argument</tt> if
     *        <tt>spec.method() != zgbsv_specification::zgbsvx</tt>.
     */
    bsmbsm_solver_zgbsvx(const suzerain_bsmbsm&     bsmbsm,
                         const zgbsv_specification& spec,
                         const int                  nrhs);

    /** @return \c equed from any prior suzerain_lapack_zgbsvx() usage. */
    char equed() const
    {
        return equed_;
    }

    /** @return \c rcond from any prior suzerain_lapack_zgbsvx() usage. */
    double rcond() const
    {
        return rcond_;
    }

    /** @return <tt>1/cond</tt> from any prior suzerain_lapack_zgbsvx() usage.*/
    double cond() const
    {
        return 1 / rcond_;
    }

    /** @return <tt>r</tt> from any prior suzerain_lapack_zgbsvx() usage. */
    rcwork_type::ConstColXpr r()     const { return rcwork_.col(0); }

    /** @return <tt>c</tt> from any prior suzerain_lapack_zgbsvx() usage. */
    rcwork_type::ConstColXpr c()     const { return rcwork_.col(1); }

    /** @return <tt>rwork</tt> from any prior suzerain_lapack_zgbsvx() usage. */
    rcwork_type::ConstColXpr rwork() const { return rcwork_.col(2); }

    /** @return <tt>ferr</tt> from any prior suzerain_lapack_zgbsvx() usage. */
    err_type::ConstRowXpr ferr() const { return err_.row(0); }

    /** @return <tt>berr</tt> from any prior suzerain_lapack_zgbsvx() usage. */
    err_type::ConstRowXpr berr() const { return err_.row(1); }

    /** Type for tracking statistics on various quantities after each solve. */
    typedef running_statistics<double, 5> stats_type;

    /** Human-readable names of the statistics which are tracked. */
    static const char * const stats_names[stats_type::static_size];

    /** Tracks statistics on various quantities after each solve. */
    stats_type stats;

    virtual std::vector<std::string> summarize_statistics() const;

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

/**
 * Solve BSMBSM problems using suzerain_lapackext_zcgbsvx().
 */
class bsmbsm_solver_zcgbsvx : public bsmbsm_solver
{
    /** Type of the contiguous storage housing siter and diter. */
    typedef Matrix<int, 2, Dynamic, ColMajor> iter_type;

    /** Type of the contiguous storage housing tolsc and res. */
    typedef Matrix<double, 2, Dynamic, ColMajor> tolscres_type;

public:

    /**
     * Construct an instance for the given problem size,
     * solver specification, and maximum number of right hand sides.
     *
     * @param bsmbsm Defines the BSMBSM problem.
     * @param spec   Defines the solver behavior.
     * @param nrhs   Maximum number of right hand sides per #solve usage.
     *
     * @throw <tt>std::invalid_argument</tt> if
     *        <tt>spec.method() != zgbsv_specification::zcgbsvx</tt>.
     */
    bsmbsm_solver_zcgbsvx(const suzerain_bsmbsm&     bsmbsm,
                          const zgbsv_specification& spec,
                          const int                  nrhs);

    virtual bsmbsm_solver& supplied_PAPT();

    /** @return \c afrob from any prior suzerain_lapackext_zcgbsvx() usage. */
    double afrob() const
    {
        return afrob_;
    }

    /** @return \c siter from any prior suzerain_lapackext_zcgbsvx() usage. */
    iter_type::ConstRowXpr siter() const { return iter_.row(0); }

    /** @return \c diter from any prior suzerain_lapackext_zcgbsvx() usage. */
    iter_type::ConstRowXpr diter() const { return iter_.row(1); }

    /** @return \c tolsc from any prior suzerain_lapackext_zcgbsvx() usage. */
    tolscres_type::ConstRowXpr tolsc() const { return tolscres_.row(0); }

    /** @return \c res from any prior suzerain_lapackext_zcgbsvx() usage. */
    tolscres_type::ConstRowXpr res()   const { return tolscres_.row(1); }

    /** Type for tracking statistics on various quantities after each solve. */
    typedef running_statistics<double, 8> stats_type;

    /** Human-readable names of the statistics which are tracked. */
    static const char * const stats_names[stats_type::static_size];

    /** Tracks statistics on various quantities after each solve. */
    stats_type stats;

    virtual std::vector<std::string> summarize_statistics() const;

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
