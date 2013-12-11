//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_PERFECT_QUANTITIES_HPP
#define SUZERAIN_PERFECT_QUANTITIES_HPP

/** @file
 * Sampling logistics for mean quantity profiles.
 */

#include <suzerain/common.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/support/esio_fwd.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class operator_tools;

namespace perfect {

// Forward declarations
class scenario_definition;

/**
 * Encapsulate the mean quantities detailed in the "Sampling logistics" section
 * of <tt>writeups/perfectgas.tex</tt>.
 *
 * Samples of each quantity are made available through two-dimensional,
 * column-major arrays.  The row index iterates over wall-normal B-spline
 * coefficients and the column index iterates over tensor indices.  Scalars
 * have only a single tensor index.  Vector quantities have three indices
 * corresponding to the streamwise x, wall-normal y, and spanwise z directions.
 * Symmetric tensors (for example, \f$\overline{\mu{}S}\f$}) have six entries
 * corresponding to the <tt>xx</tt>, <tt>xy</tt>, <tt>xz</tt>, <tt>yy</tt>,
 * <tt>yz</tt>, and <tt>zz</tt> indices.  Rank one triple products (for example
 * \f$\overline{\rho{}u\otimes{}u\otimes{}u}\f$) have ten entries corresponding
 * to the <tt>xxx</tt>, <tt>xxy</tt>, <tt>xxz</tt>, <tt>xyy</tt>, <tt>xyz</tt>,
 * <tt>xzz</tt>, <tt>yyy</tt>, <tt>yyz</tt>, <tt>yzz</tt>, and <tt>zzz</tt>
 * indices.
 *
 * \internal Many memory overhead and implementation consistency issues have
 * been traded for the headache of reading Boost.Preprocessor-based logic.  So
 * it goes.
 */
class quantities
{
public:

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed in wave
 * space.
 */
#define SUZERAIN_PERFECT_QUANTITIES_WAVE                   \
    ((rho,                      1)) /* scalar           */ \
    ((rho_u,                    3)) /* vector           */ \
    ((rho_E,                    1)) /* scalar           */

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed in physical
 * space.
 */
#define SUZERAIN_PERFECT_QUANTITIES_PHYSICAL         \
    ((E,                 1))  /* scalar           */ \
    ((T,                 1))  /* scalar           */ \
    ((a,                 1))  /* scalar           */ \
    ((h0,                1))  /* scalar           */ \
    ((H0,                1))  /* scalar           */ \
    ((mu,                1))  /* scalar           */ \
    ((nu,                1))  /* scalar           */ \
    ((u,                 3))  /* vector           */ \
    ((sym_grad_u,        6))  /* symmetric tensor */ \
    ((sym_rho_grad_u,    6))  /* symmetric tensor */ \
    ((grad_T,            3))  /* vector           */ \
    ((rho_grad_T,        3))  /* vector           */ \
    ((tau_colon_grad_u,  1))  /* scalar           */ \
    ((tau,               6))  /* symmetric tensor */ \
    ((tau_u,             3))  /* vector           */ \
    ((p_div_u,           1))  /* scalar           */ \
    ((rho_u_u,           6))  /* symmetric tensor */ \
    ((rho_u_u_u,        10))  /* symmetric tensor */ \
    ((rho_T_u,           3))  /* vector           */ \
    ((rho_mu,            1))  /* scalar           */ \
    ((mu_S,              6))  /* symmetric tensor */ \
    ((mu_div_u,          1))  /* scalar           */ \
    ((mu_grad_T,         3))  /* vector           */

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed
 * through implicit forcing.
 */
#define SUZERAIN_PERFECT_QUANTITIES_IMPLICIT         \
    ((SrhoE,             1))  /* scalar           */ \
    ((Srhou,             3))  /* vector           */ \
    ((Srho,              1))  /* scalar           */ \
    ((Srhou_dot_u,       1))  /* scalar           */ \
    ((f,                 3))  /* vector           */ \
    ((f_dot_u,           1))  /* scalar           */ \
    ((qb,                1))  /* scalar           */ \
    ((CrhoE,             1))  /* scalar           */ \
    ((Crhou,             3))  /* vector           */ \
    ((Crho,              1))  /* scalar           */ \
    ((Crhou_dot_u,       1))  /* scalar           */

/** A Boost.Preprocessor sequence of tuples of all sampled quantities. */
#define SUZERAIN_PERFECT_QUANTITIES      \
    SUZERAIN_PERFECT_QUANTITIES_WAVE     \
    SUZERAIN_PERFECT_QUANTITIES_PHYSICAL \
    SUZERAIN_PERFECT_QUANTITIES_IMPLICIT

    /* Compile-time totals of the number of scalars sampled at each point */
    struct nscalars { enum {
#define EXTRACT(r, data, tuple) BOOST_PP_TUPLE_ELEM(2, 1, tuple)
#define SUM(s, state, x) BOOST_PP_ADD(state, x)

        wave = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_QUANTITIES_WAVE)),

        physical = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_QUANTITIES_PHYSICAL)),

        implicit = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_QUANTITIES_IMPLICIT)),

        total = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    EXTRACT,,SUZERAIN_PERFECT_QUANTITIES))

#undef EXTRACT
#undef SUM
    }; };

    /** Simulation time when mean quantities were obtained */
    real_t t;

    /** Type of the contiguous storage used to house all scalars */
    typedef Array<real_t, Dynamic, nscalars::total> storage_type;

    /** Contiguous storage used to house all means */
    storage_type storage;

    /**
     * Constructor setting <tt>this->t = NaN</tt>.
     * Caller will need to resize <tt>this->storage</tt> prior to use.
     */
    quantities();

    /**
     * Constructor setting <tt>this->t = t</tt>.
     * Caller will need to resize <tt>this->storage</tt> prior to use.
     */
    explicit quantities(real_t t);

    /**
     * Constructor setting <tt>this->t = t</tt> and preparing a zero-filled \c
     * storage containing \c Ny rows.
     */
    quantities(real_t t, storage_type::Index Ny);

    /**
     * Save quantities to a restart file.
     *
     * @return True if all quantities could be saved.  False otherwise.
     */
    bool save(const esio_handle h) const;

    /**
     * Load quantities from a restart file.  Statistics not present in the
     * restart file are considered to be all NaNs.  Member #t, which is not
     * modified by this routine, is presumably set in some other fashion.
     *
     * @return True if all quantities could be loaded.  False otherwise.
     */
    bool load(const esio_handle h);

#define OP(r, data, tuple)                                              \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple) = BOOST_PP_TUPLE_ELEM(2, 1, tuple)

    /** Compile-time offsets for each quantity within \c storage */
    struct start { enum {
        wave     = 0,                              // Start of wave block
        physical = wave     + nscalars::wave,      // Start of physical block
        implicit = physical + nscalars::physical,  // Start of implicit block
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(
                OP,,SUZERAIN_SHIFTED_SUM(SUZERAIN_PERFECT_QUANTITIES)))
    }; };

    /**
     * Provide access to contiguous subregions within storage.
     * @{
     */
    storage_type::NColsBlockXpr<nscalars::wave>::Type wave()
    { return storage.middleCols<nscalars::wave>(start::wave); }

    storage_type::NColsBlockXpr<nscalars::physical>::Type physical()
    { return storage.middleCols<nscalars::physical>(start::physical); }

    storage_type::NColsBlockXpr<nscalars::implicit>::Type implicit()
    { return storage.middleCols<nscalars::implicit>(start::implicit); }

    storage_type::ConstNColsBlockXpr<nscalars::wave>::Type wave() const
    { return storage.middleCols<nscalars::wave>(start::wave); }

    storage_type::ConstNColsBlockXpr<nscalars::physical>::Type physical() const
    { return storage.middleCols<nscalars::physical>(start::physical); }

    storage_type::ConstNColsBlockXpr<nscalars::implicit>::Type implicit() const
    { return storage.middleCols<nscalars::implicit>(start::implicit); }
    /** @} */

    /** Compile-time sizes for each quantity within \c storage */
    struct size { enum {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP,,
                    SUZERAIN_PERFECT_QUANTITIES))
    }; };

#undef OP

    // Declare a named, mutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                               \
    storage_type::NColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)()                                        \
    {                                                                         \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(    \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                     \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_PERFECT_QUANTITIES)
#undef DECLARE

    // Declare a named, immutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                                    \
    storage_type::ConstNColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)() const                                       \
    {                                                                              \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(         \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                          \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_PERFECT_QUANTITIES)
#undef DECLARE

    /**
     * A foreach operation iterating over all mutable quantities in \c storage.
     * The functor is invoked as <tt>f(std::string("foo",
     * storage_type::NColsBlockXpr<size::foo>::Type))</tt> for a quantity named
     * "foo".  Each invocation must return a <tt>bool</tt> result.  See Eigen's
     * "Writing Functions Taking Eigen Types as Parameters" for suggestions on
     * how to write a functor, especially the \c const_cast hack details
     * therein.  See <tt>boost::ref</tt> for how to use a stateful functor.
     *
     * @return True if all invocations returned \c true.  False otherwise.
     */
    template <typename BinaryFunction>
    bool foreach(BinaryFunction f) {
        bool retval = true;
#define INVOKE(r, data, tuple)                                             \
        retval &= f(::std::string(                                         \
                    BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, tuple))), \
                    this->BOOST_PP_TUPLE_ELEM(2, 0, tuple)());
        BOOST_PP_SEQ_FOR_EACH(INVOKE,,SUZERAIN_PERFECT_QUANTITIES)
#undef INVOKE
        return retval;
    }

    /**
     * A foreach operation iterating over all immutable quantities in \c
     * storage.  The functor is invoked as <tt>f(std::string("foo",
     * storage_type::NColsBlockXpr<size::foo>::Type))</tt> for a quantity named
     * "foo".  Each invocation must return a <tt>bool</tt> result.  See Eigen's
     * "Writing Functions Taking Eigen Types as Parameters" for suggestions on
     * how to write a functor.  See <tt>boost::ref</tt> for how to use a
     * stateful functor.
     *
     * @return True if all invocations returned \c true.  False otherwise.
     */
    template <typename BinaryFunction>
    bool foreach(BinaryFunction f) const {
        bool retval = true;
#define INVOKE(r, data, tuple)                                             \
        retval &= f(::std::string(                                         \
                    BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, tuple))), \
                    this->BOOST_PP_TUPLE_ELEM(2, 0, tuple)());
        BOOST_PP_SEQ_FOR_EACH(INVOKE,,SUZERAIN_PERFECT_QUANTITIES)
#undef INVOKE
        return retval;
    }

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

/**
 * Using the provided state, sample the mean quantities declared in \ref
 * quantities with the notable exceptions of those listed in \ref
 * SUZERAIN_PERFECT_QUANTITIES_IMPLICIT.  This is an expensive, collective
 * method.
 *
 * @param[in]     scenario Scenario parameters.
 * @param[in]     otool    Operator definitions in use.
 * @param[in,out] swave    Destroyed in the computation
 * @param[in]     t        Current simulation time
 *
 * @return Mean quantities as B-spline coefficients.
 */
quantities sample_quantities(
        const scenario_definition &scenario,
        const operator_tools& otool,
        contiguous_state<4,complex_t> &swave,
        const real_t t);

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_QUANTITIES_HPP
