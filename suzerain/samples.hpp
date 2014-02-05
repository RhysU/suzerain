//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_SAMPLES_HPP
#define SUZERAIN_SAMPLES_HPP

/** @file
 * Named storage for sampled quantities of interest for 3D Navier--Stokes.
 */

#include <suzerain/common.hpp>
#include <suzerain/samples.h>

namespace suzerain {

/**
 * Encapsulates sample quantities like those detailed in the "Sampling
 * logistics" section of <tt>writeups/perfectgas.tex</tt>.
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
class samples
{
public:

#define NCOMPONENTS(r, data, tuple)                         \
        BOOST_PP_SEQ_SIZE(BOOST_PP_TUPLE_ELEM(2, 1, tuple))

    /* Compile-time totals of the number of scalars sampled at each point */
    struct nscalars { enum {
#define SUM(s, state, x) BOOST_PP_ADD(state, x)

        wave = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    NCOMPONENTS,,SUZERAIN_SAMPLES_WAVE)),

        physical = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    NCOMPONENTS,,SUZERAIN_SAMPLES_PHYSICAL)),

        implicit = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    NCOMPONENTS,,SUZERAIN_SAMPLES_IMPLICIT)),

        total = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0, BOOST_PP_SEQ_TRANSFORM(
                    NCOMPONENTS,,SUZERAIN_SAMPLES))

#undef SUM
    }; };

    /** Simulation time when samples were obtained */
    real_t t;

    /** Type of the contiguous storage used to house all scalars */
    typedef Array<real_t, Dynamic, nscalars::total, ColMajor> storage_type;

    /** Contiguous storage used to house all means */
    storage_type storage;

    /**
     * Constructor setting <tt>this->t = NaN</tt>.
     * Caller will need to resize <tt>this->storage</tt> prior to use.
     */
    samples();

    /**
     * Constructor setting <tt>this->t = t</tt>.
     * Caller will need to resize <tt>this->storage</tt> prior to use.
     */
    explicit samples(real_t t);

    /**
     * Constructor setting <tt>this->t = t</tt> and preparing a zero-filled \c
     * storage containing \c Ny rows.
     */
    samples(real_t t, storage_type::Index Ny);

#define ASSIGN(r, data, tuple)                                              \
        BOOST_PP_TUPLE_ELEM(2, 0, tuple) = BOOST_PP_TUPLE_ELEM(2, 1, tuple)

#define QUANTITIES_NCOMPONENTS(r, data, tuple)                                          \
        (BOOST_PP_TUPLE_ELEM(2, 0, tuple) BOOST_PP_COMMA() NCOMPONENTS(r, data, tuple))

    /** Compile-time offsets for each quantity within \c storage */
    struct start { enum {
        wave     = 0,                              // Start of wave block
        physical = wave     + nscalars::wave,      // Start of physical block
        implicit = physical + nscalars::physical,  // Start of implicit block
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(
                ASSIGN,,SUZERAIN_SHIFTED_SUM(BOOST_PP_SEQ_TRANSFORM(
                        QUANTITIES_NCOMPONENTS,,SUZERAIN_SAMPLES))))
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
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(ASSIGN,, BOOST_PP_SEQ_TRANSFORM(
                QUANTITIES_NCOMPONENTS,,SUZERAIN_SAMPLES)))
    }; };

#undef QUANTITIES_NCOMPONENTS

#undef ASSIGN

#undef NCOMPONENTS

    // Declare a named, mutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                               \
    storage_type::NColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)()                                        \
    {                                                                         \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(    \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                     \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_SAMPLES)
#undef DECLARE

    // Declare a named, immutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                                    \
    storage_type::ConstNColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)() const                                       \
    {                                                                              \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(         \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                          \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,SUZERAIN_SAMPLES)
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
        BOOST_PP_SEQ_FOR_EACH(INVOKE,,SUZERAIN_SAMPLES)
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
        BOOST_PP_SEQ_FOR_EACH(INVOKE,,SUZERAIN_SAMPLES)
#undef INVOKE
        return retval;
    }

};

} // end namespace suzerain

#endif // SUZERAIN_SAMPLES_HPP
