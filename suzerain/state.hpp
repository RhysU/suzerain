//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// state.hpp: Interfaces for describing/manipulating mutable state vectors
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef __SUZERAIN_STATE_HPP
#define __SUZERAIN_STATE_HPP

#include <suzerain/common.hpp>

/** @file
 * Provide interfaces for describing and manipulating mutable state vectors.
 */

namespace suzerain
{

// FIXME Document

template< typename Element >
class IStateBase
{
public:
    /**
     * Construct an instance holding the given amount of state information.
     *
     * @param variable_count Number of individual variables or
     *                       pieces of state at a given position.
     * @param vector_length  Number of positions per full state vector.
     * @param vector_count   Number of independent state vectors to store.
     */
    template< typename I1, typename I2, typename I3 >
    IStateBase(I1 variable_count, I2 vector_length, I3 vector_count)
        : variable_count(boost::numeric_cast<int>(variable_count)),
          vector_length(boost::numeric_cast<int>(vector_length)),
          vector_count(boost::numeric_cast<int>(vector_count)) {}

    /**
     * Construct an instance with the same amount of state information
     * as another.
     *
     * @param other instance to mimic in shape.
     */
    IStateBase(const IStateBase& other)
        : variable_count(other.variable_count),
          vector_length(other.vector_length),
          vector_count(other.vector_count) {}

    /** Destructor is virtual as appropriate for an abstract base class */
    virtual ~IStateBase() {}

    /**
     * Is \c this instance's shape "conformant" with <tt>other</tt>'s?
     *
     * @param other another instance to compare against.
     *
     * @return True if all of #variable_count, #vector_length, and
     *         #vector_count are identical.  False otherwise.
     */
    virtual bool isConformant(const IStateBase& other) const
    {
        return    variable_count == other.variable_count
               && vector_length  == other.vector_length
               && vector_count   == other.vector_count;
    }

    /**
     * Scale all state information by the given scale factor, i.e.
     * \f$\mbox{this} \leftarrow{} \mbox{factor}\times\mbox{this}\f$.
     *
     * @param factor Scale factor to use.
     */
    virtual void scale(const Element &factor) = 0;

    /** Number of state variables present at each position in a state vector */
    const int variable_count;

    /** Number of positions within each state vector */
    const int vector_length;

    /** Number of state vectors maintained within this instance */
    const int vector_count;
};

// FIXME Update documentation
/**
 * An interface describing a state space consisting of one or more state
 * vectors containing one or more state variables.  IState instances of the
 * same subtype and shape may be scaled by type \c FPT and added together.
 * Implementations are expected to make the underlying information available
 * through the Boost.MultiArray concept.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
 *      Boost.MultiArray</a> for more information on the MultiArray concept.
 */
template<
    typename Element,
    typename Storage,
    typename CompatibleStorage = Storage
>
class IState : public virtual IStateBase<Element>
{
    BOOST_STATIC_ASSERT(
            (Storage::dimensionality == CompatibleStorage::dimensionality));

public:
    /**
     * Construct an instance holding the given amount of state information.
     *
     * @param variable_count Number of individual variables or
     *                       pieces of state at a given position.
     * @param vector_length  Number of positions per full state vector.
     * @param vector_count   Number of independent state vectors to store.
     */
    template< typename I1, typename I2, typename I3 >
    IState(I1 variable_count, I2 vector_length, I3 vector_count)
        : IStateBase<Element>(variable_count, vector_length, vector_count) {}

    /**
     * Copy construct an instance with the same amount of state information
     * as another.
     *
     * @param other instance to mimic in shape.
     */
    template< typename T, typename U >
    IState(const IState<Element,T,U>& other)
        : IStateBase<Element>(other) {}

    /** Destructor is virtual as appropriate for an abstract base class */
    virtual ~IState() {}

    /**
     * Accumulate scaled state information by computing \f$\mbox{this}
     * \leftarrow{} \mbox{this} + \mbox{factor}\times\mbox{other}\f$.
     *
     * @param factor Scale factor to apply to <tt>other</tt>'s
     *               state information.
     * @param other Another state instance to scale and add to \c this.
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant.
     */
    virtual void addScaled(
            const Element &factor,
            const IState<Element,CompatibleStorage,Storage>& other)
            throw(std::bad_cast, std::logic_error) = 0;

    /**
     * Assign to this instance the state information from another instance
     * holding compatibly stored elements.
     *
     * @param other instance with state to be deep copied.
     *
     * @return *this
     */
    virtual void assign(
            const IState<Element,CompatibleStorage,Storage>& other)
            throw(std::bad_cast, std::logic_error) = 0;

    /**
     * Exchange <tt>this</tt>'s storage with <tt>other</tt>'s storage by moving
     * the relevant data without incurring a lot of memory overhead.  In
     * contrast with potentially optimized swap semantics, this method should
     * swap data instead of playing tricks with underlying pointers.
     *
     * @param other instance with which to exchange data.
     * @throw std::bad_cast if \c that does not have a compatible type.
     * @throw std::logic_error if \c that is not conformant.
     */
    virtual void exchange(
            IState<Element,CompatibleStorage,Storage>& other)
            throw(std::bad_cast, std::logic_error) = 0;
};

} // namespace suzerain

#endif // __SUZERAIN_STATE_HPP
