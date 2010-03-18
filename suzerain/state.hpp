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
// state.hpp: Class to manage mutable state vectors
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef __SUZERAIN_STATE_HPP
#define __SUZERAIN_STATE_HPP

#include <suzerain/common.hpp>

/** @file
 * Provides an interface for an abstract state vector collection.
 */

namespace suzerain
{

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
template< typename Storage, typename CompatibleStorage = Storage >
class IState
{
    BOOST_STATIC_ASSERT((Storage::dimensionality == 3));
    BOOST_STATIC_ASSERT((CompatibleStorage::dimensionality == 3));
    BOOST_STATIC_ASSERT((boost::is_same<typename Storage::element_type, typename CompatibleStorage::element_type>::value));

public:
    typedef typename Storage::element_type element_type;

    /** Number of state variables present at each position in a state vector */
    const int variable_count;

    /** Number of positions within each state vector */
    const int vector_length;

    /** Number of state vectors maintained within this instance */
    const int vector_count;

    /**
     * Construct an instance holding the given amount of state information.
     *
     * @param variable_count Number of individual variables or
     *                       pieces of state at a given position.
     * @param vector_length  Number of positions per full state vector.
     * @param vector_count   Number of independent state vectors to store.
     */
    template< typename Integer1, typename Integer2, typename Integer3 >
    IState(Integer1 variable_count,
           Integer2 vector_length,
           Integer3 vector_count)
        :   variable_count(boost::numeric_cast<int>(variable_count)),
            vector_length(boost::numeric_cast<int>(vector_length)),
            vector_count(boost::numeric_cast<int>(vector_count)) {}

    /**
     * Construct an instance with the same amount of state information
     * as another.
     *
     * @param istate instance to mimic in shape.
     */
    template< typename T, typename U >
    IState(const IState<T, U>& istate)
        : variable_count(istate.variable_count),
          vector_length(istate.vector_length),
          vector_count(istate.vector_count) {}

    /** Destructor is virtual as appropriate for an abstract base class */
    virtual ~IState() {}

    /**
     * Is \c this instance's shape "conformant" with <tt>other</tt>'s?
     *
     * @param other another instance to compare against.
     *
     * @return True if all of #variable_count, #vector_length, and
     *         #vector_count are identical.  False otherwise.
     */
    virtual bool isConformant(
            const IState<CompatibleStorage,Storage>& other) const
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
    virtual void scale(const element_type &factor) = 0;

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
    virtual void addScaled(const element_type &factor,
                           const IState<CompatibleStorage,Storage>& other)
                           throw(std::bad_cast, std::logic_error) = 0;

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param other instance to copy.
     *
     * @return *this
     */
    virtual void copy(const IState<CompatibleStorage,Storage>& other)
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
    virtual void exchange(IState<CompatibleStorage,Storage>& other)
                          throw(std::bad_cast, std::logic_error) = 0;
};

} // namespace suzerain

#endif // __SUZERAIN_STATE_HPP
