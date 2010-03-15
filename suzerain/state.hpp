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
#include <suzerain/blas_et_al.hpp>

/** @file
 * Provides an interface and implementations for an abstract state concept.
 */

namespace suzerain
{

/**
 * An interface describing a state space consisting of one or more state
 * vectors containing one or more state variables.  IState instances of the
 * same subtype and shape may be scaled by type \c FPT and added together.
 * Implementations are expected to make the underlying information available
 * through the Boost.MultiArray concept.
 *
 * When the template parameter \c Interleaved is true, state variables are
 * interleaved with one another.  That is, variables are the fastest index.
 * When \c Interleaved is false, each state variable is stored in a contiguous
 * block of memory.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
 *      Boost.MultiArray</a> for more information on the MultiArray concept.
 */
template< typename FPT, bool Interleaved = true >
class IState
{
public:

    /**
     * Are state variables interleaved with one another?  When true, variables
     * are the fastest index.  When false, each state variable is stored in a
     * separate, contiguous block of memory.
     */
    static const bool interleaved = Interleaved;

    /**
     * A signed integral type used to index into MultiArrays.
     * Provided as a convenience typedef from boost::multi_array::types.
     **/
    typedef boost::multi_array_types::index index;

    /**
     * An unsigned integral type used to express a MultiArray's shape.
     * Provided as a convenience typedef from boost::multi_array::types.
     **/
    typedef boost::multi_array_types::size_type size_type;

    /**
     * A signed integral type used to express index offsets for MultiArrays.
     * Provided as a convenience typedef from boost::multi_array::types.
     **/
    typedef boost::multi_array_types::difference_type difference_type;

    /** Number of state variables present at each position in a state vector */
    const size_type variable_count;

    /** Number of positions within each state vector */
    const size_type vector_length;

    /** Number of state vectors maintained within this instance */
    const size_type vector_count;

    /**
     * Construct an instance holding the given amount of state information.
     *
     * @param variable_count Number of individual variables or
     *                       pieces of state at a given position.
     * @param vector_length  Number of positions per full state vector.
     * @param vector_count   Number of independent state vectors to store.
     */
    IState(size_type variable_count,
           size_type vector_length,
           size_type vector_count)
        :   variable_count(variable_count),
            vector_length(vector_length),
            vector_count(vector_count) {}

    /**
     * Construct an instance with the same amount of state information
     * as another.
     *
     * @param other instance to mimic in shape.
     */
    IState(const IState& other)
        : variable_count(other.variable_count),
          vector_length(other.vector_length),
          vector_count(other.vector_count) {}

    /**
     * Construct an instance with the same amount of state information
     * as another.  The other instance varies in how it is interleaved.
     *
     * @param other instance to mimic in shape.
     */
    IState(const IState<FPT,!Interleaved>& other)
        : variable_count(other.variable_count),
          vector_length(other.vector_length),
          vector_count(other.vector_count) {}

    /** Virtual destructor appropriate for an abstract base class */
    virtual ~IState() {}

    /**
     * Is \c this instance's shape "conformant" with <tt>other</tt>'s?
     * Interleaving of state variables does not influence the comparison,
     * but subclasses are free to change this behavior.
     *
     * @param other another instance to compare against.
     *
     * @return True if all of #variable_count, #vector_length, and
     *         #vector_count are identical.  False otherwise.
     */
    virtual bool isConformant(const IState& other) const
    {
        return isConformantHelper(other);
    }

    /**
     * Is \c this instance's shape "conformant" with <tt>other</tt>'s?
     * Interleaving of state variables does not influence the comparison,
     * but subclasses are free to change this behavior.
     *
     * @param other another instance to compare against.
     *
     * @return True if all of #variable_count, #vector_length, and
     *         #vector_count are identical.  False otherwise.
     */
    virtual bool isConformant(const IState<FPT,!Interleaved>& other) const
    {
        return isConformantHelper(other);
    }

    /**
     * Scale all state information by the given scale factor, i.e.
     * \f$\mbox{this} \leftarrow{} \mbox{factor}\times\mbox{this}\f$.
     *
     * @param factor Scale factor to use.
     */
    virtual void scale(const FPT factor) = 0;

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
    virtual void addScaled(const FPT factor,
                           const IState<FPT,Interleaved>& other)
                           throw(std::bad_cast, std::logic_error) = 0;

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
    virtual void addScaled(const FPT factor,
                           const IState<FPT,!Interleaved>& other)
                           throw(std::bad_cast, std::logic_error) = 0;

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     *
     * @return *this
     */
    virtual IState& operator=(const IState<FPT,Interleaved>& that)
                              throw(std::bad_cast, std::logic_error) = 0;

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     *
     * @return *this
     */
    virtual IState& operator=(const IState<FPT,!Interleaved>& that)
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
    virtual void exchange(IState<FPT,Interleaved>& other)
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
    virtual void exchange(IState<FPT,!Interleaved>& other)
                          throw(std::bad_cast, std::logic_error) = 0;

    // TODO Add IState<FPT>::swap(IState<FPT>&) to public interface

private:

    /**
     * Is \c this instance's shape equivalent to <tt>other</tt>'s?
     * Interleaving of state variables does not influence the comparison.
     *
     * @param other another instance to compare against.
     *
     * @return True if all of #variable_count, #vector_length, and
     *         #vector_count are identical.  False otherwise.
     */
    template<bool OtherInterleaved>
    bool isConformantHelper(const IState<FPT,OtherInterleaved>& other) const {
        return    variable_count == other.variable_count
               && vector_length  == other.vector_length
               && vector_count   == other.vector_count;
    }
};

/**
 * An implementation of IState<FPT> for real-valued state information.
 */
template< typename FPT, bool Interleaved = true >
class RealState : public IState<FPT,Interleaved>
{
public:
    /**
     * Are state variables interleaved with one another?  When true, variables
     * are the fastest index.  When false, each state variable is stored in a
     * separate, contiguous block of memory.
     */
    static const bool interleaved = Interleaved;

    /**
     * A signed integral type used to index into MultiArrays.
     * Provided as a convenience typedef from IState<FPT,Interleaved>::index.
     **/
    typedef typename IState<FPT,Interleaved>::index index;

    /**
     * An unsigned integral type used to express a MultiArray's shape.
     * Provided as a convenience typedef from IState<FPT,Interleaved>::size_type.
     **/
    typedef typename IState<FPT,Interleaved>::size_type size_type;

    /**
     * A signed integral type used to express index offsets for MultiArrays.
     * Provided as a convenience typedef from
     * IState<FPT,Interleaved>::difference_type.
     **/
    typedef typename IState<FPT,Interleaved>::difference_type difference_type;

    /**
     * Construct an instance holding the given amount of real-valued state
     * information.
     *
     * @param variable_count Number of individual variables or
     *                       pieces of state at a given position.
     * @param vector_length  Number of positions per full state vector.
     * @param vector_count   Number of independent state vectors to store.
     * @throw std::bad_alloc if a memory allocation error occurs.
     */
    RealState(size_type variable_count,
              size_type vector_length,
              size_type vector_count)
              throw(std::bad_alloc);

    /**
     * Construct an instance with the same amount of state information
     * as another.  The copy constructor is explicit because instances maintain
     * non-trivial amounts of information; unintentional copying would be
     * problematic.
     *
     * @param other instance to mimic in shape.
     * @throw std::bad_alloc if a memory allocation error occurs.
     */
    explicit RealState(const RealState& other)
                       throw(std::bad_alloc);

    /**
     * Construct an instance with the same amount of state information
     * as another.  The copy constructor is explicit because instances maintain
     * non-trivial amounts of information; unintentional copying would be
     * problematic.
     *
     * @param other instance to mimic in shape.
     * @throw std::bad_alloc if a memory allocation error occurs.
     */
    explicit RealState(const RealState<FPT,!Interleaved>& other)
                       throw(std::bad_alloc);

    /**
     * Scale all state information by the given scale factor, i.e.
     * \f$\mbox{this} \leftarrow{} \mbox{factor}\times\mbox{this}\f$.
     *
     * @param factor Scale factor to use.
     */
    virtual void scale(const FPT factor);

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
    virtual void addScaled(const FPT factor,
                           const IState<FPT,Interleaved>& other)
                           throw(std::bad_cast, std::logic_error);

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
    virtual void addScaled(const FPT factor,
                           const IState<FPT,!Interleaved>& other)
                           throw(std::bad_cast, std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::logic_error if \c other is not conformant.
     */
    RealState& operator=(const RealState& that)
                         throw(std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::logic_error if \c other is not conformant.
     */
    RealState& operator=(const RealState<FPT,!Interleaved>& that)
                         throw(std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant.
     */
    virtual RealState& operator=(const IState<FPT,Interleaved>& that)
                                 throw(std::bad_cast, std::logic_error)
    {
        return operator=(dynamic_cast<const RealState&>(that));
    }

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant.
     */
    virtual RealState& operator=(const IState<FPT,!Interleaved>& that)
                                 throw(std::bad_cast, std::logic_error)
    {
        return operator=(dynamic_cast<const RealState<FPT,!Interleaved>&>(that));
    }

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
    virtual void exchange(IState<FPT,Interleaved>& other)
                          throw(std::bad_cast, std::logic_error);

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
    virtual void exchange(IState<FPT,!Interleaved>& other)
                          throw(std::bad_cast, std::logic_error);

    // TODO Implement swap using https://svn.boost.org/trac/boost/ticket/1045

protected:
    /**
     * Provides MultiArray storage order information for the template
     * instantiation.
     **/
    static boost::general_storage_order<3> storage_order();

public:
    /** Type of the Boost.MultiArray exposed through RealState#data. */
    typedef typename boost::multi_array< FPT, 3,
            typename ::suzerain::blas::allocator<FPT>::type > state_type;

    /**
     * A three-dimensional Boost.MultiArray of real-valued state kept
     * zero-indexed in column-major order.  The indices are over
     * #variable_count, #vector_length, and #vector_count respectively.
     * The information is contiguous.
     */
    state_type data;
};

template< typename FPT, bool Interleaved >
boost::general_storage_order<3> RealState<FPT,Interleaved>::storage_order()
{
    if (Interleaved) {
        const size_type ordering[3]  = { 0, 1, 2 }; // Real variables fastest
        const bool      ascending[3] = { true, true, true };
        boost::general_storage_order<3> result(ordering, ascending);
        return result;
    } else {
        const size_type ordering[3]  = { 2, 0, 1 }; // Real variables slowest
        const bool      ascending[3] = { true, true, true };
        boost::general_storage_order<3> result(ordering, ascending);
        return result;
    }
}

template< typename FPT, bool Interleaved >
RealState<FPT,Interleaved>::RealState(
        typename IState<FPT,Interleaved>::size_type variable_count,
        typename IState<FPT,Interleaved>::size_type vector_length,
        typename IState<FPT,Interleaved>::size_type vector_count)
throw(std::bad_alloc)
    : IState<FPT,Interleaved>(variable_count, vector_length, vector_count),
      data(boost::extents[variable_count][vector_length][vector_count],
           RealState<FPT,Interleaved>::storage_order())
{
    // NOP
}

template< typename FPT, bool Interleaved >
RealState<FPT,Interleaved>::RealState(
        const RealState& other)
throw(std::bad_alloc)
    : IState<FPT,Interleaved>(other),
      data(boost::extents[other.variable_count]
                         [other.vector_length]
                         [other.vector_count],
           RealState<FPT,Interleaved>::storage_order())
{
    std::memcpy(data.data(), other.data.data(),
                sizeof(state_type::element)*data.num_elements());
}

template< typename FPT, bool Interleaved >
RealState<FPT,Interleaved>::RealState(
        const RealState<FPT,!Interleaved>& other)
throw(std::bad_alloc)
    : IState<FPT,Interleaved>(other),
      data(boost::extents[other.variable_count]
                         [other.vector_length]
                         [other.vector_count],
           RealState<FPT,Interleaved>::storage_order())
{
    data = other.data; // TODO Ensure deep assignment has adequate performance
}

template< typename FPT, bool Interleaved >
void RealState<FPT,Interleaved>::scale(
        const FPT factor)
{
    if (factor == FPT(0)) {
        memset(data.data(), 0,
               data.num_elements()*sizeof(state_type::element));
    } else if (factor == FPT(1)) {
        // NOP
    } else {
        suzerain::blas::scal(data.num_elements(), factor, data.data(), 1);
    }
}

template< typename FPT, bool Interleaved >
void RealState<FPT,Interleaved>::addScaled(
        const FPT factor,
        const IState<FPT,Interleaved>& other)
throw(std::bad_cast,
      std::logic_error)
{
    if (!isConformant(other))
        throw std::logic_error("Nonconformant other in addScaled");

    const RealState& o = dynamic_cast<const RealState&>(other);

    suzerain::blas::axpy(data.num_elements(), factor,
                         o.data.data(), 1,
                         data.data(), 1);
}

template< typename FPT, bool Interleaved >
void RealState<FPT,Interleaved>::addScaled(
        const FPT factor,
        const IState<FPT,!Interleaved>& other)
throw(std::bad_cast,
      std::logic_error)
{
    if (!isConformant(other))
        throw std::logic_error("Nonconformant other in addScaled");

    const RealState<FPT,!Interleaved>& o
        = dynamic_cast<const RealState<FPT,!Interleaved>&>(other);

    for (int i = 0; i < this->variable_count; ++i) {
        suzerain::blas::axpy(this->vector_length*this->vector_count, factor,
                             &o.data[i][0][0], o.data.strides()[1],
                             &data[i][0][0], data.strides()[1]);
    }
}

template< typename FPT, bool Interleaved >
RealState<FPT,Interleaved>& RealState<FPT,Interleaved>::operator=(
        const RealState& that)
throw(std::logic_error)
{
    if (this != &that) {
        if (!isConformant(that))
            throw std::logic_error("Nonconformant that in operator=");

        std::memcpy(data.data(), that.data.data(),
                    sizeof(state_type::element)*data.num_elements());
    }

    return *this;
}

template< typename FPT, bool Interleaved >
RealState<FPT,Interleaved>& RealState<FPT,Interleaved>::operator=(
        const RealState<FPT,!Interleaved>& that)
throw(std::logic_error)
{
    if (!isConformant(that))
        throw std::logic_error("Nonconformant that in operator=");

    this->data = that.data;

    return *this;
}

template< typename FPT, bool Interleaved >
void RealState<FPT,Interleaved>::exchange(
        IState<FPT,Interleaved>& other)
throw(std::bad_cast,
      std::logic_error)
{
    if (this != &other) {
        if (!isConformant(other))
            throw std::logic_error("Nonconformant other in exchange");

        RealState& o = dynamic_cast<RealState&>(other);

        suzerain::blas::swap(data.num_elements(),
                             o.data.data(), 1,
                             data.data(), 1);
    }
}

template< typename FPT, bool Interleaved >
void RealState<FPT,Interleaved>::exchange(
        IState<FPT,!Interleaved>& other)
throw(std::bad_cast,
      std::logic_error)
{
    if (!isConformant(other))
        throw std::logic_error("Nonconformant other in exchange");

    RealState<FPT,!Interleaved>& o
        = dynamic_cast<RealState<FPT,!Interleaved>&>(other);

    for (int i = 0; i < this->variable_count; ++i) {
        suzerain::blas::swap(this->vector_length*this->vector_count,
                             &o.data[i][0][0], o.data.strides()[1],
                             &data[i][0][0], data.strides()[1]);
    }
}

/**
 * An implementation of IState<FPT> for complex-valued state information.
 * \c FPT is the real scalar type used underneath the complex type.
 */
template< typename FPT, bool Interleaved = true >
class ComplexState : public IState<FPT,Interleaved>
{
public:
    /**
     * Are state variables interleaved with one another?  When true, variables
     * are the fastest index.  When false, each state variable is stored in a
     * separate, contiguous block of memory.
     */
    static const bool interleaved = Interleaved;

    /**
     * A signed integral type used to index into MultiArrays.
     * Provided as a convenience typedef from IState<FPT,Interleaved>::index.
     **/
    typedef typename IState<FPT,Interleaved>::index index;

    /**
     * An unsigned integral type used to express a MultiArray's shape.
     * Provided as a convenience typedef from IState<FPT,Interleaved>::size_type.
     **/
    typedef typename IState<FPT,Interleaved>::size_type size_type;

    /**
     * A signed integral type used to express index offsets for MultiArrays.
     * Provided as a convenience typedef from
     * IState<FPT,Interleaved>::difference_type.
     **/
    typedef typename IState<FPT,Interleaved>::difference_type difference_type;

    /**
     * Construct an instance holding the given amount of complex-valued state
     * information.
     *
     * @param variable_count Number of individual variables or
     *                       pieces of state at a given position.
     * @param vector_length  Number of positions per full state vector.
     * @param vector_count   Number of independent state vectors to store.
     * @throw std::bad_alloc if a memory allocation error occurs
     */
    ComplexState(size_type variable_count,
                 size_type vector_length,
                 size_type vector_count)
                 throw(std::bad_alloc);

    /**
     * Construct an instance with the same amount of state information as
     * another.  The copy constructor is explicit because instances maintain
     * non-trivial amounts of information; unintentional copying would be
     * problematic.
     *
     * @param other instance to mimic in shape.
     * @throw std::bad_alloc if a memory allocation error occurs.
     */
    explicit ComplexState(const ComplexState& other)
                          throw(std::bad_alloc);

    /**
     * Construct an instance with the same amount of state information as
     * another.  The copy constructor is explicit because instances maintain
     * non-trivial amounts of information; unintentional copying would be
     * problematic.
     *
     * @param other instance to mimic in shape.
     * @throw std::bad_alloc if a memory allocation error occurs.
     */
    explicit ComplexState(const ComplexState<FPT,!Interleaved>& other)
                          throw(std::bad_alloc);

    /**
     * Scale all state information by the given scale factor, i.e.
     * \f$\mbox{this} \leftarrow{} \mbox{factor}\times\mbox{this}\f$.
     *
     * @param factor Scale factor to use.
     */
    virtual void scale(const FPT factor);

    /**
     * Accumulate scaled state information by computing \f$\mbox{this}
     * \leftarrow{} \mbox{this} + \mbox{factor}\times\mbox{other}\f$.
     *
     * @param factor Scale factor to apply to <tt>other</tt>'s
     *               state information.
     * @param other Another state instance to scale and add to \c this.
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant in shape.
     */
    virtual void addScaled(const FPT factor,
                           const IState<FPT,Interleaved>& other)
                           throw(std::bad_cast, std::logic_error);

    /**
     * Accumulate scaled state information by computing \f$\mbox{this}
     * \leftarrow{} \mbox{this} + \mbox{factor}\times\mbox{other}\f$.
     *
     * @param factor Scale factor to apply to <tt>other</tt>'s
     *               state information.
     * @param other Another state instance to scale and add to \c this.
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant in shape.
     */
    virtual void addScaled(const FPT factor,
                           const IState<FPT,!Interleaved>& other)
                           throw(std::bad_cast, std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::logic_error if \c other is not conformant.
     */
    ComplexState& operator=(const ComplexState& that)
                            throw(std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::logic_error if \c other is not conformant.
     */
    ComplexState& operator=(const ComplexState<FPT,!Interleaved>& that)
                            throw(std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant.
     */
    virtual ComplexState& operator=(const IState<FPT,Interleaved>& that)
                                    throw(std::bad_cast, std::logic_error)
    {
        return operator=( dynamic_cast<const ComplexState&>(that));
    }

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant.
     */
    virtual ComplexState& operator=(const IState<FPT,!Interleaved>& that)
                                    throw(std::bad_cast, std::logic_error)
    {
        return operator=(
                dynamic_cast<const ComplexState<FPT,!Interleaved>&>(that));
    }

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
    virtual void exchange(IState<FPT,Interleaved>& other)
                          throw(std::bad_cast, std::logic_error);

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
    virtual void exchange(IState<FPT,!Interleaved>& other)
                          throw(std::bad_cast, std::logic_error);

    // TODO Implement swap using https://svn.boost.org/trac/boost/ticket/1045

protected:
    /**
     * Provides MultiArray complex component storage order information for the
     * template instantiation.
     **/
    static boost::general_storage_order<4> components_storage_order();

    /**
     * Provides MultiArray complex-valued storage order information for the
     * template instantiation.
     **/
    static boost::general_storage_order<3> storage_order();

public:
    /** Type of the Boost.MultiArray exposed through #components. */
    typedef typename boost::multi_array< FPT, 4,
            typename ::suzerain::blas::allocator<FPT>::type > components_type;

    /** Type of the Boost.MultiArray exposed through #real and #imag. */
    typedef typename boost::array_view_gen<
        components_type,3>::type component_type;

    /** Type of the Boost.MultiArray exposed through #data. */
    typedef typename boost::multi_array_ref<std::complex<FPT>, 3> state_type;

    /**
     * A four-dimensional Boost.MultiArray of real-valued state kept
     * zero-indexed in column-major order.  The indices are over the real and
     * imaginary component, #variable_count, #vector_length, and
     * #vector_count respectively.  The information is contiguous.
     */
    components_type components;

    /**
     * A three-dimensional Boost.MultiArray view of the real portion of the
     * complex-valued state zero-indexed in column-major order.  The indices
     * are over #variable_count, #vector_length, and #vector_count
     * respectively.  The data is not contiguous in memory.
     */
    component_type real;

    /**
     * A three-dimensional Boost.MultiArray view of the imaginary portion of
     * the complex-valued state zero-indexed in column-major order.  The
     * indices are over #variable_count, #vector_length, and #vector_count
     * respectively.  The data is not contiguous in memory.
     */
    component_type imag;

    /**
     * A three-dimensional Boost.MultiArray view of complex-valued state kept
     * zero-indexed in column-major order.  The indices are over
     * #variable_count, #vector_length, and #vector_count respectively.
     * The information is contiguous.
     */
    state_type data;
};

template< typename FPT, bool Interleaved >
boost::general_storage_order<4>
ComplexState<FPT,Interleaved>::components_storage_order()
{
    if (Interleaved) {
        const size_type ordering[4]  = { 0, 1, 2, 3 }; // Complex variables fastest
        const bool      ascending[4] = { true, true, true, true };
        boost::general_storage_order<4> result(ordering, ascending);
        return result;
    } else {
        const size_type ordering[4]  = { 0, 3, 1, 2 }; // Complex variables slowest
        const bool      ascending[4] = { true, true, true };
        boost::general_storage_order<4> result(ordering, ascending);
        return result;
    }
}

template< typename FPT, bool Interleaved >
boost::general_storage_order<3>
ComplexState<FPT,Interleaved>::storage_order()
{
    if (Interleaved) {
        const size_type ordering[3]  = { 0, 1, 2 }; // Complex variables fastest
        const bool      ascending[3] = { true, true, true };
        boost::general_storage_order<3> result(ordering, ascending);
        return result;
    } else {
        const size_type ordering[3]  = { 2, 0, 1 }; // Complex variables slowest
        const bool      ascending[3] = { true, true, true };
        boost::general_storage_order<3> result(ordering, ascending);
        return result;
    }
}

template< typename FPT, bool Interleaved >
ComplexState<FPT,Interleaved>::ComplexState(
        typename IState<FPT,Interleaved>::size_type variable_count,
        typename IState<FPT,Interleaved>::size_type vector_length,
        typename IState<FPT,Interleaved>::size_type vector_count)
throw(std::bad_alloc)
    : IState<FPT,Interleaved>(variable_count, vector_length, vector_count),
      components(
           boost::extents[2][variable_count][vector_length][vector_count],
           ComplexState<FPT,Interleaved>::components_storage_order()),
      real(components[boost::indices[0]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      imag(components[boost::indices[1]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      data(reinterpret_cast<std::complex<FPT> *>(components.data()),
           boost::extents[variable_count][vector_length][vector_count],
           ComplexState<FPT,Interleaved>::storage_order())
{
    // NOP
}

template< typename FPT, bool Interleaved >
ComplexState<FPT,Interleaved>::ComplexState(
        const ComplexState& other)
throw(std::bad_alloc)
    : IState<FPT,Interleaved>(other),
      components(boost::extents[2]
                               [other.variable_count]
                               [other.vector_length]
                               [other.vector_count],
                 ComplexState<FPT,Interleaved>::components_storage_order()),
      real(components[boost::indices[0]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      imag(components[boost::indices[1]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      data(reinterpret_cast<std::complex<FPT> *>(components.data()),
           boost::extents[other.variable_count]
                         [other.vector_length]
                         [other.vector_count],
           ComplexState<FPT,Interleaved>::storage_order())
{
    std::memcpy(data.data(), other.data.data(),
                sizeof(state_type::element)*data.num_elements());
}


template< typename FPT, bool Interleaved >
ComplexState<FPT,Interleaved>::ComplexState(
        const ComplexState<FPT,!Interleaved>& other)
throw(std::bad_alloc)
    : IState<FPT,Interleaved>(other),
      components(boost::extents[2]
                               [other.variable_count]
                               [other.vector_length]
                               [other.vector_count],
                 ComplexState<FPT,Interleaved>::components_storage_order()),
      real(components[boost::indices[0]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      imag(components[boost::indices[1]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      data(reinterpret_cast<std::complex<FPT> *>(components.data()),
           boost::extents[other.variable_count]
                         [other.vector_length]
                         [other.vector_count],
           ComplexState<FPT,Interleaved>::storage_order())
{
    components = other.components;
}

template< typename FPT, bool Interleaved >
void ComplexState<FPT,Interleaved>::scale(const FPT factor)
{
    if (factor == FPT(0)) {
        memset(data.data(), 0,
               data.num_elements()*sizeof(state_type::element));
    } else if (factor == FPT(1)) {
        // NOP
    } else {
        suzerain::blas::scal(
                components.num_elements(), factor, components.data(), 1);
    }
}

template< typename FPT, bool Interleaved >
void ComplexState<FPT,Interleaved>::addScaled(
        const FPT factor,
        const IState<FPT,Interleaved>& other)
throw(std::bad_cast,
      std::logic_error)
{
    if (!isConformant(other))
        throw std::logic_error("Nonconformant other in addScaled");

    const ComplexState<FPT,Interleaved>& o
        = dynamic_cast<const ComplexState<FPT,Interleaved>&>(other);

    suzerain::blas::axpy(components.num_elements(),
                         factor, o.components.data(), 1,
                         components.data(), 1);
}

template< typename FPT, bool Interleaved >
void ComplexState<FPT,Interleaved>::addScaled(
        const FPT factor,
        const IState<FPT,!Interleaved>& other)
throw(std::bad_cast,
      std::logic_error)
{
    if (!isConformant(other))
        throw std::logic_error("Nonconformant other in addScaled");

    const ComplexState<FPT,!Interleaved>& o
        = dynamic_cast<const ComplexState<FPT,!Interleaved>&>(other);

    for (int i = 0; i < this->variable_count; ++i) {
        suzerain::blas::axpy(this->vector_length*this->vector_count, factor,
                             &o.data[i][0][0], o.data.strides()[1],
                             &data[i][0][0], data.strides()[1]);
    }
}

template< typename FPT, bool Interleaved >
ComplexState<FPT,Interleaved>& ComplexState<FPT,Interleaved>::operator=(
        const ComplexState& that)
throw(std::logic_error)
{
    if (this != &that) {
        if (!isConformant(that))
            throw std::logic_error("Nonconformant that in operator=");

        std::memcpy(data.data(), that.data.data(),
                    sizeof(state_type::element)*data.num_elements());
    }

    return *this;
}

template< typename FPT, bool Interleaved >
ComplexState<FPT,Interleaved>& ComplexState<FPT,Interleaved>::operator=(
        const ComplexState<FPT,!Interleaved>& that)
throw(std::logic_error)
{
    if (!isConformant(that))
        throw std::logic_error("Nonconformant that in operator=");

    this->data = that.data;

    return *this;
}

template< typename FPT, bool Interleaved >
void ComplexState<FPT,Interleaved>::exchange(
        IState<FPT,Interleaved>& other)
throw(std::bad_cast,
      std::logic_error)
{
    if (this != &other) {
        if (!isConformant(other))
            throw std::logic_error("Nonconformant other in exchange");

        ComplexState<FPT,Interleaved>& o
            = dynamic_cast<ComplexState<FPT,Interleaved>&>(other);

        suzerain::blas::swap(components.num_elements(),
                             o.components.data(), 1,
                             components.data(), 1);
    }
}

template< typename FPT, bool Interleaved >
void ComplexState<FPT,Interleaved>::exchange(
        IState<FPT,!Interleaved>& other)
throw(std::bad_cast,
      std::logic_error)
{
    if (!isConformant(other))
        throw std::logic_error("Nonconformant other in exchange");

    ComplexState<FPT,!Interleaved>& o
        = dynamic_cast<ComplexState<FPT,!Interleaved>&>(other);

    for (int i = 0; i < this->variable_count; ++i) {
        suzerain::blas::swap(this->vector_length*this->vector_count,
                             &o.data[i][0][0], o.data.strides()[1],
                             &data[i][0][0], data.strides()[1]);
    }
}

} // namespace suzerain

#endif // __SUZERAIN_STATE_HPP
