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
// state_impl.hpp: Implementations of the IState interface
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef __SUZERAIN_STATE_IMPL_HPP
#define __SUZERAIN_STATE_IMPL_HPP

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/storage.hpp>
#include <suzerain/state.hpp>

/** @file
 * Provides implementations of the IState interface.
 */

namespace suzerain
{

template<
    typename Element,
    typename Allocator = std::allocator<Element>
>
class InterleavedState
    : public boost::noncopyable,
      public IState<typename suzerain::storage::Interleaved<Element> >
{
public:
    typedef typename suzerain::storage::Interleaved<Element>
        interleaved_storage;

    template< typename Integer1,
              typename Integer2,
              typename Integer3,
              typename Integer4 >
    InterleavedState(Integer1 variable_count,
                     Integer2 vector_length,
                     Integer3 vector_count,
                     Integer4 min_contiguous_block = 0,
                     const Allocator &allocator = Allocator() );

    virtual ~InterleavedState();

    virtual void scale(const Element &factor);

    virtual void addScaled(
            const Element &factor,
            const IState<interleaved_storage>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void copy(
            const IState<interleaved_storage>& other)
            throw(std::bad_cast, std::logic_error);

    virtual void exchange(
            IState<interleaved_storage>& other)
            throw(std::bad_cast, std::logic_error);

private:
    typedef typename boost::multi_array_ref<
                Element, interleaved_storage::dimensionality> root_type_;

    Allocator allocator_;
    typename Allocator::size_type block_size_;
    typename Allocator::pointer raw_;
    root_type_ root_;
};

template< typename Element, typename Allocator >
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4 >
InterleavedState<Element,Allocator>::InterleavedState(
        Integer1 variable_count,
        Integer2 vector_length,
        Integer3 vector_count,
        Integer4 min_contiguous_block,
        const Allocator &allocator)
    : IState<interleaved_storage>(variable_count, vector_length, vector_count),
      allocator_(allocator),
      block_size_(std::max(
           variable_count*vector_length*vector_count, min_contiguous_block)),
      raw_(allocator_.allocate(block_size_)),
      root_(raw_,
            boost::extents[variable_count][variable_length][vector_count],
            interleaved_storage::element_storage_order())
{
    // NOP
}

template< typename Element, typename Allocator >
InterleavedState<Element,Allocator>::~InterleavedState()
{
    allocator_.deallocate(raw_, block_size_);
}


// FIXME
#ifdef FIXME_BLOCK_DISABLED

/**
 * An implementation of IState<FPT> for real-valued state information.
 */
template< typename FPT, bool Interleaved = true >
class RealState : public IState<FPT,Interleaved>
{
public:

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
    const bool ascending[3] = { true, true, true };
    if (Interleaved) {
        const size_type ordering[3] = { 0, 1, 2 }; // Real variables fastest
        boost::general_storage_order<3> result(ordering, ascending);
        return result;
    } else {
        const size_type ordering[3] = { 1, 2, 0 }; // Real variables slowest
        boost::general_storage_order<3> result(ordering, ascending);
        return result;
    }
}

template< typename FPT, bool Interleaved >
RealState<FPT,Interleaved>::RealState(
        size_type variable_count,
        size_type vector_length,
        size_type vector_count)
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
    suzerain::blas::copy(
            data.num_elements(), other.data.data(), 1, data.data(), 1);
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
throw(std::bad_cast, std::logic_error)
{
    if (!isConformant(other)) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

    const RealState& o = dynamic_cast<const RealState&>(other);

    suzerain::blas::axpy(data.num_elements(), factor,
                         o.data.data(), 1,
                         data.data(), 1);
}

template< typename FPT, bool Interleaved >
void RealState<FPT,Interleaved>::addScaled(
        const FPT factor,
        const IState<FPT,!Interleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (!isConformant(other)) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

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
        if (!isConformant(that)) throw std::logic_error(
                std::string("Nonconformant that in ") + __PRETTY_FUNCTION__);

        suzerain::blas::copy(
                data.num_elements(), that.data.data(), 1, data.data(), 1);
    }

    return *this;
}

template< typename FPT, bool Interleaved >
RealState<FPT,Interleaved>& RealState<FPT,Interleaved>::operator=(
        const RealState<FPT,!Interleaved>& that)
throw(std::logic_error)
{
    if (!isConformant(that)) throw std::logic_error(
            std::string("Nonconformant that in ") + __PRETTY_FUNCTION__);

    this->data = that.data;

    return *this;
}

template< typename FPT, bool Interleaved >
void RealState<FPT,Interleaved>::exchange(
        IState<FPT,Interleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (this != &other) {
        if (!isConformant(other)) throw std::logic_error(
                std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

        RealState& o = dynamic_cast<RealState&>(other);

        suzerain::blas::swap(data.num_elements(),
                             o.data.data(), 1,
                             data.data(), 1);
    }
}

template< typename FPT, bool Interleaved >
void RealState<FPT,Interleaved>::exchange(
        IState<FPT,!Interleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (!isConformant(other)) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

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
    const bool ascending[4] = { true, true, true, true };
    if (Interleaved) {
        const size_type ordering[4] = { 0, 1, 2, 3 }; // Complex variables fastest
        boost::general_storage_order<4> result(ordering, ascending);
        return result;
    } else {
        const size_type ordering[4] = { 0, 2, 3, 1 }; // Complex variables slowest
        boost::general_storage_order<4> result(ordering, ascending);
        return result;
    }
}

template< typename FPT, bool Interleaved >
boost::general_storage_order<3>
ComplexState<FPT,Interleaved>::storage_order()
{
    const bool ascending[3] = { true, true, true };
    if (Interleaved) {
        const size_type ordering[3] = { 0, 1, 2 }; // Complex variables fastest
        boost::general_storage_order<3> result(ordering, ascending);
        return result;
    } else {
        const size_type ordering[3] = { 1, 2, 0 }; // Complex variables slowest
        boost::general_storage_order<3> result(ordering, ascending);
        return result;
    }
}

template< typename FPT, bool Interleaved >
ComplexState<FPT,Interleaved>::ComplexState(
        size_type variable_count,
        size_type vector_length,
        size_type vector_count)
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
    suzerain::blas::copy(components.num_elements(),
                         other.components.data(), 1, components.data(), 1);
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
    // TODO Ensure deep assignment has adequate performance
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
throw(std::bad_cast, std::logic_error)
{
    if (!isConformant(other)) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

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
throw(std::bad_cast, std::logic_error)
{
    if (!isConformant(other)) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

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
        if (!isConformant(that)) throw std::logic_error(
                std::string("Nonconformant that in ") + __PRETTY_FUNCTION__);

        suzerain::blas::copy(components.num_elements(),
                             that.components.data(), 1, components.data(), 1);
    }

    return *this;
}

template< typename FPT, bool Interleaved >
ComplexState<FPT,Interleaved>& ComplexState<FPT,Interleaved>::operator=(
        const ComplexState<FPT,!Interleaved>& that)
throw(std::logic_error)
{
    if (!isConformant(that)) throw std::logic_error(
            std::string("Nonconformant that in ") + __PRETTY_FUNCTION__);

    // TODO Ensure deep assignment has adequate performance
    this->components = that.components;

    return *this;
}

template< typename FPT, bool Interleaved >
void ComplexState<FPT,Interleaved>::exchange(
        IState<FPT,Interleaved>& other)
throw(std::bad_cast, std::logic_error)
{
    if (this != &other) {
        if (!isConformant(other)) throw std::logic_error(
                std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

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
throw(std::bad_cast, std::logic_error)
{
    if (!isConformant(other)) throw std::logic_error(
            std::string("Nonconformant other in ") + __PRETTY_FUNCTION__);

    ComplexState<FPT,!Interleaved>& o
        = dynamic_cast<ComplexState<FPT,!Interleaved>&>(other);

    for (int i = 0; i < this->variable_count; ++i) {
        suzerain::blas::swap(this->vector_length*this->vector_count,
                             &o.data[i][0][0], o.data.strides()[1],
                             &data[i][0][0], data.strides()[1]);
    }
}

#endif // FIXME_BLOCK_DISABLED

} // namespace suzerain

#endif // __SUZERAIN_STATE_IMPL_HPP
