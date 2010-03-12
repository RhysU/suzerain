/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * state.hpp: Class to manage mutable state vectors
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
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
 * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
 *      Boost.MultiArray</a> for more information on the MultiArray concept.
 */
template< typename FPT >
class IState
{
public:
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

    /** Virtual destructor appropriate for an abstract base class */
    virtual ~IState() {}

    /**
     * Is \c this instance's shape "conformant" with <tt>other</tt>'s?
     *
     * @param other another instance to compare against.
     *
     * @return True if all of #variable_count, #vector_length, and
     *         #vector_count are identical.  False otherwise.
     */
    virtual bool isConformant(const IState &other) const
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
                           const IState<FPT> &other)
                           throw(std::bad_cast, std::logic_error) = 0;

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     *
     * @return *this
     */
    virtual IState& operator=(const IState<FPT>& that)
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
    virtual void exchange(IState<FPT>& other)
                          throw(std::bad_cast, std::logic_error) = 0;

    // TODO Add IState<FPT>::swap(IState<FPT>&) to public interface
};

/**
 * An implementation of IState<FPT> for real-valued state information.
 */
template< typename FPT >
class RealState : public IState<FPT>
{
public:
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
    explicit RealState(typename IState<FPT>::size_type variable_count,
                       typename IState<FPT>::size_type vector_length,
                       typename IState<FPT>::size_type vector_count)
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
                           const IState<FPT> &other)
                           throw(std::bad_cast, std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::logic_error if \c other is not conformant.
     */
    RealState& operator=(const RealState& that) throw(std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant.
     */
    virtual RealState& operator=(const IState<FPT>& that)
                                 throw(std::bad_cast, std::logic_error)
    {
        return operator=(dynamic_cast<const RealState<FPT>&>(that));
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
    virtual void exchange(IState<FPT>& other)
                          throw(std::bad_cast, std::logic_error);

    // TODO Implement swap using https://svn.boost.org/trac/boost/ticket/1045

protected:
    /**
     * Raw, aligned storage of all state information stored in column-major
     * storage over indices IState<FPT>#variable_count,
     * IState<FPT>#vector_length, and IState<FPT>#vector_count.
     **/
    boost::shared_array<double> raw;

public:
    /** Type of the Boost.MultiArray exposed through RealState#data. */
    typedef typename boost::multi_array_ref<FPT, 3> state_type;

    /**
     * A three-dimensional Boost.MultiArray of real-valued state kept
     * zero-indexed in column-major order.  The indices are over
     * #variable_count, #vector_length, and #vector_count respectively.
     * The information is contiguous.
     */
    state_type data;
};

template< typename FPT >
RealState<FPT>::RealState(typename IState<FPT>::size_type variable_count,
                          typename IState<FPT>::size_type vector_length,
                          typename IState<FPT>::size_type vector_count)
throw(std::bad_alloc)
    : IState<FPT>(variable_count, vector_length, vector_count),
      raw(reinterpret_cast<FPT *>(suzerain_blas_malloc(
                  sizeof(state_type::element)
                * variable_count*vector_length*vector_count)),
          std::ptr_fun(free)),
      data(raw.get(),
           boost::extents[variable_count][vector_length][vector_count],
           boost::fortran_storage_order())
{
    if (!raw) throw std::bad_alloc();
}

template< typename FPT >
RealState<FPT>::RealState(const RealState<FPT> &other)
throw(std::bad_alloc)
    : IState<FPT>(other),
      raw(reinterpret_cast<FPT *>(suzerain_blas_malloc(
                sizeof(state_type::element)
              * other.variable_count
              * other.vector_length
              * other.vector_count)),
          std::ptr_fun(free)),
      data(raw.get(),
           boost::extents[other.variable_count]
                         [other.vector_length]
                         [other.vector_count],
           boost::fortran_storage_order())
{
    if (!raw) throw std::bad_alloc();
    memcpy(data.data(), other.data.data(),
           sizeof(state_type::element)*data.num_elements());
}

template< typename FPT >
void RealState<FPT>::scale(const FPT factor)
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

template< typename FPT >
void RealState<FPT>::addScaled(const FPT factor,
                               const IState<FPT> &other)
throw(std::bad_cast,
      std::logic_error)
{
    if (!isConformant(other))
        throw std::logic_error("Nonconformant other in addScaled");

    const RealState<FPT> &o = dynamic_cast<const RealState<FPT>&>(other);

    suzerain::blas::axpy(data.num_elements(),
                         factor, o.data.data(), 1,
                         data.data(), 1);
}

template< typename FPT >
RealState<FPT>& RealState<FPT>::operator=(const RealState<FPT> &that)
throw(std::logic_error)
{
    if (this != &that) {
        if (!isConformant(that))
            throw std::logic_error("Nonconformant that in operator=");

        memcpy(data.data(), that.data.data(),
               sizeof(state_type::element)*data.num_elements());
    }

    return *this;
}

template< typename FPT >
void RealState<FPT>::exchange(IState<FPT> &other)
throw(std::bad_cast,
      std::logic_error)
{
    if (this != &other) {
        if (!isConformant(other))
            throw std::logic_error("Nonconformant other in exchange");

        RealState<FPT> &o = dynamic_cast<RealState<FPT>&>(other);

        suzerain::blas::swap(data.num_elements(),
                             o.data.data(), 1,
                             data.data(), 1);
    }
}

/**
 * An implementation of IState<FPT> for complex-valued state information.
 * \c FPT is the real scalar type used underneath the complex type.
 */
template< typename FPT >
class ComplexState : public IState<FPT>
{
public:
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
    explicit ComplexState(typename IState<FPT>::size_type variable_count,
                          typename IState<FPT>::size_type vector_length,
                          typename IState<FPT>::size_type vector_count)
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
                           const IState<FPT> &other)
                           throw(std::bad_cast,
                                 std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::logic_error if \c other is not conformant.
     */
    ComplexState& operator=(const ComplexState& that) throw(std::logic_error);

    /**
     * Assign to this instance the state information from another instance.
     *
     * @param that instance to copy.
     * @return *this
     * @throw std::bad_cast if \c other does not have a compatible type.
     * @throw std::logic_error if \c other is not conformant.
     */
    virtual ComplexState& operator=(const IState<FPT>& that)
                                    throw(std::bad_cast, std::logic_error)
    {
        return operator=(dynamic_cast<const ComplexState<FPT>&>(that));
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
    virtual void exchange(IState<FPT>& other)
                          throw(std::bad_cast, std::logic_error);

    // TODO Implement swap using https://svn.boost.org/trac/boost/ticket/1045

protected:
    /**
     * Raw, aligned storage of all state information stored in column-major
     * storage over \c FPT real-valued components (i.e. real and imaginary),
     * #variable_count, #vector_length, and #vector_count.
     **/
    boost::shared_array<FPT> raw;

public:
    /** Type of the Boost.MultiArray exposed through #components. */
    typedef typename boost::multi_array_ref<FPT, 4> components_type;

    /** Type of the Boost.MultiArray exposed through #real and #imag. */
    typedef typename boost::array_view_gen<components_type,3>::type component_type;

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

template< typename FPT >
ComplexState<FPT>::ComplexState(typename IState<FPT>::size_type variable_count,
                                typename IState<FPT>::size_type vector_length,
                                typename IState<FPT>::size_type vector_count)
throw(std::bad_alloc)
    : IState<FPT>(variable_count, vector_length, vector_count),
      raw(reinterpret_cast<FPT *>(suzerain_blas_malloc(
                  sizeof(state_type::element)
                * variable_count*vector_length*vector_count)),
          std::ptr_fun(free)),
      components(raw.get(),
           boost::extents[2][variable_count][vector_length][vector_count],
           boost::fortran_storage_order()),
      real(components[boost::indices[0]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      imag(components[boost::indices[1]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      data(reinterpret_cast<std::complex<FPT> *>(raw.get()),
           boost::extents[variable_count][vector_length][vector_count],
           boost::fortran_storage_order())
{
    if (!raw) throw std::bad_alloc();
}

template< typename FPT >
ComplexState<FPT>::ComplexState(const ComplexState<FPT> &other)
throw(std::bad_alloc)
    : IState<FPT>(other),
      raw(reinterpret_cast<FPT *>(suzerain_blas_malloc(
                   sizeof(state_type::element)
                 * other.variable_count
                 * other.vector_length
                 * other.vector_count)),
          std::ptr_fun(free)),
      components(raw.get(),
           boost::extents[2]
                         [other.variable_count]
                         [other.vector_length]
                         [other.vector_count],
           boost::fortran_storage_order()),
      real(components[boost::indices[0]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      imag(components[boost::indices[1]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]
                                    [boost::multi_array_types::index_range()]]),
      data(reinterpret_cast<std::complex<FPT> *>(raw.get()),
           boost::extents[other.variable_count]
                         [other.vector_length]
                         [other.vector_count],
           boost::fortran_storage_order())
{
    if (!raw) throw std::bad_alloc();
    memcpy(data.data(), other.data.data(),
           sizeof(state_type::element)*data.num_elements());
}

template< typename FPT >
void ComplexState<FPT>::scale(const FPT factor)
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

template< typename FPT >
void ComplexState<FPT>::addScaled(const FPT factor,
                                  const IState<FPT> &other)
throw(std::bad_cast,
      std::logic_error)
{
    if (!isConformant(other))
        throw std::logic_error("Nonconformant other in addScaled");

    const ComplexState<FPT> &o = dynamic_cast<const ComplexState<FPT>&>(other);

    suzerain::blas::axpy(components.num_elements(),
                         factor, o.components.data(), 1,
                         components.data(), 1);
}

template< typename FPT >
ComplexState<FPT>& ComplexState<FPT>::operator=(const ComplexState<FPT> &that)
throw(std::logic_error)
{
    if (this != &that) {
        if (!isConformant(that))
            throw std::logic_error("Nonconformant that in operator=");

        memcpy(data.data(), that.data.data(),
               sizeof(state_type::element)*data.num_elements());
    }

    return *this;
}

template< typename FPT >
void ComplexState<FPT>::exchange(IState<FPT> &other)
throw(std::bad_cast,
      std::logic_error)
{
    if (this != &other) {
        if (!isConformant(other))
            throw std::logic_error("Nonconformant other in exchange");

        ComplexState<FPT> &o = dynamic_cast<ComplexState<FPT>&>(other);

        suzerain::blas::swap(components.num_elements(),
                             o.components.data(), 1,
                             components.data(), 1);
    }
}

} // namespace suzerain

#endif // __SUZERAIN_STATE_HPP
