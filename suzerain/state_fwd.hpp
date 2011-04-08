//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2011 The PECOS Development Team
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
// state_fwd.hpp: Interfaces for manipulating mutable state vectors
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef __SUZERAIN_STATE_FWD_HPP
#define __SUZERAIN_STATE_FWD_HPP

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/contiguous_memory.hpp>
#include <suzerain/mpl.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/storage.hpp>

// TODO Document better
// TODO Allow InterleavedState/NoninterleavedState to interoperate

/** @file
 * Provide interfaces for describing and manipulating mutable state vectors.
 * Many of the classes here are designed using the Curiously Recurring
 * Template Pattern (CRTP) to provide static polymorphism following
 * the style used in Eigen.
 *
 * @see Either <a href="https://secure.wikimedia.org/wikipedia/en/wiki/Curiously_recurring_template_pattern"> Wikipedia</a> or <a href="https://secure.wikimedia.org/wikibooks/en/wiki/More_C%2B%2B_Idioms/Curiously_Recurring_Template_Pattern"> More C++ Idioms</a>
 * for more information on the Curiously Recurring Template Pattern.
 */

namespace suzerain
{

// Mandatory forward declarations
template<typename Derived> struct StateBase;
template<std::size_t NumDims, typename Element, typename Allocator>
    class NoninterleavedState;
template<std::size_t NumDims, typename Element, typename Allocator>
    class InterleavedState;

/** Implementation details for StateBase */
namespace detail
{

// Forward declaration for StateTraits
template<typename T> struct StateTraits;

/** Traits providing basic type details for NoninterleavedState */
template<std::size_t NumDims, typename Element, typename Allocator>
struct StateTraits<NoninterleavedState<NumDims,Element,Allocator> >
{
    typedef Element element;
    typedef boost::multi_array_types::index index;
    typedef boost::multi_array_types::size_type size_type;
    static const size_type dimensionality = NumDims;

private:
    StateTraits();
};

/** Traits providing basic type details for InterleavedState */
template<std::size_t NumDims, typename Element, typename Allocator>
struct StateTraits<InterleavedState<NumDims,Element,Allocator> >
{
    typedef Element element;
    typedef boost::multi_array_types::index index;
    typedef boost::multi_array_types::size_type size_type;
    static const size_type dimensionality = NumDims;

private:
    StateTraits();
};

} // end namespace detail

/**
 * Base class for all of Suzerain's state implementations.  Static polymorphism
 * via the Curiously Recurring Template Pattern is used to look up the
 * appropriate type-specific functionality.  Subclasses are usable with
 * functionality found in suzerain::timestepper::lowstorage.
 */
template<class Derived>
class StateBase
{
public:
    typedef typename detail::StateTraits<Derived>::element   element;
    typedef typename detail::StateTraits<Derived>::index     index;
    typedef typename detail::StateTraits<Derived>::size_type size_type;
    static const size_type dimensionality
        = detail::StateTraits<Derived>::dimensionality;

/** @name Static polymorphism helpers to simplify casting */
/**@{*/

    /** @returns a reference to the derived object */
    Derived& derived() {
        return *static_cast<Derived*>(this);
    }

    /** @returns a const reference to the derived object */
    const Derived& derived() const {
        return *static_cast<const Derived*>(this);
    }

/**@}*/

/** @name Queries supported by subclasses */
/**@{*/

    /** @returns each dimension's extents as a length #dimensionality list */
    const size_type* shape() const {
        return derived().shape();
    }

/**@}*/

/** @name Operations supported by subclasses */
/**@{*/

    /**
     * Scale all state information by the given scale factor, i.e.
     * \f$\mbox{this} \leftarrow{} \mbox{factor}\times\mbox{this}\f$.
     *
     * @param factor Scale factor to use.
     */
    void scale(const element& factor) {
        return derived().scale(factor);
    }

    /**
     * Is \c this instance isomorphic with <tt>other</tt>'s?
     * More concretely, do they have the same dimensionality and shape?
     *
     * @param other another instance against which to compare.
     *
     * @return True if other's shape matches this instance's.
     *         False otherwise.
     */
    template<class OtherDerived>
    bool isIsomorphic(const StateBase<OtherDerived>& other) const {
        if (dimensionality == StateBase<OtherDerived>::dimensionality) {
            return std::equal(shape(), shape() + dimensionality, other.shape());
        } else {
            return false;
        }
    }

    /**
     * Accumulate scaled state information by computing \f$\mbox{this}
     * \leftarrow{} \mbox{this} + \mbox{factor}\times\mbox{other}\f$.
     *
     * @param factor Scale factor to apply to <tt>other</tt>'s
     *               state information.
     * @param other Another state instance to scale and add to \c this.
     * @throw std::logic_error if \c other is not isomorphic.
     */
    template<class OtherDerived>
    void addScaled(const element& factor,
                   const StateBase<OtherDerived>& other) {
        return derived().addScaled(factor, other.derived());
    }

    /**
     * Assign to this instance the state information from another instance
     * holding compatibly stored elements.
     *
     * @param other instance with state to be deep copied.
     * @throw std::logic_error if \c other is not isomorphic.
     */
    template<class OtherDerived>
    void assign(const StateBase<OtherDerived> &other) {
        return derived().assign(other.derived());
    }

    /**
     * Exchange <tt>this</tt>'s storage with <tt>other</tt>'s storage by moving
     * the relevant data without incurring a lot of memory overhead.  In
     * contrast with potentially optimized swap semantics, this method should
     * swap data instead of playing tricks with underlying pointers.
     *
     * @param other instance with which to exchange data.
     * @throw std::logic_error if \c other is not isomorphic.
     */
    template<class OtherDerived>
    void exchange(StateBase<OtherDerived>& other) {
        return derived().exchange(other.derived());
    }

/**@}*/

};


/**
 * A state class storing elements with the first index \em slowest followed by
 * the other indices in Fortran order.  The indexing is indented to store some
 * number of non-interleaved scalars in the first index with dimensionality
 * <tt>NumDims</tt> spread across the remaining indices.
 *
 * @tparam NumDims    Number of dimensions
 * @tparam Element    Type of elements to store
 * @tparam Allocator  Allocator for element storage
 */
template<
    std::size_t NumDims,
    typename Element,
    typename Allocator = typename suzerain::blas::allocator<Element>::type
>
class NoninterleavedState
    : public StateBase<NoninterleavedState<NumDims,Element,Allocator> >,
      public ContiguousMemory<Element,Allocator>,
      public suzerain::multi_array::ref<Element, NumDims>
{
public:

/** @name Declarations bringing in information from public base classes */
/**@{*/
    typedef typename suzerain::multi_array::ref<Element,NumDims> multi_array_type;
    typedef typename multi_array_type::const_reference      const_reference;
    typedef typename multi_array_type::difference_type      difference_type;
    typedef typename multi_array_type::element              element;
    typedef typename multi_array_type::index                index;
    typedef typename multi_array_type::reference            reference;
    typedef typename multi_array_type::size_type            size_type;
    typedef typename multi_array_type::value_type           value_type;
    static const size_type dimensionality = multi_array_type::dimensionality;
    using multi_array_type::shape;
/**@}*/

    /** The storage ordering in use */
    typedef typename storage::noninterleaved<NumDims> storage_type;

    template<typename ExtentList>
    explicit NoninterleavedState(const ExtentList& sizes);

    template<typename ExtentList, typename MinStrideList>
    explicit NoninterleavedState(const ExtentList& sizes,
                                 const MinStrideList& minstrides);

    NoninterleavedState(const NoninterleavedState& other);

    void scale(const Element& factor);

    void addScaled(const Element& factor,
                   const NoninterleavedState& other);
    void assign(const NoninterleavedState &other);

    void exchange(NoninterleavedState &other);

private:
    // Disable assignment operators
    const multi_array_type& operator=( const multi_array_type& );
    const NoninterleavedState& operator=( const NoninterleavedState& );
};

/**
 * A state class storing elements in Fortran order.  The class is named
 * InterleavedState to contrast it with NoninterleavedState's storage ordering.
 *
 * @tparam NumDims    Number of dimensions
 * @tparam Element    Type of elements to store
 * @tparam Allocator  Allocator for element storage
 */
template<
    std::size_t NumDims,
    typename Element,
    typename Allocator = typename suzerain::blas::allocator<Element>::type
>
class InterleavedState
    : public StateBase<InterleavedState<NumDims,Element,Allocator> >,
      public ContiguousMemory<Element,Allocator>,
      public suzerain::multi_array::ref<Element, NumDims>
{
public:

/** @name Declarations bringing in information from public base classes */
/**@{*/
    typedef typename suzerain::multi_array::ref<Element,NumDims> multi_array_type;
    typedef typename multi_array_type::const_reference      const_reference;
    typedef typename multi_array_type::difference_type      difference_type;
    typedef typename multi_array_type::element              element;
    typedef typename multi_array_type::index                index;
    typedef typename multi_array_type::reference            reference;
    typedef typename multi_array_type::size_type            size_type;
    typedef typename multi_array_type::value_type           value_type;
    static const size_type dimensionality = multi_array_type::dimensionality;
    using multi_array_type::shape;
/**@}*/

    /** The storage ordering in use */
    typedef typename storage::interleaved<NumDims> storage_type;

    template<typename ExtentList>
    explicit InterleavedState(
            const ExtentList& sizes,
            size_type min_total_contiguous_count = 0);

    InterleavedState(const InterleavedState& other);

    void scale(const Element& factor);

    void addScaled(const Element& factor,
                   const InterleavedState& other);

    void assign(const InterleavedState &other);

    void exchange(InterleavedState &other);

private:
    // Disable assignment operators
    const multi_array_type& operator=( const multi_array_type& );
    const InterleavedState& operator=( const InterleavedState& );
};

} // namespace suzerain

#endif // __SUZERAIN_STATE_FWD_HPP
