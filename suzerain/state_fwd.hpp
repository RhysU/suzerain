//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#ifndef SUZERAIN_STATE_FWD_HPP
#define SUZERAIN_STATE_FWD_HPP

/** @file
 * Provide interfaces for describing and manipulating mutable state vectors.
 * Many of the classes here are designed using the Curiously Recurring
 * Template Pattern (CRTP) to provide static polymorphism following
 * the style used in Eigen.
 *
 * @see Both <a
 * href="https://secure.wikimedia.org/wikipedia/en/wiki/Curiously_recurring_template_pattern">
 * Wikipedia</a> and <a
 * href="https://secure.wikimedia.org/wikibooks/en/wiki/More_C%2B%2B_Idioms/Curiously_Recurring_Template_Pattern">
 * More C++ Idioms</a> for more information on the Curiously Recurring Template
 * Pattern.
 */

// TODO Document better

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/mpl.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/shared_range.hpp>
#include <suzerain/storage.hpp>

namespace suzerain {

// Mandatory forward declarations
template<typename Derived> struct state_base;
template<std::size_t Dim, typename Element> class contiguous_state;
template<std::size_t Dim, typename Element> class interleaved_state;

/** Implementation details for state_base */
namespace detail {

// Forward declaration for state_traits
template<typename T> struct state_traits;

/** Traits providing basic type details for contiguous_state */
template<std::size_t Dim, typename Element>
struct state_traits<contiguous_state<Dim,Element> >
{
    typedef Element element;
    typedef boost::multi_array_types::index index;
    typedef boost::multi_array_types::size_type size_type;
    static const size_type dimensionality = Dim;

private:
    state_traits();
};

/** Traits providing basic type details for interleaved_state */
template<std::size_t Dim, typename Element>
struct state_traits<interleaved_state<Dim,Element> >
{
    typedef Element element;
    typedef boost::multi_array_types::index index;
    typedef boost::multi_array_types::size_type size_type;
    static const size_type dimensionality = Dim;

private:
    state_traits();
};

} // end namespace detail

/**
 * Base class for all of Suzerain's state implementations.  Static polymorphism
 * via the Curiously Recurring Template Pattern is used to look up the
 * appropriate type-specific functionality.  Subclasses are usable with
 * functionality found in \ref lowstorage.
 */
template<class Derived>
class state_base
{
public:
    typedef typename detail::state_traits<Derived>::element   element;
    typedef typename detail::state_traits<Derived>::index     index;
    typedef typename detail::state_traits<Derived>::size_type size_type;
    static const size_type dimensionality
        = detail::state_traits<Derived>::dimensionality;

/** @name Static polymorphism helpers to simplify casting */
/**@{*/

    /** @returns a reference to the derived object */
    Derived& derived()
    {
        return *static_cast<Derived*>(this);
    }

    /** @returns a const reference to the derived object */
    const Derived& derived() const
    {
        return *static_cast<const Derived*>(this);
    }

/**@}*/

/** @name Queries supported by subclasses */
/**@{*/

    /**
     * @returns each dimension's extents as a list of length
     * state_base::dimensionality.
     */
    const size_type* shape() const
    {
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
    void scale(const element& factor)
    {
        return derived().scale(factor);
    }

    /**
     * Is \c this instance isomorphic with <tt>other</tt>'s?
     * More concretely, do they have the same dimensionality and shape?
     *
     * Implementation works on all <tt>state_base<Derived></tt> subclasses
     * and everything adhering to the <tt>Boost.MultiArray</tt> concept.
     *
     * @param other another object instance against which to compare.
     *
     * @return True if other's shape matches this instance's.
     *         False otherwise.
     */
    template<class OtherType>
    bool is_isomorphic(const OtherType& other) const
    {
        if (dimensionality == OtherType::dimensionality) {
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
    void add_scaled(const element& factor,
                    const state_base<OtherDerived>& other)
    {
        return derived().add_scaled(factor, other.derived());
    }

    /**
     * Assign to this instance the state information from another instance
     * holding compatibly stored elements.
     *
     * @param other instance with state to be deep copied.
     * @throw std::logic_error if \c other is not isomorphic.
     */
    template<class OtherDerived>
    void assign_from(const state_base<OtherDerived>& other)
    {
        return derived().assign_from(other.derived());
    }

    /**
     * Exchange <tt>this</tt>'s storage with <tt>other</tt>'s storage by moving
     * the relevant data while incurring minimal memory overhead.  In contrast
     * with potentially optimized swap semantics, this method might swap data
     * instead of playing tricks with underlying pointers.
     *
     * @param other instance with which to exchange data.
     * @throw std::logic_error if \c other is not isomorphic.
     */
    template<class OtherDerived>
    void exchange(state_base<OtherDerived>& other)
    {
        return derived().exchange(other.derived());
    }

/**@}*/

};


/**
 * A state class storing elements with the first index \em slowest followed by
 * the other indices in Fortran order.  The indexing is indented to store some
 * number of non-interleaved scalars in the first index with dimensionality \c
 * Dim spread across the remaining indices.
 *
 * @tparam Dim     Number of dimensions
 * @tparam Element Type of elements to store
 *
 * @see #storage_order_type for more details on the storage used.
 */
template< std::size_t Dim, typename Element >
class contiguous_state
    : public  state_base<contiguous_state<Dim,Element> >,
      private shared_range<Element>,
      public  suzerain::multi_array::ref<Element, Dim>
{
public:

/** @name Declarations bringing in information from public base classes */
/**@{*/
    typedef typename suzerain::multi_array::ref<Element,Dim> multi_array_type;
    typedef typename multi_array_type::const_iterator        const_iterator;
    typedef typename multi_array_type::const_reference       const_reference;
    typedef typename multi_array_type::difference_type       difference_type;
    typedef typename multi_array_type::element               element;
    typedef typename multi_array_type::index                 index;
    typedef typename multi_array_type::iterator              iterator;
    typedef typename multi_array_type::reference             reference;
    typedef typename multi_array_type::size_type             size_type;
    typedef typename multi_array_type::value_type            value_type;
    typedef typename suzerain::shared_range<Element>         shared_range_type;
    static const size_type dimensionality = multi_array_type::dimensionality;
    using multi_array_type::shape;
/**@}*/

    /** The storage ordering in use */
    typedef typename storage::contiguous<Dim> storage_order_type;

    template<typename ExtentList>
    explicit contiguous_state(const ExtentList& sizes);

    template<typename ExtentList>
    contiguous_state(const shared_range_type& storage,
                     const ExtentList& sizes);

    template<typename ExtentList, typename MinStrideList>
    contiguous_state(const ExtentList& sizes,
                     const MinStrideList& minstrides);

    template<typename ExtentList, typename MinStrideList>
    contiguous_state(const shared_range_type& storage,
                     const ExtentList& sizes,
                     const MinStrideList& minstrides);

    contiguous_state(const contiguous_state& other);

    void scale(const Element& factor);

    void add_scaled(const Element& factor,
                    const contiguous_state& other);

    void add_scaled(const Element& factor,
                    const multi_array_type& other);

    void assign_from(const contiguous_state& other);

    void assign_from(const multi_array_type& other);

    void exchange(contiguous_state& other);

    void exchange(multi_array_type& other);

    const shared_range_type& range() const
    {
        return reinterpret_cast<const shared_range_type&>(*this);
    }

    // Clarify ambiguities due to private inheritance from shared_range_type
    using multi_array_type::operator[];
    using multi_array_type::begin;
    using multi_array_type::end;
    using multi_array_type::size;

private:
    // Disable assignment operators
    const multi_array_type& operator=( const multi_array_type& );
    const contiguous_state& operator=( const contiguous_state& );
};

/**
 * A state class storing elements in Fortran order \em except that the
 * first two indices are exchanged.
 *
 * @tparam Dim     Number of dimensions
 * @tparam Element Type of elements to store
 *
 * @see #storage_order_type for more details on the storage used.
 */
template< std::size_t Dim, typename Element >
class interleaved_state
    : public  state_base<interleaved_state<Dim,Element> >,
      private shared_range<Element>,
      public  suzerain::multi_array::ref<Element, Dim>
{
public:

/** @name Declarations bringing in information from public base classes */
/**@{*/
    typedef typename suzerain::multi_array::ref<Element,Dim> multi_array_type;
    typedef typename multi_array_type::const_iterator        const_iterator;
    typedef typename multi_array_type::const_reference       const_reference;
    typedef typename multi_array_type::difference_type       difference_type;
    typedef typename multi_array_type::element               element;
    typedef typename multi_array_type::index                 index;
    typedef typename multi_array_type::iterator              iterator;
    typedef typename multi_array_type::reference             reference;
    typedef typename multi_array_type::size_type             size_type;
    typedef typename multi_array_type::value_type            value_type;
    typedef typename suzerain::shared_range<Element>         shared_range_type;
    static const size_type dimensionality = multi_array_type::dimensionality;
    using multi_array_type::shape;
/**@}*/

    /** The storage ordering in use */
    typedef typename storage::interleaved<Dim> storage_order_type;

    template<typename ExtentList>
    explicit interleaved_state(const ExtentList& sizes,
                               size_type min_total_contiguous_count = 0);

    template<typename ExtentList>
    interleaved_state(shared_range_type storage,
                      const ExtentList& sizes,
                      size_type min_total_contiguous_count = 0);

    interleaved_state(const interleaved_state& other);

    void scale(const Element& factor);

    void add_scaled(const Element& factor,
                    const interleaved_state& other);

    void add_scaled(const Element& factor,
                    const multi_array_type& other);

    void assign_from(const interleaved_state& other);

    void assign_from(const multi_array_type& other);

    void exchange(interleaved_state& other);

    void exchange(multi_array_type& other);

    const shared_range_type& range() const
    {
        return reinterpret_cast<const shared_range_type&>(*this);
    }

    // Clarify ambiguities due to private inheritance from shared_range_type
    using multi_array_type::operator[];
    using multi_array_type::begin;
    using multi_array_type::end;
    using multi_array_type::size;

private:
    // Disable assignment operators
    const multi_array_type& operator=( const multi_array_type& );
    const interleaved_state& operator=( const interleaved_state& );
};

} // namespace suzerain

#endif // SUZERAIN_STATE_FWD_HPP
