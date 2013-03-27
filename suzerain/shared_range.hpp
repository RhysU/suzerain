// Copyright Olaf van der Spek and Rhys Ulerich (2011).
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See http://www.boost.org/LICENSE_1_0.txt)

#ifndef SUZERAIN_SHARED_RANGE_HPP
#define SUZERAIN_SHARED_RANGE_HPP

/*! \file
 *  \brief Defines \c shared_range and related functions.
 */

#include <algorithm>
#include <memory>
#include <boost/range/distance.hpp>
#include <boost/range/difference_type.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/value_type.hpp>
#include <boost/shared_array.hpp>
#include <boost/swap.hpp>
#include <boost/version.hpp>

namespace suzerain {

// IMPLEMENTATION COMMENTS:
//
// Olaf van der Spek suggested this concept to boost-users
// (http://boost.2283326.n4.nabble.com/RFC-A-better-shared-array-td3895154.html)
// and provided an initial implementation at http://pastebin.com/wkdLVqM1.
//
// Per discussion with Peter Dimov in that thread,
// (shared_array/iterator_range) state was used for more generality instead of
// the more space efficient (shared_array/size_t) state choice.  The assumption
// is that someone doing lots of lightweight "range work" is explicitly
// handling resource management while someone doing someone doing lots of
// lightweight "array work" is similarly explicitly handling known lengths.
//
// More of the shared_array stored as shared_range.p_ (e.g. use_count(), get())
// was deliberately NOT exposed as I could not reconcile those semantics with a
// desired to have shared_range truly be an iterator_range (and therefore
// support construction and assignment from a generic ForwardRange).
//
// Full shared_ptr-like Allocator support (in particular, simultaneous header
// and storage allocation) seems to require
// https://svn.boost.org/trac/boost/ticket/5924 which (small world) Olaf also
// requested.
//
// See http://www.boost.org/doc/libs/release/libs/range/doc/html/range/upgrade.html
// for version-to-version changes in Boost.Range that impact the
// implementation.  In particular, pre-1.43 there are many BOOST_ASSERT(
// !is_singular() ) calls that interfere with the empty range cases below.
// Workarounds use preprocessor keyed on BOOST_VERSION to mitigate problems.
// Several places in the implementation rely on <algorithm> rather than Range
// algorithms to improve version-independence.
//
// No enclosing namespace is assumed in the implementation.  Hence ::std:: and
// ::boost:: appears where one might expect std:: and boost::.  The
// implementation therefore can be #included into any namespace.

/*!
 * \brief A \c shared_range extends the functionality of \c iterator_range by
 * adding smart pointer-like resource ownership semantics.
 *
 * Semantically, \c shared_range resembles a <tt>shared_array<T></tt> which
 * "knows" its own length.  Instances assume ownership of pointers provided to
 * constructors and guarantee resources will be deleted when the last \c
 * shared_range pointing to them is destroyed or reset.
 */
template< class T >
class shared_range
    : public ::boost::iterator_range<T*>
{
public: // typedefs

    //! Base class type
    typedef typename ::boost::iterator_range<T*>::iterator_range iterator_range;

    //! Encapsulated value type (to mimic shared_array, shared_ptr).
    typedef T element_type;

public: // construction

    //! Default constructor produces an empty range instance.
    shared_range()
        : iterator_range(),
          p_()
    {}

    //! Constructor from a <tt>new[]</tt>-allocated pointer and size.
    //  See shared_array(T*) for resource management details.
    shared_range(T *b, ::std::size_t sz)
        : iterator_range(b, b + sz),
          p_(b)
    {}

    //! Constructor from a pointer and size using a custom deleter.
    //  See shared_array(T*,D) for resource management details.
    template< class D >
    shared_range(T *b, ::std::size_t sz, D d)
        : iterator_range(b, b + sz),
          p_(b, d)
    {}

    //! Constructor from the half closed, contiguous range <tt>[b,e)</tt>.
    //  See shared_array(T*) for resource management details of \c b.
    shared_range(T *b, T *e)
        : iterator_range(b, e),
          p_(b)
    {}

    //! Constructor from the half closed, contiguous range <tt>[b,e)</tt>.
    //  See shared_array(T*, D) for resource management details of \c b.
    template< class D >
    shared_range(T *b, T *e, D d)
        : iterator_range(b, e),
          p_(b, d)
    {}

    //! Constructor from an existing \c shared_array with specified size.
    //  See shared_array(shared_array const &) for resource management details.
    shared_range(const ::boost::shared_array<T>& p, ::std::size_t sz)
        : iterator_range(p.get(), p.get() + sz),
          p_(p)
    {}

    //! Constructor from a Range providing no ownership semantics.
    //  Provided for compatibility with \c iterator_range.
    template< class ForwardRange >
    explicit shared_range( ForwardRange& r )
        : iterator_range(r),
          p_()  // resource not managed
    {}

    //! Constructor from a Range providing no ownership semantics.
    //  Provided for compatibility with \c iterator_range.
    template< class ForwardRange >
    explicit shared_range( const ForwardRange& r )
        : iterator_range(r),
          p_()  // resource not managed
    {}

public: // copy construction

    //! Construct an instance sharing ownership with \c o.
    shared_range(shared_range& o)
        : iterator_range(o),
          p_(o.p_)
    {}

    //! Construct an instance sharing ownership with \c o.
    shared_range(const shared_range& o)
        : iterator_range(o),
          p_(o.p_)
    {}

public: // copy assignment, assignment

    //! Releases current range ownership to acquires shared ownership of \c o.
    shared_range& operator=(shared_range o) // pass-by-value for elision
    {
        this->swap(o);
        return *this;
    }

    //! Release current range ownership to assign from a Range.
    //  This assignment provides no ownership semantics.
    //  Provided for compatibility with \c iterator_range.
    template< class ForwardRange >
    shared_range& operator=( ForwardRange& r )
    {
        this->iterator_range::operator=(r);
        p_.reset();  // resource not managed
        return *this;
    }

    //! Assign from a Range and provided no ownership semantics.
    //  This assignment provides no ownership semantics.
    //  Provided for compatibility with \c iterator_range.
    template< class ForwardRange >
    shared_range& operator=( const ForwardRange& r )
    {
        this->iterator_range::operator=(r);
        p_.reset();  // resource not managed
        return *this;
    }

public: // swap, reset

    //! Exchange the contents of this instance with \c o.
    void swap(shared_range& o)
    {
#if BOOST_VERSION >= 104300
        // It would be nice if iterator_range provided a member swap...
        ::boost::swap(static_cast<iterator_range&>(*this),
                      static_cast<iterator_range&>(o));
#else
        // Perform a "logical" swap without triggering is_singular() asserts.
        switch (this->is_singular() + (o.is_singular() << 1)) {
        case 0:  // neither singular
            ::boost::swap(static_cast<iterator_range&>(*this),
                          static_cast<iterator_range&>(o));
            break;
        case 1:  // only this singular so copy o and make o singular
            this->iterator_range::operator=(o);
            o.iterator_range::advance_begin(::boost::distance(o));
            break;
        case 2:  // only o singular so copy this and make this singular
            o.iterator_range::operator=(*this);
            iterator_range::advance_begin(::boost::distance(*this));
            break;
        case 3:  // both singular so do nothing
            break;
        }
#endif
        p_.swap(o.p_);
    }

    //! Release ownership of any currently owned range.
    void reset()
    {
#if BOOST_VERSION >= 104300
        shared_range().swap(*this);
#else
        if (!iterator_range::is_singular()) {
            iterator_range::advance_begin(::boost::distance(*this));
        }
        p_.reset();
#endif
    }

    //! Construct a new instance per the constructor with the same signature
    //  and replace this instance with the new one.
    void reset(T *b, ::std::size_t sz)
    {
        shared_range(b, sz).swap(*this);
    }

    //! Construct a new instance per the constructor with the same signature
    //  and replace this instance with the new one.
    template< class D >
    void reset(T *b, ::std::size_t sz, D d)
    {
        shared_range(b, sz, d).swap(*this);
    }

    //! Construct a new instance per the constructor with the same signature
    //  and replace this instance with the new one.
    void reset(T *b, T *e)
    {
        shared_range(b, e).swap(*this);
    }

    //! Construct a new instance per the constructor with the same signature
    //  and replace this instance with the new one.
    template< class D >
    void reset(T *b, T *e, D d)
    {
        shared_range(b, e, d).swap(*this);
    }

public: // shared_array-like functionality

    //! Does this instance uniquely maintain ownership of any resources?
    bool unique() const { return p_ && p_.unique(); }

public: // iterator_range methods (updated for workarounds or covariant return)

// These methods circumnavigate earlier iterator_range is_singular() asserts
#if BOOST_VERSION < 104300
    bool empty() const
    {
        return iterator_range::is_singular() || iterator_range::empty();
    }

    operator typename iterator_range::unspecified_bool_type() const
    {
        return empty() ? 0 : &iterator_range::end;
    }
#endif

    //! Mutate instance by advancing \c begin() by \c n.
    shared_range& advance_begin(typename iterator_range::difference_type n)
    {
        iterator_range::advance_begin(n);
        return *this;
    }

    //! Mutate instance by advancing \c end() by \c n.
    shared_range& advance_end(typename iterator_range::difference_type n)
    {
        iterator_range::advance_end(n);
        return *this;
    }

private:

    //! Performs shared resource ownership semantics.
    ::boost::shared_array<T> p_;

};

//! Free swap delegating to the member swap implementation.
template< class T >
inline void swap(shared_range<T>& a, shared_range<T>& b)
{
    a.swap(b);
}

//! Are \c l and \c r semantically equivalent?
//  To check if \c l and \c r are identical, use <tt>l.equal(r)</tt>.
template< class T >
inline bool operator==(const shared_range<T>& l, const shared_range<T>& r)
{
    if (l.equal(r)) return true;  // Short circuit full range comparison

    typedef typename shared_range<T>::iterator_range iterator_range;
    return    static_cast<const iterator_range&>(l)
           == static_cast<const iterator_range&>(r);
}

//! Are \c l and \c r not semantically equivalent?
//  To check if \c l and \c r are not identical, use <tt>!l.equal(r)</tt>.
template< class T >
inline bool operator!=(const shared_range<T>& l, const shared_range<T>& r)
{
    return !(l == r);
}

//! Is \c l semantically less than \c r?
template< class T >
inline bool operator<(const shared_range<T>& l, const shared_range<T>& r)
{
    typedef typename shared_range<T>::iterator_range iterator_range;
    return   static_cast<const iterator_range&>(l)
           < static_cast<const iterator_range&>(r);
}

//! Allocate and manage ownership a <tt>new T[]</tt>-ed array of size \c sz.
template< class T >
shared_range<T> make_shared_range(::std::size_t sz)
{
    shared_range<T> t(new T[sz], sz);
    return t;
}

namespace {

// Preserves allocation size information for use in Allocator::deallocate.
// Necessary to use a custom allocator in conjunction with shared_array.
template<class Allocator>
struct allocator_deleter : public Allocator
{
    typename Allocator::size_type sz;

    allocator_deleter(typename Allocator::size_type sz) : sz(sz) {}

    void operator()(typename Allocator::pointer p)
    {
        if (p) {
            for (typename Allocator::size_type i = 0; i < sz; ++i)
                Allocator::destroy(p + i);

            Allocator::deallocate(p, sz);
        }
    }
};

}

//! Use \c a to allocate a default-constructed \c shared_range of size \c sz.
template< class A >
shared_range<typename A::value_type> allocate_shared_range(
        A a,
        ::std::size_t sz)
{
    typedef typename A::value_type value_type;
    shared_range<value_type> t(a.allocate(sz), sz, allocator_deleter<A>(sz));
    ::std::uninitialized_fill(t.begin(), t.end(), value_type());
    return t;
}

//! Use \c a to allocate a copy-constructed \c shared_range of size \c sz.
//  Argument \c val is the copy construction source.  Allocator \c a
//  is \em not used to allocate the shared_range details, only the storage
//  used to contain the range.
template< class A, class V >
shared_range<typename A::value_type> allocate_shared_range(
        A a,
        ::std::size_t sz,
        const V& val,
        typename A::const_pointer hint = 0)
{
    shared_range<typename A::value_type> t(
            a.allocate(sz, hint), sz, allocator_deleter<A>(sz));
    ::std::uninitialized_fill(t.begin(), t.end(), val);
    return t;
}

//! Clone a ForwardRange.
template< class ForwardRange >
shared_range<typename ::boost::range_value<ForwardRange>::type>
clone_shared_range(const ForwardRange &r)
{
    // TODO Employ uninitialized_copy without making allocator assumptions?
    // See http://stackoverflow.com/questions/7822925 for details/concerns.
    // Would become a single pass operation at the cost of a custom deleter.
    const ::std::size_t sz = ::boost::distance(r);
    typedef typename ::boost::range_value<ForwardRange>::type value_type;
    shared_range<value_type> t = make_shared_range<value_type>(sz);
    ::std::copy(::boost::begin(r), ::boost::end(r), t.begin());
    return t;
}

//! Clone a ForwardRange using allocator \c a to obtain storage.
//  Allocator \c a is \em not used to allocate the shared_range details, only
//  the storage used to contain the range.
template< class A, class ForwardRange >
shared_range<typename A::value_type> clone_shared_range(
        A a, const ForwardRange &r, typename A::const_pointer hint = 0)
{
    const typename A::size_type sz = ::boost::distance(r);
    shared_range<typename A::value_type> t(
            a.allocate(sz, hint), sz, allocator_deleter<A>(sz));
    ::std::uninitialized_copy(::boost::begin(r), ::boost::end(r), t.begin());
    return t;
}

} // namespace suzerain

#endif /* SUZERAIN_SHARED_RANGE_HPP */
