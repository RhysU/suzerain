//--------------------------------------------------------------------------
//
// Code adapted from Lai Shiaw San Kent's article "C++ Standard
// Allocator, An Introduction and Implementation" from The Code Project
// (http://www.codeproject.com/KB/cpp/allocator.aspx).
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// except where that copyright conflicts with the rights of
// Lai Shiaw San Kent and Sean Kent.
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_ALLOCATOR_HPP
#define SUZERAIN_ALLOCATOR_HPP

/** @file
 * Allocator implementations, including aligned allocators
 */

#include <limits>
#include <boost/utility.hpp>

// FIXME Policy-based framework lacks allocator<void> specializations
// See http://msdn.microsoft.com/en-us/library/w3b7688x.aspx for information.

namespace suzerain {

/**
 * A generic set of object traits for the policy driven allocator.  Necessary
 * since some classes may overload <tt>operator&</tt>.  Classes which overload
 * this operator should provide a specialization of allocator_object_traits.
 *
 * @see Lai Shiaw San Kent's article <a
 * href="http://www.codeproject.com/KB/cpp/allocator.aspx">C++ Standard
 * Allocator, An Introduction and Implementation</a> for more information.
 */
template<typename T>
class allocator_object_traits {
public:
    /** Convert an allocator_object_traits<T> to allocator_object_traits<U> */
    template<typename U>
    struct rebind {
        typedef allocator_object_traits<U> other;
    };

    // Make explicit all constructors, including copy constructors
    explicit allocator_object_traits() {}
    ~allocator_object_traits() {}
    explicit allocator_object_traits(allocator_object_traits  const&) {}
    template <typename U>
    explicit allocator_object_traits(allocator_object_traits<U> const&) {}

    /** Address-related methods */
    //@{
    T* address(T& r) { return boost::addressof(r); }
    T const* address(T const& r) { return boost::addressof(r); }
    //@}

    /** Copy construction and destruction helpers */
    //@{
    static void construct(T* p, const T& t) { new(p) T(t); }
    static void destroy(T* p) { p->~T(); }
    //@}
};

/**
 * A free-store based allocation policy.  Uses global <tt>operator new</tt> and
 * global <tt>operator delete</tt> to handle allocate and deallocate
 * invocations.
 *
 * @see Lai Shiaw San Kent's article <a
 * href="http://www.codeproject.com/KB/cpp/allocator.aspx">C++ Standard
 * Allocator, An Introduction and Implementation</a> for more information.
 */
template<typename T>
class allocator_freestore_policy {
public:
    // Typedefs following std::allocator contract
    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    /** Rebind to obtain allocator_freestore_policy<U> */
    template<typename U>
    struct rebind {
        typedef allocator_freestore_policy<U> other;
    };

    // Make explicit all constructors, including copy constructors
    explicit allocator_freestore_policy() {}
    ~allocator_freestore_policy() {}
    explicit allocator_freestore_policy(allocator_freestore_policy const&) {}
    template <typename U>
    explicit allocator_freestore_policy(allocator_freestore_policy<U> const&) {}

    /** Free-store based memory allocation */
    //@{
    pointer allocate(size_type cnt,
                     typename std::allocator<void>::const_pointer = 0)
    {
        return reinterpret_cast<pointer>(::operator new(cnt * sizeof (T)));
    }
    void deallocate(pointer p, size_type) { ::operator delete(p); }
    //@}

    /** max_size method following std::allocator contract */
    size_type max_size() const
    {
        return std::numeric_limits<size_type>::max() / sizeof(T);
    }
};

/** Provide information about compatibility between different allocators. */
//@{
/**
 * Memory can be deallocated from other allocator_freestore_policy
 * instantiations.
 **/
template<typename T, typename T2>
inline bool operator==(allocator_freestore_policy<T> const&,
                       allocator_freestore_policy<T2> const&)
{
    return true;
}

/** Memory cannot be deallocated by other allocator types */
template<typename T, typename OtherAllocator>
inline bool operator==(allocator_freestore_policy<T> const&,
                       OtherAllocator const&)
{
    return false;
}
//@}

/**
 * A policy-driven allocator implementation.
 *
 * @see Lai Shiaw San Kent's article <a
 * href="http://www.codeproject.com/KB/cpp/allocator.aspx">C++ Standard
 * Allocator, An Introduction and Implementation</a> for more information.
 */
template<typename T,
         typename Policy = allocator_freestore_policy<T>,
         typename Traits = allocator_object_traits<T> >
class allocator : public Policy, public Traits {
public:
    // Typedefs required per std::allocator contract
    typedef typename Policy::size_type size_type;
    typedef typename Policy::difference_type difference_type;
    typedef typename Policy::pointer pointer;
    typedef typename Policy::const_pointer const_pointer;
    typedef typename Policy::reference reference;
    typedef typename Policy::const_reference const_reference;
    typedef typename Policy::value_type value_type;

    /** Provides required rebind templated typedef */
    template<typename U>
    struct rebind {
        typedef allocator<U,typename Policy::template rebind<U>::other> other;
    };

    // Make explicit all constructors, including copy constructors
    explicit allocator() {}
    ~allocator() {}
    allocator(allocator const& rhs) : Policy(rhs), Traits(rhs) {}
    template <typename U>
    explicit allocator(allocator<U> const&) {}
    template <typename U, typename P, typename T2>
    allocator(allocator<U, P, T2> const& rhs) : Policy(rhs), Traits(rhs) {}

    /** Allocation and deallocation uses the policy template parameter */
    //@{
    pointer allocate(
        size_type cnt,
        typename std::allocator<void>::const_pointer hint = 0)
    {
        return Policy::allocate(cnt, hint);
    }

    void deallocate(pointer p, size_type cnt)
    {
        Policy::deallocate(p, cnt);
    }
    //@}
};

/**
 * Determines if memory from another allocator
 * can be deallocated from this one
 */
//@{
template<typename T, typename P, typename Tr>
inline bool operator==(allocator<T, P, Tr> const& lhs,
                       allocator<T, P, Tr> const& rhs)
{
    return operator==(static_cast<P const&>(lhs), static_cast<P const&>(rhs));
}

template<typename T,  typename P,  typename Tr,
         typename T2, typename P2, typename Tr2>
inline bool operator==(allocator<T,  P,  Tr>  const& lhs,
                       allocator<T2, P2, Tr2> const& rhs)
{
    return operator==(static_cast<P const&>(lhs), static_cast<P2 const&>(rhs));
}

template<typename T, typename P, typename Tr, typename OtherAllocator>
inline bool operator==(allocator<T, P, Tr> const& lhs,
                       OtherAllocator const& rhs)
{
    return operator==(static_cast<P const&>(lhs), rhs);
}

template<typename T, typename P, typename Tr>
inline bool operator!=(allocator<T, P, Tr> const& lhs,
                       allocator<T, P, Tr> const& rhs)
{
    return !operator==(lhs, rhs);
}

template<typename T,  typename P,  typename Tr,
         typename T2, typename P2, typename Tr2>
inline bool operator!=(allocator<T, P, Tr> const& lhs,
                       allocator<T2, P2, Tr2> const& rhs)
{
    return !operator==(lhs, rhs);
}

template<typename T, typename P, typename Tr, typename OtherAllocator>
inline bool operator!=(allocator<T, P, Tr> const& lhs,
                       OtherAllocator const& rhs)
{
    return !operator==(lhs, rhs);
}
//@}

} // namespace suzerain

#endif // SUZERAIN_ALLOCATOR_HPP
