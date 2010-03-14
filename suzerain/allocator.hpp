//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
// except where that copyright may conflict with Lai Shiaw San Kent and Sean
// Kent.  Code adapted from Lai Shiaw San Kent's article "C++ Standard
// Allocator, An Introduction and Implementation" from The Code Project
// (http://www.codeproject.com/KB/cpp/allocator.aspx).
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// allocator.hpp: STL Allocator implementations, including aligned allocators
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef __SUZERAIN_ALLOCATOR_HPP
#define __SUZERAIN_ALLOCATOR_HPP

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
    T* address(T& r) { return &r; }
    T const* address(T const& r) { return &r; }
    //@}

    /** Copy construction and destruction helpers */
    //@{
    static void construct(T* p, const T& t) { new(p) T(t); }
    static void destroy(T* p) { p->~T(); }
    //@}
};

/**
 * A free-store based allocation policy.  Uses <tt>::operator new</tt> and
 * <tt>::operator delete</tt> to handle allocate and deallocate invocations.
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
        typedef allocator<U,typename Policy::rebind<U>::other> other;
    };

    // Make explicit all constructors, including copy constructors
    explicit allocator() {}
    ~allocator() {}
    allocator(allocator const& rhs):Traits(rhs), Policy(rhs) {}
    template <typename U>
    explicit allocator(allocator<U> const&) {}
    template <typename U, typename P, typename T2>
    allocator(allocator<U, P, T2> const& rhs) : Traits(rhs), Policy(rhs) {}

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
    return operator==(static_cast<P&>(lhs), static_cast<P&>(rhs));
}

template<typename T,  typename P,  typename Tr,
         typename T2, typename P2, typename Tr2>
inline bool operator==(allocator<T,  P,  Tr>  const& lhs,
                       allocator<T2, P2, Tr2> const& rhs)
{
    return operator==(static_cast<P&>(lhs), static_cast<P2&>(rhs));
}

template<typename T, typename P, typename Tr, typename OtherAllocator>
inline bool operator==(allocator<T, P, Tr> const& lhs,
                       OtherAllocator const& rhs)
{
    return operator==(static_cast<P&>(lhs), rhs);
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

#endif // __SUZERAIN_ALLOCATOR_HPP
