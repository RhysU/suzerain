///////////////////////////////////////////////////////////////////////////////
//
// StlAllocatorTestPolicy.hpp
//
// Copyright (C) 2003 Pete Isensee (PKIsensee@msn.com).
// All rights reserved worldwide.
//
// Permission to copy, modify, reproduce or redistribute this source code is
// granted provided the above copyright notice is retained in the resulting
// source code.
//
// This software is provided "as is" and without any express or implied
// warranties.
//
// Description:
//
//    Where allocator testing policy is defined

#ifndef STL_ALLOCATOR_TEST_POLICY_HPP
#define STL_ALLOCATOR_TEST_POLICY_HPP

#include <functional>
#include <utility>

namespace StlAllocatorTestbed
{

//-----------------------------------------------------------------------------

// The following template policies allow certain portions of the testbed
// to be enabled or disabled. For instance, most allocators can be default
// constructed. However, some custom allocators cannot be default constructed.

template< typename Allocator >
struct Policy
{
    // Construct an allocator using the default constructor
    // Parameter is not used, but required for specialization
    static Allocator Construct( const Allocator& )
    {
        return Allocator();
    };

    // Construct a container using the container default constructor
    // Parameter is not used, but required for specialization
    template< typename Container >
    static Container GetDefault( const Allocator& )
    {
        return Container();
    };

    // Construct a set using the default set constructor
    // Separate function since specializations may be required
    template< typename Set >
    static Set GetDefaultSet( const Allocator& )
    {
        return Set();
    };

    // Construct a map using the default map constructor
    // Separate function since specializations may be required
    template< typename Map >
    static Map GetDefaultMap( const Allocator& )
    {
        return Map();
    };

    // Construct an adapter (stack, queue) using the default adapter constructor
    // Separate function since specializations may be required
    template< typename Adapter, typename Container >
    static Adapter GetDefaultAdapter( const Allocator& )
    {
        return Adapter();
    }

    // Construct a priority queue using the default priority queue constructor
    // Separate function since specializations may be required
    template< typename PriorityQueue, typename Container >
    static PriorityQueue GetDefaultPriorityQueue( const Allocator& )
    {
        return PriorityQueue();
    }

    // Construct a container by specifying the allocator
    template< typename Container >
    static Container GetCopied( const Allocator& a )
    {
        return Container( a );
    }

    // Construct a set by specifying the allocator
    template< typename Set >
    static Set GetCopiedSet( const Allocator& a )
    {
        // Sets require a comparison functor prior to the allocator
        // object in the container ctor
        return Set( std::less< typename Allocator::value_type >(), a );
    }

    // Construct a map by specifying the allocator
    template< typename Map >
    static Map GetCopiedMap( const Allocator& aIn )
    {
        // map allocators are pair allocators
        typedef const typename Allocator::value_type const_value_type;
        typedef std::pair<const_value_type,int> Pair;
        typedef typename Allocator::template rebind< Pair >::other Alloc;
        Alloc a( aIn );

        // Maps require a comparison functor prior to the allocator
        // object in the map ctor.
        return Map( std::less< typename Allocator::value_type >(), a );
    }

    // Construct an adapter (stack, queue) by specifying the allocator
    template< typename Adapter, typename Container >
    static Adapter GetCopiedAdapter( const Allocator& a )
    {
        // Adapters take the container in the constructor
        return Adapter( Container( a ) );
    }

    // Construct a priority queue by specifying the allocator
    template< typename PriorityQueue, typename Container >
    static PriorityQueue GetCopiedPriorityQueue( const Allocator& a )
    {
        // Priority queues require a comparison functor prior to the
        // allocator object in the container ctor.
        return PriorityQueue( std::less< typename Allocator::value_type >(), Container( a ) );
    }

};

// Any allocator that requires constructor parameters can use this macro
// immediately following the inclusion of this header to disable the
// portion of the testbed that constructs allocators using the default
// constructor. Example: SET_DEFAULT_CONSTRUCTABLE_OFF( MyAllocator )

#define SET_DEFAULT_CONSTRUCTABLE_OFF( x )              \
template< typename T >                                  \
struct StlAllocatorTestbed::Policy< x<T> >              \
{                                                       \
    static x<T> Construct( const x<T>& a )              \
    {                                                   \
        return x<T>(a);                                 \
    };                                                  \
                                                        \
    template< typename Container >                      \
    static Container GetDefault( const x<T>& a )        \
    {                                                   \
        return GetCopied< Container >( a );             \
    };                                                  \
                                                        \
    template< typename Set >                            \
    static Set GetDefaultSet( const x<T>& a )           \
    {                                                   \
        return GetCopiedSet< Set >( a );                \
    };                                                  \
                                                        \
    template< typename Map >                            \
    static Map GetDefaultMap( const x<T>& a )           \
    {                                                   \
        return GetCopiedMap< Map >( a );                \
    };                                                  \
                                                        \
    template< typename Adapter, typename Container >    \
    static Adapter GetDefaultAdapter( const x<T>& a )   \
    {                                                   \
        return GetCopiedAdapter< Adapter, Container >( a ); \
    }                                                   \
                                                        \
    template< typename PriorityQueue, typename Container > \
    static PriorityQueue GetDefaultPriorityQueue( const x<T>& a ) \
    {                                                   \
        return GetCopiedPriorityQueue< PriorityQueue, Container >( a ); \
    }                                                   \
                                                        \
    template< typename Container >                      \
    static Container GetCopied( const x<T>& a )         \
    {                                                   \
        return Container( a );                          \
    }                                                   \
                                                        \
    template< typename Set >                            \
    static Set GetCopiedSet( const x<T>& a )            \
    {                                                   \
        return Set( std::less< T >(), a );              \
    }                                                   \
                                                        \
    template< typename Map >                            \
    static Map GetCopiedMap( const x<T>& aIn )          \
    {                                                   \
        typedef typename std::pair< const T, int > Pair;         \
        typedef typename x<T>::template rebind< Pair >::other Alloc; \
        Alloc a( aIn );                                 \
        return Map( std::less< T >(), a );              \
    }                                                   \
                                                        \
    template< typename Adapter, typename Container >    \
    static Adapter GetCopiedAdapter( const x<T>& a )    \
    {                                                   \
        return Adapter( Container( a ) );               \
    }                                                   \
                                                        \
    template< typename PriorityQueue, typename Container > \
    static PriorityQueue GetCopiedPriorityQueue( const x<T>& a ) \
    {                                                   \
        return PriorityQueue( std::less< T >(), Container( a ) ); \
    }                                                   \
};

} // end StlAllocatorTestbed namespace

#endif // STL_ALLOCATOR_TEST_POLICY_HPP
