///////////////////////////////////////////////////////////////////////////////
//
// StlAllocatorTestMembers.hpp
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
//    Tests all required allocator member functions.

#ifndef STL_ALLOCATOR_TEST_MEMBERS_HPP
#define STL_ALLOCATOR_TEST_MEMBERS_HPP

#include "StlAllocatorTestOutput.hpp"
#include <cstring>

namespace StlAllocatorTestbed
{

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestMemberFunctionPresence( const Allocator& a )
{
    // Typedefs
    typename Allocator::size_type       st;
    typename Allocator::difference_type dt;
    typename Allocator::pointer         p;
    typename Allocator::const_pointer   cp;
    typename Allocator::value_type      v;
    typename Allocator::reference       r = v;   // Note: must be initialized
    typename Allocator::const_reference cr = v;  // Note: must be initialized

    // allocator()
    Allocator b( Policy< Allocator >::Construct( a ) );
    USED( b );

    // allocator( const allocator& )
    Allocator c(a);

    // allocator( const allocator<U>& ) and rebind
    typename Allocator::template rebind< C >::other d( c );

    // Although operator=() may be implemented by an allocator, it is
    // not required by the C++ Standard, nor is it to be used by
    // container implementations.

    // address( reference )
    p = c.address( r );

    // address( const_reference )
    cp = c.address( cr );

    // difference_type
    dt = p - cp;

    // max_size()
    st = c.max_size();

    // allocate( size_type )
    p = c.allocate( 1 );

    // deallocate( pointer, size_type )
    c.deallocate( p, 1 );

    // allocate( size_type, allocator<void>::const_pointer )
    p = c.allocate( 1, NULL );

    // construct( pointer, const T& )
    c.construct( p, v );

    // destroy( pointer )
    c.destroy( p );

    // clean up
    c.deallocate( p, 1 );

    // operator ==( const A&, const A& )
    bool f1 = ( a == c );

    // operator ==( const A&, const B& )
    bool f2 = ( a == d );

    // operator !=( const A&, const A& )
    bool f3 = ( a != c );

    // operator !=( const A&, const B& )
    bool f4 = ( a != d );

    // Force the compiler to believe that all variables are important
    // so it doesn't optimize away any code
    USED( a );
    USED( c );
    USED( d );

    USED( f1 );
    USED( f2 );
    USED( f3 );
    USED( f4 );

    USED( st );
    USED( dt );
    USED( p  );
    USED( cp );
    USED( v  );
    USED( r  );
    USED( cr );

    // Implicit ~allocator() happens here
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestMemberFunctions( const Allocator& a )
{
    // Typedefs
    typename Allocator::pointer         p = NULL;
    typename Allocator::value_type      v;
    typename Allocator::reference       r = v;   // Note: must be initialized
    typename Allocator::const_reference cr = v;  // Note: must be initialized

    // allocator()
    Allocator b( Policy< Allocator >::Construct( a ) );
    VERIFY( b == Policy< Allocator >::Construct( a ) );

    // allocator( const allocator& )
    // The C++ Standard requires that Alloc<T> a(b) must give the post-condition
    // that Alloc<U>(a) == b
    Allocator c(a);
    VERIFY( c == a );

    // allocator( const allocator<U>& ) and rebind
    typename Allocator::template rebind< C >::other d( c );
    typename Allocator::template rebind< typename Allocator::value_type >::other e( d );
    VERIFY( c == e );

    // address( reference )
    VERIFY( c.address( r ) == &v );

    // address( const_reference )
    VERIFY( c.address( cr ) == &v );

    // max_size()
    VERIFY( c.max_size() > 0 );

    // allocate( size_type )
    try
    {
        p = c.allocate( 1 );
    }
    catch(...)
    {
        BOOST_FAIL( "Allocation Failure" );
    }
    VERIFY( p != NULL );

    // Write to allocated memory
    memset( p, 'x', sizeof( typename Allocator::value_type ) );

    // deallocate( pointer, size_type )
    c.deallocate( p, 1 );

    // allocate( size_type, allocator<void>::const_pointer )
    try
    {
        p = c.allocate( 1, NULL );
    }
    catch(...)
    {
        BOOST_FAIL( "Allocation Failure" );
    }
    VERIFY( p != NULL );

    // Write to allocated memory
    memset( p, 'x', sizeof( typename Allocator::value_type ) );

    // construct( pointer, const T& )
    c.construct( p, v );
    VERIFY( *p == v );

    // destroy( pointer )
    c.destroy( p );

    // clean up
    c.deallocate( p, 1 );

    // allocate multiple items
    try
    {
        p = c.allocate( 2, NULL );
    }
    catch(...)
    {
        BOOST_FAIL( "Allocation Failure" );
    }
    VERIFY( p != NULL );

    // Write to allocated memory
    memset( p, 'x', sizeof( typename Allocator::value_type ) * 2 );

    // construct( pointer, const T& )
    c.construct( p, v );
    VERIFY( *p == v );

    c.construct( p + 1, v );
    VERIFY( *(p + 1) == v );

    // destroy( pointer )
    c.destroy( p );
    c.destroy( p + 1 );

    // clean up
    c.deallocate( p, 2 );

    // operator ==( const A&, const A& )
    VERIFY( a == c );

    // operator ==( const A&, const B& )
    bool f1 = ( a == d );

    // operator !=( const A&, const A& )
    VERIFY( !( a != c ) );

    // operator !=( const A&, const B& )
    bool f2 = ( a != d );

    // Force the compiler to believe that all variables are important
    // so it doesn't optimize away any code
    USED( a );
    USED( c );
    USED( d );

    USED( f1 );
    USED( f2 );

    USED( p  );
    USED( v  );
    USED( r  );
    USED( cr );

    // Implicit ~allocator() happens here
}

} // end StlAllocatorTestbed namespace

#endif // STL_ALLOCATOR_TEST_MEMBERS_HPP
