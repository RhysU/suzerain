///////////////////////////////////////////////////////////////////////////////
//
// StlAllocatorTestTypes.hpp
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
//    Types used for testing allocators. Includes some simple types as well as
//    objects designed to put maximum stress on allocators. All types can be
//    constructed from integers, which eases testing containers that only support
//    unique elements, like set and map.

#ifndef STL_ALLOCATOR_TEST_TYPES_HPP
#define STL_ALLOCATOR_TEST_TYPES_HPP

#include <string>
#include <complex>
#include "StlAllocatorTestOutput.hpp"

#pragma warning(push,disable:2017)
namespace StlAllocatorTestbed
{

//-----------------------------------------------------------------------------

// This test class contains plain old data types (PODs) only

class C
{
    char   c;
    short  s;
    int    i;
    long   l;
    float  f;
    double d;

public:

    C()
    :
        c( '0' ),
        s( 0 ),
        i( 0 ),
        l( 0L ),
        f( 0.0f ),
        d( 0.0 )
    {
        Verify( *this );
    }

    C( const C& x )
    :
        c( x.c ),
        s( x.s ),
        i( x.i ),
        l( x.l ),
        f( x.f ),
        d( x.d )
    {
        Verify( *this );
    }

    explicit C( int in )
    :
        c( '0' ),
        s( 0 ),
        i( in ),
        l( 0L ),
        f( 0.0f ),
        d( 0.0 )
    {
        Verify( *this );
    }

    C& operator=( const C& x )
    {
        Verify( x );
        Verify( *this );

        C temp( x ); // Exceptional C++ copy assignment swap trick; Item 13
        std::swap( c, temp.c );
        std::swap( s, temp.s );
        std::swap( i, temp.i );
        std::swap( l, temp.l );
        std::swap( f, temp.f );
        std::swap( d, temp.d );

        return *this;
    }

    ~C()
    {
        Verify( *this );
    }

    // Comparable so it can be used in sets and maps
    bool operator<( const C& x ) const
    {
        Verify( x );
        Verify( *this );
        return i < x.i;
    }

    // Hashable so it can be used in hash containers
    operator size_t() const
    {
        Verify( *this );
        return size_t( &c );
    }

    bool operator==( const C& x ) const
    {
#pragma warning(push,disable:1572)
        return c == x.c &&
               s == x.s &&
               i == x.i &&
               l == x.l &&
               f == x.f &&
               d == x.d;
#pragma warning(pop)
    }

    void Verify( const C& x ) const
    {
#pragma warning(push,disable:1572)
        VERIFY( x.c == '0'  );
        VERIFY( x.s == 0    );
        VERIFY( x.i >= 0    );  // assumes outside callers never exceed 100
        VERIFY( x.i <= 100  );
        VERIFY( x.l == 0L   );
        VERIFY( x.f == 0.0f );
        VERIFY( x.d == 0.0  );
#pragma warning(pop)
    }

};

//-----------------------------------------------------------------------------

// This test class has non-trival constructors, assignment op and destructor

class D
{
    char* p;

public:

    D()
    :
        p( new char( 'p' ) )
    {
        Verify( *this );
    }

    D( const D& d )
    :
        p( new char( *d.p ) )
    {
        Verify( d );
    }

    explicit D( int i )
    :
#pragma warning(push,disable:2259)
        p( new char( static_cast<char>(i) ) )
#pragma warning(pop)
    {
        Verify( *this );
    }

    D& operator=( const D& d )
    {
        Verify( d );
        Verify( *this );
        D temp( d ); // Exceptional C++ copy assignment swap trick; Item 13
        std::swap( p, temp.p );
        return *this;
    }

    ~D()
    {
        Verify( *this );
        delete p;
    }

    // Comparable so it can be used in sets and maps
    bool operator<( const D& d ) const
    {
        Verify( d );
        Verify( *this );
        return *p < *d.p;
    }

    // Hashable so it can be used in hash containers
    operator size_t() const
    {
        Verify( *this );
        return size_t( p );
    }

    bool operator==( const D& d ) const
    {
        return *p == *d.p;
    }

    void Verify( const D& d ) const
    {
        VERIFY( *d.p == 'p' || // assumes outside callers never exceed 100
                ( *d.p >= 0 && *d.p <= 100 ) );
    }

};

//-----------------------------------------------------------------------------

// A collection of a variety of types

class E
{
    std::string         s;
    std::complex<float> c;
    int                 i;
    double              d;

public:

    E()
    :
        s( "0123456789abcxyz" ), // enough to exceed the small string optimization
        c( 0.1f, 2.3f ),
        i( 0x12345678 ),
        d( 123.456789 )
    {
        Verify( *this );
    }

    E( const E& e )
    :
        s( e.s ), // enough to exceed the small string optimization
        c( e.c ),
        i( e.i ),
        d( e.d )
    {
        Verify( *this );
    }

    explicit E( int in )
    :
        s( "0123456789abcxyz" ), // enough to exceed the small string optimization
        c( 0.1f, 2.3f ),
        i( in ),
        d( 123.456789 )
    {
        Verify( *this );
    }

    E& operator=( const E& e )
    {
        Verify( e );
        Verify( *this );
        E temp( e ); // Exceptional C++ copy assignment swap trick; Item 13
        std::swap( s, temp.s );
        std::swap( c, temp.c );
        std::swap( i, temp.i );
        std::swap( d, temp.d );
        return *this;
    }

    ~E()
    {
        Verify( *this );
    }

    // Comparable so it can be used in sets and maps
    bool operator<( const E& e ) const
    {
        Verify( e );
        Verify( *this );
        return i < e.i;
    }

    // Hashable so it can be used in hash containers
    operator size_t() const
    {
        Verify( *this );

        // Use the string pointer, since it's guaranteed unique
        return size_t( s.c_str() );
    }

    bool operator==( const E& e ) const
    {
#pragma warning(push,disable:1572)
        return s == e.s &&
               c == e.c &&
               i == e.i &&
               d == e.d;
#pragma warning(pop)
    }

    void Verify( const E& e ) const
    {
        VERIFY( e.s == "0123456789abcxyz" );
        VERIFY( e.c == std::complex<float>( 0.1f, 2.3f ) );
        VERIFY( e.i == 0x12345678 || // assumes outside callers never exceed 100
                ( e.i >= 0 && e.i <= 100 ) );
#pragma warning(push,disable:1572)
        VERIFY( e.d == 123.456789 );
#pragma warning(pop)
    }

};


} // end StlAllocatorTestbed namespace
#pragma warning(pop)

#endif // STL_ALLOCATOR_TEST_TYPES_HPP
