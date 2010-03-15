///////////////////////////////////////////////////////////////////////////////
//
// StlAllocatorTestContainer.hpp
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
// Container tests for allocators. Typical tests involve:
//
//      Adding an element to the container (invokes allocate and construct)
//      Removing an element from the container (invokes destroy and deallocate)
//      Reserving elements (invokes allocate)
//      Swapping containers (invokes deallocate and destroy)
//      Extracting allocator object (invokes container get_allocator)

#if !defined( STL_ALLOCATOR_TEST_CONTAINER_HPP )
#define STL_ALLOCATOR_TEST_CONTAINER_HPP

#include <vector>
#include <deque>
#include <list>
#include <set>
#include <map>
#include <string>
#include <stack>
#include <queue>

#include "StlAllocatorTestOutput.hpp"

namespace StlAllocatorTestbed
{

//-----------------------------------------------------------------------------

template< typename Allocator, typename Container >
void TestVector( const Allocator& a, Container& c )
{
    c.push_back( typename Allocator::value_type() );
    c.clear();
    c.push_back( typename Allocator::value_type() );
    c.reserve( 100 );
    VERIFY( c[0] == typename Allocator::value_type() );
    c.clear();
    USED( a );
    USED( c );
}

//-----------------------------------------------------------------------------

template< typename Allocator, typename Container >
void TestDeque( const Allocator& a, Container& c )
{
    c.push_back( typename Allocator::value_type() );
    c.clear();
    c.assign( 100, typename Allocator::value_type() );
    VERIFY( c.front() == typename Allocator::value_type() );
    c.clear();
    USED( a );
    USED( c );
}

//-----------------------------------------------------------------------------

template< typename Allocator, typename Container >
void TestList( const Allocator& a, Container& c )
{
    c.push_back( typename Allocator::value_type() );
    c.clear();
    c.insert( c.begin(), 100, typename Allocator::value_type() );
    VERIFY( c.front() == typename Allocator::value_type() );
    c.clear();
    USED( a );
    USED( c );
}

//-----------------------------------------------------------------------------

template< typename Allocator, typename Container >
void TestSet( const Allocator& a, Container& c )
{
    c.insert( typename Allocator::value_type() );
    c.clear();
    for( int i = 0; i < 100; ++i )
        c.insert( typename Allocator::value_type( i ) );
    c.insert( typename Allocator::value_type( 0 ) );
    VERIFY( c.find( typename Allocator::value_type( 0 ) ) == c.begin() );
    c.clear();
    USED( a );
    USED( c );
}

//-----------------------------------------------------------------------------

template< typename Allocator, typename Container >
void TestMap( const Allocator& a, Container& c )
{
    c.insert( std::make_pair( typename Allocator::value_type(), 1 ) );
    c.clear();
    for( int i = 0; i < 100; ++i )
        c.insert( std::make_pair( typename Allocator::value_type( i ), i ) );
    c.insert( std::make_pair( typename Allocator::value_type( 0 ), 0 ) );
    VERIFY( c.find( typename Allocator::value_type( 0 ) ) == c.begin() );
    c.clear();
    USED( a );
    USED( c );
}

//-----------------------------------------------------------------------------

template< typename Allocator, typename Container >
void TestString( const Allocator& a, Container& c )
{
    c.push_back( typename Allocator::value_type() );
    c.clear();
    c.assign( 100, typename Allocator::value_type( 1 ) );
    VERIFY( c.find( typename Allocator::value_type( 1 ) ) == 0 );
    c.clear();
    USED( a );
    USED( c );
}

//-----------------------------------------------------------------------------

template< typename Value, typename Allocator, typename Container >
void TestStack( const Value& v, const Allocator&, Container& c )
{
    c.push( v );
    c.pop();
    for( int i = 0; i < 100; ++i )
        c.push( typename Allocator::value_type( i ) );
#if( !defined(_MSC_VER) || _MSC_VER < 1310 )
    VERIFY( c.top() == typename Allocator::value_type( 99 ) );
#endif
    for( int i = 0; i < 100; ++i )
        c.pop();
    USED( c );
}

#if( !defined(_MSC_VER) || _MSC_VER >= 1310 )

// Specialize for bool, since we can't convert from c.top() to bool with
// vector<bool> (see More Exceptional C++ Item 6)
template< typename Value, typename Allocator, typename Container >
void TestStack( const bool& v, const Allocator&, Container& c )
{
    c.push( v );
    c.pop();
    for( int i = 0; i < 100; ++i )
        c.push( v );
    for( int i = 0; i < 100; ++i )
        c.pop();
    USED( c );
}

#endif

//-----------------------------------------------------------------------------

template< typename Allocator, typename Container >
void TestQueue( const Allocator&, Container& c )
{
    c.push( typename Allocator::value_type() );
    c.pop();
    for( int i = 0; i < 100; ++i )
        c.push( typename Allocator::value_type( i ) );
    VERIFY( c.front() == typename Allocator::value_type( 0 ) );
    for( int i = 0; i < 100; ++i )
        c.pop();
    USED( c );
}

//-----------------------------------------------------------------------------

template< typename Allocator, typename Container >
void TestPriorityQueue( const Allocator&, Container& c )
{
    c.push( typename Allocator::value_type() );
    c.pop();
    for( int i = 0; i < 100; ++i )
        c.push( typename Allocator::value_type( i ) );
    for( int i = 0; i < 100; ++i )
        c.pop();
    USED( c );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestVector( const Allocator& a )
{
    typedef std::vector< typename Allocator::value_type, Allocator > Vector;

    Vector v1( Policy< Allocator >::template GetDefault< Vector >( a ) );
    Vector v2( Policy< Allocator >::template GetCopied< Vector >( a ) );

    TestVector( a, v1 );
    TestVector( a, v2 );

    v1.swap( v2 );

    TestVector( a, v1 );
    TestVector( a, v2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestDeque( const Allocator& a )
{
    typedef std::deque< typename Allocator::value_type, Allocator > Deque;

    Deque d1( Policy< Allocator >::template GetDefault< Deque >( a ) );
    Deque d2( Policy< Allocator >::template GetCopied< Deque >( a ) );

    TestDeque( a, d1 );
    TestDeque( a, d2 );

    d1.swap( d2 );

    TestDeque( a, d1 );
    TestDeque( a, d2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestList( const Allocator& a )
{
    typedef std::list< typename Allocator::value_type, Allocator > List;

    List l1( Policy< Allocator >::template GetDefault< List >( a ) );
    List l2( Policy< Allocator >::template GetCopied< List >( a ) );

    TestList( a, l1 );
    TestList( a, l2 );

    l1.swap( l2 );

    TestList( a, l1 );
    TestList( a, l2 );

    l1.splice( l1.end(), l2 );

    TestList( a, l1 );
    TestList( a, l2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestSet( const Allocator& a )
{
    typedef std::set< typename Allocator::value_type,
                      std::less< typename Allocator::value_type >,
                      Allocator > Set;

    Set s1( Policy< Allocator >::template GetDefaultSet< Set >( a ) );
    Set s2( Policy< Allocator >::template GetCopiedSet< Set >( a ) );

    TestSet( a, s1 );
    TestSet( a, s2 );

    s1.swap( s2 );

    TestSet( a, s1 );
    TestSet( a, s2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestMultiset( const Allocator& a )
{
    typedef std::multiset< typename Allocator::value_type,
                           std::less< typename Allocator::value_type >,
                           Allocator > Multiset;

    Multiset s1( Policy< Allocator >::template GetDefaultSet< Multiset >( a ) );
    Multiset s2( Policy< Allocator >::template GetCopiedSet< Multiset >( a ) );

    TestSet( a, s1 );
    TestSet( a, s2 );

    s1.swap( s2 );

    TestSet( a, s1 );
    TestSet( a, s2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestMap( const Allocator& a )
{
    // map allocators are pair allocators
    typedef std::pair< const typename Allocator::value_type, int > Pair;
    typedef typename Allocator::template rebind< Pair >::other Alloc;

    typedef std::map< typename Allocator::value_type, int,
                      std::less< typename Allocator::value_type >,
                      Alloc > Map;

    Map m1( Policy< Allocator >::template GetDefaultMap< Map >( a ) );
    Map m2( Policy< Allocator >::template GetCopiedMap< Map >( a ) );

    TestMap( a, m1 );
    TestMap( a, m2 );

    m1.swap( m2 );

    TestMap( a, m1 );
    TestMap( a, m2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestMultimap( const Allocator& a )
{
    // map allocators are pair allocators
    typedef std::pair< const typename Allocator::value_type, int > Pair;
    typedef typename Allocator::template rebind< Pair >::other Alloc;

    typedef std::multimap< typename Allocator::value_type, int,
                           std::less< typename Allocator::value_type >,
                           Alloc > Multimap;

    Multimap m1( Policy< Allocator >::template GetDefaultMap< Multimap >( a ) );
    Multimap m2( Policy< Allocator >::template GetCopiedMap< Multimap >( a ) );

    TestMap( a, m1 );
    TestMap( a, m2 );

    m1.swap( m2 );

    TestMap( a, m1 );
    TestMap( a, m2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestString( const Allocator& aIn )
{
    // Strings are made of chars or wide chars
    typedef typename Allocator::template rebind< char >::other CharAlloc;
    typedef typename Allocator::template rebind< wchar_t >::other WCharAlloc;

    typedef std::basic_string< char,
                               std::char_traits< char >,
                               CharAlloc > String;
    typedef std::basic_string< wchar_t,
                               std::char_traits< wchar_t >,
                               WCharAlloc > WString;

    // string
    CharAlloc a( aIn );

    String s1( Policy< Allocator >::template GetDefault< String >( a ) );
    String s2( Policy< Allocator >::template GetCopied< String >( a ) );

    TestString( a, s1 );
    TestString( a, s2 );

    s1.swap( s2 );

    TestString( a, s1 );
    TestString( a, s2 );

    // wstring
    WCharAlloc wa( aIn );
    WString w1( Policy< Allocator >::template GetDefault< WString >( wa ) );
    WString w2( Policy< Allocator >::template GetCopied< WString >( wa ) );

    TestString( wa, w1 );
    TestString( wa, w2 );

    w1.swap( w2 );

    TestString( wa, w1 );
    TestString( wa, w2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestStack( const Allocator& a )
{
    // A stack can be made out of deque or vector
    typedef std::deque< typename Allocator::value_type, Allocator > Deque;
    typedef std::vector< typename Allocator::value_type, Allocator > Vector;

    typedef std::stack< typename Allocator::value_type, Deque > DStack;
    typedef std::stack< typename Allocator::value_type, Vector > VStack;

    DStack d1( Policy< Allocator >::template GetDefaultAdapter< DStack, Deque >( a ) );
    DStack d2( Policy< Allocator >::template GetCopiedAdapter< DStack, Deque >( a ) );
    VStack v1( Policy< Allocator >::template GetDefaultAdapter< VStack, Vector >( a ) );
    VStack v2( Policy< Allocator >::template GetCopiedAdapter< VStack, Vector >( a ) );

    typename Allocator::value_type v;

    // Test deque stacks
    TestStack< typename Allocator::value_type >( v, a, d1 );
    TestStack< typename Allocator::value_type >( v, a, d2 );

    // Test vector stacks
    TestStack< typename Allocator::value_type >( v, a, v1 );
    TestStack< typename Allocator::value_type >( v, a, v2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestQueue( const Allocator& a )
{
    // A queue can be made out of deque or list
    typedef std::deque< typename Allocator::value_type, Allocator > Deque;
    typedef std::list< typename Allocator::value_type, Allocator > List;

    typedef std::queue< typename Allocator::value_type, Deque > DQueue;
    typedef std::queue< typename Allocator::value_type, List > LQueue;

    DQueue d1( Policy< Allocator >::template GetDefaultAdapter< DQueue, Deque >( a ) );
    DQueue d2( Policy< Allocator >::template GetCopiedAdapter< DQueue, Deque >( a ) );
    LQueue l1( Policy< Allocator >::template GetDefaultAdapter< LQueue, List >( a ) );
    LQueue l2( Policy< Allocator >::template GetCopiedAdapter< LQueue, List >( a ) );

    // Test deque queues
    TestQueue( a, d1 );
    TestQueue( a, d2 );

    // Test list queues
    TestQueue( a, l1 );
    TestQueue( a, l2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestPriorityQueue( const Allocator& a )
{
    // A priority queue can be made out of vector or deque
    typedef std::vector< typename Allocator::value_type, Allocator > Vector;
    typedef std::deque< typename Allocator::value_type, Allocator > Deque;

    typedef std::priority_queue< typename Allocator::value_type, Vector > VQueue;
    typedef std::priority_queue< typename Allocator::value_type, Deque > DQueue;

    VQueue v1( Policy< Allocator >::template GetDefaultPriorityQueue< VQueue, Vector >( a ) );
    VQueue v2( Policy< Allocator >::template GetCopiedPriorityQueue< VQueue, Vector >( a ) );
    DQueue d1( Policy< Allocator >::template GetDefaultPriorityQueue< DQueue, Deque >( a ) );
    DQueue d2( Policy< Allocator >::template GetCopiedPriorityQueue< DQueue, Deque >( a ) );

    // Test vector queues
    TestPriorityQueue( a, v1 );
    TestPriorityQueue( a, v2 );

    // Test deque queues
    TestPriorityQueue( a, d1 );
    TestPriorityQueue( a, d2 );
}

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestWithContainers( const Allocator& a )
{
    TestVector( a );
    TestDeque( a );
    TestList( a );
    TestSet( a );
    TestMultiset( a );
    TestMap( a );
    TestMultimap( a );
    TestString( a );
    TestStack( a );
    TestQueue( a );
    TestPriorityQueue( a );
}

} // end StlAllocatorTestbed namespace

#endif // STL_ALLOCATOR_TEST_CONTAINER_HPP
