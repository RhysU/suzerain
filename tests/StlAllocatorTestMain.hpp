///////////////////////////////////////////////////////////////////////////////
//
// StlAllocatorTestMain.hpp
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
//    This module is designed for testing custom allocators. Simply pass your
//    custom allocator to the TestStlAllocator() function, stand back and let
//    the show begin.
//
// Details:
//
//    TestStlAllocator() verifies that the incoming allocator meets the
//    requirements of the C++ Standard. It does this in three ways:
//
//    1) Ensuring that all allocator functions are available. If a
//       function is not available, the code will fail to compile.
//    2) Ensuring that all allocator functions work as advertised.
//       If the functions fail, error messages are sent to standard output.
//    3) Ensuring that the allocator works properly will all standard
//       containers using a variety of stored objects.

#ifndef STL_ALLOCATOR_TEST_MAIN_HPP
#define STL_ALLOCATOR_TEST_MAIN_HPP

#include "StlAllocatorTestTypes.hpp"
#include "StlAllocatorTestMembers.hpp"
#include "StlAllocatorTestContainer.hpp"

namespace StlAllocatorTestbed
{

//-----------------------------------------------------------------------------

template< typename Allocator >
void TestAlloc( const Allocator& a )
{
    // Verify functionality with a variety of types.
    // Conditional complication can selectively reduce footprint.

#if !defined(STL_ALLOCATOR_TEST_PART) || STL_ALLOCATOR_TEST_PART == 0
    StlAllocatorTestbed::TestMemberFunctionPresence( a );

    StlAllocatorTestbed::TestMemberFunctions( a );

    { // tests vector<bool> specialization
        typename Allocator::template rebind< bool >::other b( a );
        StlAllocatorTestbed::TestMemberFunctions( b );
    }

    { // POD types
        typename Allocator::template rebind< C >::other c( a );
        StlAllocatorTestbed::TestMemberFunctions( c );
    }

    { // allocated internal type
        typename Allocator::template rebind< D >::other d( a );
        StlAllocatorTestbed::TestMemberFunctions( d );
    }

    { // complicated type
        typename Allocator::template rebind< E >::other e( a );
        StlAllocatorTestbed::TestMemberFunctions( e );
    }

    { // simple test
        typename Allocator::template rebind< int >::other i( a );
        StlAllocatorTestbed::TestMemberFunctions( i );
    }

    { // void pointers
        typename Allocator::template rebind< void* >::other v( a );
        StlAllocatorTestbed::TestMemberFunctions( v );
    }

    { // regular pointers
        typename Allocator::template rebind< char* >::other p( a );
        StlAllocatorTestbed::TestMemberFunctions( p );
    }

    // Disabled!  See http://msdn.microsoft.com/en-us/library/w3b7688x.aspx
    // void specialization
    // typename Allocator::template rebind< void >::other x( a );

    StlAllocatorTestbed::TestWithContainers( a );
#endif

#if !defined(STL_ALLOCATOR_TEST_PART) || STL_ALLOCATOR_TEST_PART == 1
    { // tests vector<bool> specialization
        typename Allocator::template rebind< bool >::other b( a );
        StlAllocatorTestbed::TestWithContainers( b );
    }
#endif

#if !defined(STL_ALLOCATOR_TEST_PART) || STL_ALLOCATOR_TEST_PART == 2
    { // POD types
        typename Allocator::template rebind< C >::other c( a );
        StlAllocatorTestbed::TestWithContainers( c );
    }
#endif

#if !defined(STL_ALLOCATOR_TEST_PART) || STL_ALLOCATOR_TEST_PART == 3
    { // allocated internal type
        typename Allocator::template rebind< D >::other d( a );
        StlAllocatorTestbed::TestWithContainers( d );
    }
#endif

#if !defined(STL_ALLOCATOR_TEST_PART) || STL_ALLOCATOR_TEST_PART == 4
    { // complicated type
        typename Allocator::template rebind< E >::other e( a );
        StlAllocatorTestbed::TestWithContainers( e );
    }
#endif

#if !defined(STL_ALLOCATOR_TEST_PART) || STL_ALLOCATOR_TEST_PART == 5
    { // simple test
        typename Allocator::template rebind< int >::other i( a );
        StlAllocatorTestbed::TestWithContainers( i );
    }
#endif

#if !defined(STL_ALLOCATOR_TEST_PART) || STL_ALLOCATOR_TEST_PART == 6
    { // void pointers
        typename Allocator::template rebind< void* >::other v( a );
        StlAllocatorTestbed::TestWithContainers( v );
    }
#endif

#if !defined(STL_ALLOCATOR_TEST_PART) || STL_ALLOCATOR_TEST_PART == 7
    { // regular pointers
        typename Allocator::template rebind< char* >::other p( a );
        StlAllocatorTestbed::TestWithContainers( p );
    }
#endif

}

} // end StlAllocatorTestbed namespace

#endif // STL_ALLOCATOR_TEST_MAIN_HPP
