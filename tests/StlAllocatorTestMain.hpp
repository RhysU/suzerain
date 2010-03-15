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
    // Verify functionality with a variety of types
    typename Allocator::template rebind< bool >::other   b( a ); // tests vector<bool> specialization
    typename Allocator::template rebind< C >::other      c( a ); // POD types
    typename Allocator::template rebind< D >::other      d( a ); // allocated internal type
    typename Allocator::template rebind< E >::other      e( a ); // complicated type
    typename Allocator::template rebind< int >::other    i( a ); // simple test
    typename Allocator::template rebind< void* >::other  v( a ); // void pointers
    typename Allocator::template rebind< char* >::other  p( a ); // regular pointers
    // FIXME Disabled!  See http://msdn.microsoft.com/en-us/library/w3b7688x.aspx
    // typename Allocator::template rebind< void >::other   x( a ); // must provide this specialization

    TestMemberFunctionPresence( a ); // in StlAllocatorTestMembers.hpp

    TestMemberFunctions( a ); // in StlAllocatorTestMembers.hpp
    TestMemberFunctions( b );
    TestMemberFunctions( c );
    TestMemberFunctions( d );
    TestMemberFunctions( e );
    TestMemberFunctions( i );
    TestMemberFunctions( v );
    TestMemberFunctions( p );

    TestWithContainers( a ); // in StlAllocatorTestContainer.hpp
    TestWithContainers( b );
    TestWithContainers( c );
    TestWithContainers( d );
    TestWithContainers( e );
    TestWithContainers( i );
    TestWithContainers( v );
    TestWithContainers( p );
}

} // end StlAllocatorTestbed namespace

#endif // STL_ALLOCATOR_TEST_MAIN_HPP
