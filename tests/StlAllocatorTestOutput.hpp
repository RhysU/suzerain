///////////////////////////////////////////////////////////////////////////////
//
// StlAllocatorTestOutput.hpp
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
//    Error and output functions

#ifndef STL_ALLOCATOR_TEST_OUTPUT_HPP
#define STL_ALLOCATOR_TEST_OUTPUT_HPP

//-----------------------------------------------------------------------------

// Forces the compiler to believe that 'a' cannot be optimized away
#define USED(a) ::StlAllocatorTestbed::Used( &(a) )

// Verifies the expression using the Boost.Test framework
#define VERIFY(b) BOOST_CHECK(b);

namespace StlAllocatorTestbed
{

//-----------------------------------------------------------------------------

template< typename T >
const T* Used( T *t )
{
    return t;
}

} // end StlAllocatorTestbed

#endif // STL_ALLOCATOR_TEST_OUTPUT_HPP
