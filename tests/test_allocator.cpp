///////////////////////////////////////////////////////////////////////////////
//
// StlAllocatorTestMain.cpp
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
//-----------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/allocator.hpp>
#include <boost/test/included/unit_test.hpp>

#pragma warning(disable:383 1418)
#include "StlAllocatorTestPolicy.hpp" // The default policy for testing

// Use DEFAULT_CONSTRUCTABLE_OFF if the allocator cannot be default constructed
// SET_DEFAULT_CONSTRUCTABLE_OFF( MyAlloc )

#include "StlAllocatorTestMain.hpp" // Allocator test code


// Testing std::allocator increases confidence Allocator tests are well-formed
BOOST_AUTO_TEST_SUITE( std_allocator )

BOOST_AUTO_TEST_CASE( allocator_test_bed )
{
    // Although some allocator type must be specified (like int), the testing
    // functions internally use allocator::rebind to test a wide variety of
    // other types.

    std::allocator<int> a;
    StlAllocatorTestbed::TestAlloc( a );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( suzerain_allocator )

BOOST_AUTO_TEST_CASE( allocator_freestore_policy )
{
    // FIXME Test suzerain::allocator_freestore_policy
}

BOOST_AUTO_TEST_CASE( allocator_aligned_malloc_policy )
{
    // FIXME Test suzerain::allocator_aligned_malloc_policy
}

BOOST_AUTO_TEST_SUITE_END()
