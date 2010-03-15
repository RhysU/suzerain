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

#ifdef __ICC
/**
 * Replacement of unknown atomic function by Johannes Singler.
 * Workaround copyright (C) 2007 by Johnannes Singler <singler@ira.uka.de>.
 * @see Intel Software Network Forum Thread <a
 * href="http://software.intel.com/en-us/forums/showthread.php?t=55544"> "Icpc
 * can't compile nothing"</a> for more information.
 **/
#define __sync_fetch_and_add(ptr,addend) _InterlockedExchangeAdd( \
        const_cast<void*>(reinterpret_cast<volatile void*>(ptr)), addend);
#endif

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/allocator.hpp>
#include <suzerain/blas_et_al.hpp>
#include <boost/test/included/unit_test.hpp>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

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

BOOST_AUTO_TEST_CASE( allocator_default_policy )
{
    // Test the default allocator policy
    suzerain::allocator<int> a;
    StlAllocatorTestbed::TestAlloc( a );
}

BOOST_AUTO_TEST_CASE( allocator_freestore_policy )
{
    suzerain::allocator<int, suzerain::allocator_freestore_policy<int> > a;
    StlAllocatorTestbed::TestAlloc( a );
}

BOOST_AUTO_TEST_CASE( allocator_blas )
{
    suzerain::blas::allocator<int>::type a; // Aligned allocator, usually
    StlAllocatorTestbed::TestAlloc( a );
}

BOOST_AUTO_TEST_SUITE_END()
