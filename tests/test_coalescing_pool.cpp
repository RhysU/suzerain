//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#include <suzerain/coalescing_pool.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <suzerain/common.hpp>

// Shorthand
using suzerain::coalescing_pool;
using boost::make_tuple;
using boost::tuple;
using std::size_t;

// Explicit instantiation to flush out any compilation errors
template class coalescing_pool<double>;

// Helper test predicates
#define CHECK_STATS(p, eu, ef, ec)                                       \
    do {                                                                 \
        size_t au, af, ac;                                               \
        size_t at = p.nblocks(au, af, ac);                               \
        tuple<size_t, size_t, size_t> actual   = make_tuple(au, af, ac); \
        tuple<size_t, size_t, size_t> expected = make_tuple(eu, ef, ec); \
        BOOST_CHECK_EQUAL(at, p.nblocks());  /* consistency */           \
        BOOST_CHECK_EQUAL(at, au + af);      /* consistency */           \
        BOOST_CHECK_EQUAL(actual, expected);                             \
    } while(0)

#define REQUIRE_STATS(p, eu, ef, ec)                                     \
    do {                                                                 \
        size_t au, af, ac;                                               \
        size_t at = p.nblocks(au, af, ac);                               \
        tuple<size_t, size_t, size_t> actual   = make_tuple(au, af, ac); \
        tuple<size_t, size_t, size_t> expected = make_tuple(eu, ef, ec); \
        BOOST_REQUIRE_EQUAL(at, p.nblocks());  /* consistency */         \
        BOOST_REQUIRE_EQUAL(at, au + af);      /* consistency */         \
        BOOST_REQUIRE_EQUAL(actual, expected);                           \
    } while(0)


// BOOST_AUTO_TEST_SUITE( suite )

BOOST_AUTO_TEST_CASE( instantiate )
{
    // Basic construction
    coalescing_pool<double> p(10);
    BOOST_CHECK(!p.any_blocks_used());
    BOOST_CHECK_EQUAL(p.blocksize(), 10U);
    BOOST_CHECK_EQUAL(p.narenas(),    0U);
    CHECK_STATS(p, 0, 0, 0);

    // Construction when blocksize is smaller than the intrusive hooks
    coalescing_pool<double> q(1);
    BOOST_CHECK_GE(q.blocksize(), 1U);
    BOOST_CHECK(!q.any_blocks_used());
}

BOOST_AUTO_TEST_CASE( simple_single_arena_reuse )
{
    coalescing_pool<double> p(10);

    coalescing_pool<double>::blocks b1;
    if (!b1) { BOOST_CHECK(true); } else { BOOST_CHECK(false); }
    b1 = p.acquire(3);
    if (b1)  { BOOST_CHECK(true); } else { BOOST_CHECK(false); }

    BOOST_CHECK(p.any_blocks_used());
    BOOST_CHECK_EQUAL(p.narenas(), 1U);
    CHECK_STATS(p, 3, 0, 0);
    BOOST_CHECK_EQUAL(b1.nelems(), 30U);

    p.release(b1);
    BOOST_CHECK(!p.any_blocks_used());
    BOOST_CHECK_EQUAL(p.narenas(), 1U);
    CHECK_STATS(p, 0, 3, 3);

    coalescing_pool<double>::blocks b2 = p.acquire(3);
    BOOST_CHECK(p.any_blocks_used());
    BOOST_CHECK_EQUAL(p.narenas(), 1U);
    CHECK_STATS(p, 3, 0, 0);
    BOOST_CHECK_EQUAL(b2.nelems(), 30U);

    p.release(b2);
    BOOST_CHECK(!p.any_blocks_used());
    BOOST_CHECK_EQUAL(p.narenas(), 1U);
    CHECK_STATS(p, 0, 3, 3);
}

BOOST_AUTO_TEST_CASE( out_of_sequence_single_arena_reuse )
{
    coalescing_pool<double> p(10, /* initial_nblock = */ 5);
    coalescing_pool<double>::blocks c1 = p.acquire(1);
    coalescing_pool<double>::blocks c2 = p.acquire(1);
    coalescing_pool<double>::blocks c3 = p.acquire(1);
    coalescing_pool<double>::blocks c4 = p.acquire(1);
    coalescing_pool<double>::blocks c5 = p.acquire(1);

    BOOST_REQUIRE(p.any_blocks_used());
    BOOST_REQUIRE_EQUAL(p.narenas(), 1U);
    REQUIRE_STATS(p, 5, 0, 0);

    // Check merging with surrounding blocks
    p.release(c1);     CHECK_STATS(p, 4, 1, 1);
    p.release(c3);     CHECK_STATS(p, 3, 2, 1);
    p.release(c2);     CHECK_STATS(p, 2, 3, 3);
    c1 = p.acquire(3); CHECK_STATS(p, 5, 0, 0);
    BOOST_REQUIRE_EQUAL(c1.end(), c4.begin());
    BOOST_REQUIRE_EQUAL(c4.end(), c5.begin());
    BOOST_REQUIRE_EQUAL(p.narenas(), 1U);

    // Check merging with preceding blocks
    p.release(c1);     CHECK_STATS(p, 2, 3, 3);
    p.release(c4);     CHECK_STATS(p, 1, 4, 4);
    c1 = p.acquire(4); CHECK_STATS(p, 5, 0, 0);
    BOOST_REQUIRE_EQUAL(c1.end(), c5.begin());
    BOOST_REQUIRE_EQUAL(p.narenas(), 1U);

    // Check merging with adjacent blocks
    p.release(c5);     CHECK_STATS(p, 4, 1, 1);
    p.release(c1);     CHECK_STATS(p, 0, 5, 5);
    BOOST_REQUIRE_EQUAL(p.narenas(), 1U);

    // Check request for all blocks succeeds
    c1 = p.acquire(5); CHECK_STATS(p, 5, 0, 0);
    p.release(c1);     CHECK_STATS(p, 0, 5, 5);
    BOOST_REQUIRE_EQUAL(p.narenas(), 1U);
}

BOOST_AUTO_TEST_CASE( changing_blocksize )
{
    coalescing_pool<double> p;
    BOOST_CHECK_EQUAL(p.blocksize(), 0U);

    // Force arena allocation at a particular blocksize
    BOOST_REQUIRE_EQUAL(p.blocksize(20), 20U);
    coalescing_pool<double>::blocks b1 = p.acquire(10);
    coalescing_pool<double>::blocks b1_copy = b1;
    if (b1)  { BOOST_CHECK(true); } else { BOOST_CHECK(false); }
    p.release(b1);
    if (!b1) { BOOST_CHECK(true); } else { BOOST_CHECK(false); }
    BOOST_REQUIRE_EQUAL(p.narenas(), 1U);

    // Check dropping the blocksize in half allows reuse
    p.blocksize(10);
    coalescing_pool<double>::blocks b2 = p.acquire(20);
    BOOST_CHECK_EQUAL(b1_copy, b2);
    p.release(b2);
    BOOST_REQUIRE_EQUAL(p.narenas(), 1U);

    // Check changing the blocksize arbitrarily purges arenas
    p.blocksize(7);
    coalescing_pool<double>::blocks b3 = p.acquire(10);
    BOOST_CHECK_NE(b1_copy, b3);
    p.release(b3);
    BOOST_REQUIRE_EQUAL(p.narenas(), 1U);
}

BOOST_AUTO_TEST_CASE( coalesce_arenas )     // Also tests scoped_blocks
{
    coalescing_pool<double> p(10);
    BOOST_CHECK_EQUAL(p.narenas(), 0U);
    CHECK_STATS(p, 0, 0, 0);

    {
        coalescing_pool<double>::scoped_blocks l1(p, 7);
        BOOST_CHECK_EQUAL(p.narenas(), 1U);
        CHECK_STATS(p, 7, 0, 0);

        {
            coalescing_pool<double>::scoped_blocks b1(p, 1);
        }
        BOOST_CHECK_EQUAL(p.narenas(), 2U);
        CHECK_STATS(p, 7, 1, 1);

        coalescing_pool<double>::scoped_blocks l2(p, 8);
        BOOST_CHECK_EQUAL(p.narenas(), 2U);
        CHECK_STATS(p, 15, 0, 0);

        coalescing_pool<double>::scoped_blocks b2(p, 2);
        BOOST_CHECK_EQUAL(p.narenas(), 3U);
        CHECK_STATS(p, 17, 0, 0);

        if (b2)  { BOOST_CHECK(true); } else { BOOST_CHECK(false); }
        b2.release();
        if (!b2) { BOOST_CHECK(true); } else { BOOST_CHECK(false); }
        BOOST_CHECK_EQUAL(p.narenas(), 3U);
        CHECK_STATS(p, 15, 2, 2);

        b2.reset(5);
        BOOST_CHECK_EQUAL(p.narenas(), 3U);
        CHECK_STATS(p, 20, 0, 0);

        b2.release();
        BOOST_CHECK_EQUAL(p.narenas(), 3U);
        CHECK_STATS(p, 15, 5, 5);

        p.drain();
        BOOST_CHECK_EQUAL(p.narenas(), 2U);
        CHECK_STATS(p, 15, 0, 0);
    }

    BOOST_CHECK_EQUAL(p.narenas(), 2U);
    CHECK_STATS(p, 0, 15, 8);
}

BOOST_AUTO_TEST_CASE( degenerate_acquire )
{
    coalescing_pool<double> p(10);
    p.anticipate(0);
    coalescing_pool<double>::blocks b1 = p.acquire(0);
    p.release(b1);
    coalescing_pool<double>::scoped_blocks b2(p, 0);
}

BOOST_AUTO_TEST_CASE( begin_end_rbegin_rend )
{
    coalescing_pool<double> p(6);
    BOOST_REQUIRE_EQUAL(p.blocksize(), 6U);
    coalescing_pool<double>::scoped_blocks b(p, 1);

    coalescing_pool<double>::iterator f = b.begin();
    *f++ = 0;
    *f++ = 1;
    *f++ = 2;
    *f++ = 3;
    *f++ = 4;
    *f++ = 5;
    BOOST_CHECK_EQUAL(f, b.end());

    coalescing_pool<double>::reverse_iterator r = b.rbegin();
    BOOST_CHECK(*r++ == 5);
    BOOST_CHECK(*r++ == 4);
    BOOST_CHECK(*r++ == 3);
    BOOST_CHECK(*r++ == 2);
    BOOST_CHECK(*r++ == 1);
    BOOST_CHECK(*r++ == 0);
    BOOST_CHECK(r == b.rend());
}

// BOOST_AUTO_TEST_SUITE_END()
