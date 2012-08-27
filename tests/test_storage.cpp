//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/storage.hpp>
#define BOOST_TEST_MAIN
#include <boost/assign.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( compute_strides_and_storage )

BOOST_AUTO_TEST_CASE( interleaved )
{
    typedef suzerain::storage::interleaved<3> storage;
    const std::size_t D = storage::dimensionality;

    const boost::array<int,D> sizes = {{2, 3, 4}};

    {
        BOOST_TEST_MESSAGE(sizes << " without minstrides");
        BOOST_CHECK_EQUAL(24U, storage::compute_storage(
                    sizes.begin()));
        boost::array<int,D> strides;
        BOOST_CHECK_EQUAL(24U, storage::compute_strides(
                    sizes.begin(), strides.begin()));
        BOOST_CHECK_EQUAL(strides, boost::assign::list_of(3)(1)(6));
    }

    {
        const boost::array<int,D> minstrides[] = { {{1, 1, 1}}, {{3, 1, 6}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(24U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(24U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(3)(1)(6));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {{1, 2, 1}}, {{6, 2, 12}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(48U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(48U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(6)(2)(12));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {{4, 1, 1}}, {{4, 1, 8}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(32U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(32U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(4)(1)(8));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {{1, 1, 24}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(96U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(96U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(3)(1)(24));
        }
    }
}

BOOST_AUTO_TEST_CASE( contiguous )
{
    typedef suzerain::storage::contiguous<3> storage;
    const std::size_t D = storage::dimensionality;

    const boost::array<int,D> sizes = {{ 2, 3, 4 }};

    {
        BOOST_TEST_MESSAGE(sizes << " without minstrides");
        BOOST_CHECK_EQUAL(24U, storage::compute_storage(
                    sizes.begin()));
        boost::array<int,D> strides;
        BOOST_CHECK_EQUAL(24U, storage::compute_strides(
                    sizes.begin(), strides.begin()));
        BOOST_CHECK_EQUAL(strides, boost::assign::list_of(12)(1)(3));
    }

    {
        const boost::array<int,D> minstrides[] = { {{1, 1, 1}}, {{12, 1, 3}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(24U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(24U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(12)(1)(3));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {{1, 2, 1}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(48U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(48U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(24)(2)(6));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {{1, 1, 4}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(32U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(32U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(16)(1)(4));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {{20, 1, 1}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(40U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(40U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(20)(1)(3));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
