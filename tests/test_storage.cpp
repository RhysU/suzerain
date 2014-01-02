//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#include <suzerain/storage.hpp>

#define BOOST_TEST_MAIN
#include <boost/assign.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

BOOST_AUTO_TEST_SUITE( compute_strides_and_storage )

BOOST_AUTO_TEST_CASE( interleaved )
{
    typedef suzerain::storage::interleaved<3> storage;
    const std::size_t D = storage::dimensionality;

    const suzerain::array<int,D> sizes = {{2, 3, 4}};

    {
        BOOST_TEST_MESSAGE(sizes << " without minstrides");
        BOOST_CHECK_EQUAL(24U, storage::compute_storage(
                    sizes.begin()));
        suzerain::array<int,D> strides;
        BOOST_CHECK_EQUAL(24U, storage::compute_strides(
                    sizes.begin(), strides.begin()));
        BOOST_CHECK_EQUAL(strides, boost::assign::list_of(3)(1)(6));
    }

    {
        const suzerain::array<int,D> minstrides[] = { {{1, 1, 1}}, {{3, 1, 6}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(24U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            suzerain::array<int,D> strides;
            BOOST_CHECK_EQUAL(24U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(3)(1)(6));
        }
    }

    {
        const suzerain::array<int,D> minstrides[] = { {{1, 2, 1}}, {{6, 2, 12}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(48U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            suzerain::array<int,D> strides;
            BOOST_CHECK_EQUAL(48U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(6)(2)(12));
        }
    }

    {
        const suzerain::array<int,D> minstrides[] = { {{4, 1, 1}}, {{4, 1, 8}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(32U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            suzerain::array<int,D> strides;
            BOOST_CHECK_EQUAL(32U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(4)(1)(8));
        }
    }

    {
        const suzerain::array<int,D> minstrides[] = { {{1, 1, 24}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(96U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            suzerain::array<int,D> strides;
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

    const suzerain::array<int,D> sizes = {{ 2, 3, 4 }};

    {
        BOOST_TEST_MESSAGE(sizes << " without minstrides");
        BOOST_CHECK_EQUAL(24U, storage::compute_storage(
                    sizes.begin()));
        suzerain::array<int,D> strides;
        BOOST_CHECK_EQUAL(24U, storage::compute_strides(
                    sizes.begin(), strides.begin()));
        BOOST_CHECK_EQUAL(strides, boost::assign::list_of(12)(1)(3));
    }

    {
        const suzerain::array<int,D> minstrides[] = { {{1, 1, 1}}, {{12, 1, 3}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(24U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            suzerain::array<int,D> strides;
            BOOST_CHECK_EQUAL(24U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(12)(1)(3));
        }
    }

    {
        const suzerain::array<int,D> minstrides[] = { {{1, 2, 1}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(48U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            suzerain::array<int,D> strides;
            BOOST_CHECK_EQUAL(48U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(24)(2)(6));
        }
    }

    {
        const suzerain::array<int,D> minstrides[] = { {{1, 1, 4}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(32U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            suzerain::array<int,D> strides;
            BOOST_CHECK_EQUAL(32U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(16)(1)(4));
        }
    }

    {
        const suzerain::array<int,D> minstrides[] = { {{20, 1, 1}} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(40U, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            suzerain::array<int,D> strides;
            BOOST_CHECK_EQUAL(40U, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(20)(1)(3));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
