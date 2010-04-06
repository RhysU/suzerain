#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/assign.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <suzerain/storage.hpp>

#pragma warning(disable:383)

BOOST_AUTO_TEST_SUITE( compute_strides_and_storage )

BOOST_AUTO_TEST_CASE( interleaved )
{
    typedef suzerain::storage::interleaved<3> storage;
    const std::size_t D = storage::dimensionality;

    const boost::array<int,D> sizes = { 2, 3, 4 };

    {
        BOOST_TEST_MESSAGE(sizes << " without minstrides");
        BOOST_CHECK_EQUAL(24, storage::compute_storage(
                    sizes.begin()));
        boost::array<int,D> strides;
        BOOST_CHECK_EQUAL(24, storage::compute_strides(
                    sizes.begin(), strides.begin()));
        BOOST_CHECK_EQUAL(strides, boost::assign::list_of(1)(2)(6));
    }

    {
        const boost::array<int,D> minstrides[] = { {1, 1, 1}, {1, 2, 6} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(24, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(24, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(1)(2)(6));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {2, 1, 1}, {2, 4, 12} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(48, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(48, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(2)(4)(12));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {1, 3, 1}, {1, 3, 8} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(36, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(36, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(1)(3)(9));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {1, 1, 24} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(96, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(96, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(1)(2)(24));
        }
    }
}

BOOST_AUTO_TEST_CASE( noninterleaved )
{
    typedef suzerain::storage::noninterleaved<3> storage;
    const std::size_t D = storage::dimensionality;

    const boost::array<int,D> sizes = { 2, 3, 4 };

    {
        BOOST_TEST_MESSAGE(sizes << " without minstrides");
        BOOST_CHECK_EQUAL(24, storage::compute_storage(
                    sizes.begin()));
        boost::array<int,D> strides;
        BOOST_CHECK_EQUAL(24, storage::compute_strides(
                    sizes.begin(), strides.begin()));
        BOOST_CHECK_EQUAL(strides, boost::assign::list_of(12)(1)(3));
    }

    {
        const boost::array<int,D> minstrides[] = { {1, 1, 1}, {12, 1, 3} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(24, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(24, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(12)(1)(3));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {1, 2, 1} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(48, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(48, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(24)(2)(6));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {1, 1, 4} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(32, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(32, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(16)(1)(4));
        }
    }

    {
        const boost::array<int,D> minstrides[] = { {20, 1, 1} };
        const std::size_t N = sizeof(minstrides)/sizeof(minstrides[0]);

        for (std::size_t i = 0; i < N; ++i) {
            BOOST_TEST_MESSAGE(sizes << " with minstrides " << minstrides[i]);
            BOOST_CHECK_EQUAL(40, storage::compute_storage(
                    sizes.begin(), minstrides[i].begin()));
            boost::array<int,D> strides;
            BOOST_CHECK_EQUAL(40, storage::compute_strides(
                    sizes.begin(), minstrides[i].begin(), strides.begin()));
            BOOST_CHECK_EQUAL(strides, boost::assign::list_of(20)(1)(3));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
