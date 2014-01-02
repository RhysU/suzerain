//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#include <suzerain/running_statistics.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

// Explicit instantiation to flush out compilation errors
template class suzerain::running_statistics<float,  1>;
template class suzerain::running_statistics<double, 5>;

// Types to be tested
typedef boost::mpl::list<float,double> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( examine, T, test_types )
{
    // Test data: Samples, one per row
    static const T data[5][3] = {
        {  2,        -3,         5},
        { -3,         5,        -7},
        {  5,        -7,        11},
        { -7,        11,       -13},
        { 11,       -13,        17},
    };

    static const std::size_t M = sizeof(data   )/sizeof(data[0]   );
    static const std::size_t N = sizeof(data[0])/sizeof(data[0][0]);

    // Test data: Running means after each sample row
    static const T avg[M][N] = {
        {     T(2),        -T(3),         T(5)},
        {  -1/T(2),         T(1),        -T(1)},
        {   4/T(3),      -5/T(3),         T(3)},
        {  -3/T(4),       3/T(2),        -T(1)},
        {   8/T(5),      -7/T(5),      13/T(5)}
    };

    // Test data: Running variance after each sample row
    static const T var[M][N] = {
        {      T(0),        T( 0),       T(  0) },
        {   25/T(2),        T(32),       T( 72) },
        {   49/T(3),    112/T( 3),       T( 84) },
        {  113/T(4),        T(65),       T(120) },
        {  244/T(5),    454/T( 5),   774/T(  5) }
    };

    // Test data: Running minimum after each sample row
    static const T min[M][N] = {
        {  2,        -3,         5},
        { -3,        -3,        -7},
        { -3,        -7,        -7},
        { -7,        -7,       -13},
        { -7,       -13,       -13}
    };

    // Test data: Running maximum after each sample row
    static const T max[M][N] = {
         { 2,        -3,         5},
         { 2,         5,         5},
         { 5,         5,        11},
         { 5,        11,        11},
         {11,        11,        17}
    };

    // Construct the accumulator
    suzerain::running_statistics<T,N> r1;
    BOOST_CHECK_EQUAL(N, r1.size());
    BOOST_TEST_PASSPOINT();

    BOOST_TEST_MESSAGE("Testing initial behavior before any samples provided");
    BOOST_CHECK_EQUAL(0U, r1.count());
    for (std::size_t i = 0; i < N; ++i) {
        BOOST_CHECK((boost::math::isnan)(r1.min(i)));
        BOOST_CHECK((boost::math::isnan)(r1.max(i)));
        BOOST_CHECK((boost::math::isnan)(r1.avg(i)));
        BOOST_CHECK((boost::math::isnan)(r1.var(i)));
        BOOST_CHECK((boost::math::isnan)(r1.std(i)));
    }

    for (std::size_t j = 0; j < M; ++j) {
        BOOST_TEST_MESSAGE("Testing behavior after sample " << (j + 1));
        const T close_enough = 25*(M + 1)*std::numeric_limits<T>::epsilon();
        r1(data[j]);
        for (std::size_t i = 0; i < N; ++i) {
            using std::sqrt;
            BOOST_CHECK_CLOSE(min[j][i], r1.min(i),       close_enough);
            BOOST_CHECK_CLOSE(max[j][i], r1.max(i),       close_enough);
            BOOST_CHECK_CLOSE(avg[j][i], r1.avg(i),       close_enough);
            BOOST_CHECK_CLOSE(var[j][i], r1.var(i),       close_enough);
            BOOST_CHECK_CLOSE(r1.std(i), sqrt(var[j][i]), close_enough);
        }
    }

    // Clear the accumulator and ensure the same behavior persists
    r1.clear();
    BOOST_CHECK_EQUAL(N, r1.size());
    BOOST_TEST_PASSPOINT();

    BOOST_TEST_MESSAGE("Testing post-reset behavior before any samples");
    BOOST_CHECK_EQUAL(0U, r1.count());
    for (std::size_t i = 0; i < N; ++i) {
        BOOST_CHECK((boost::math::isnan)(r1.min(i)));
        BOOST_CHECK((boost::math::isnan)(r1.max(i)));
        BOOST_CHECK((boost::math::isnan)(r1.avg(i)));
        BOOST_CHECK((boost::math::isnan)(r1.var(i)));
        BOOST_CHECK((boost::math::isnan)(r1.std(i)));
    }

    for (std::size_t j = 0; j < M/2; ++j) { // NB Only half used!
        BOOST_TEST_MESSAGE("Testing post-reset behavior after " << (j + 1));
        const T close_enough = 25*(M + 1)*std::numeric_limits<T>::epsilon();
        r1(data[j]);
        for (std::size_t i = 0; i < N; ++i) {
            using std::sqrt;
            BOOST_CHECK_CLOSE(min[j][i], r1.min(i),       close_enough);
            BOOST_CHECK_CLOSE(max[j][i], r1.max(i),       close_enough);
            BOOST_CHECK_CLOSE(avg[j][i], r1.avg(i),       close_enough);
            BOOST_CHECK_CLOSE(var[j][i], r1.var(i),       close_enough);
            BOOST_CHECK_CLOSE(r1.std(i), sqrt(var[j][i]), close_enough);
        }
    }

    // Copy accumulator and ensure we can resume processing
    suzerain::running_statistics<T,N> r2(r1);
    BOOST_CHECK_EQUAL(N, r2.size());
    BOOST_TEST_PASSPOINT();

    for (std::size_t j = M/2; j < M; ++j) { // NB Resuming other half!
        BOOST_TEST_MESSAGE("Testing post-copy behavior after " << (j + 1));
        const T close_enough = 25*(M + 1)*std::numeric_limits<T>::epsilon();
        r2(data[j]);
        for (std::size_t i = 0; i < N; ++i) {
            using std::sqrt;
            BOOST_CHECK_CLOSE(min[j][i], r2.min(i),       close_enough);
            BOOST_CHECK_CLOSE(max[j][i], r2.max(i),       close_enough);
            BOOST_CHECK_CLOSE(avg[j][i], r2.avg(i),       close_enough);
            BOOST_CHECK_CLOSE(var[j][i], r2.var(i),       close_enough);
            BOOST_CHECK_CLOSE(r2.std(i), sqrt(var[j][i]), close_enough);
        }
    }

    // Merge two trivial accumulators and get the results for no samples
    BOOST_TEST_MESSAGE("Testing merge behavior on trivial instances");
    r1.clear();
    r2.clear();
    r1(r2);
    for (std::size_t i = 0; i < N; ++i) {
        BOOST_CHECK((boost::math::isnan)(r1.min(i)));
        BOOST_CHECK((boost::math::isnan)(r1.max(i)));
        BOOST_CHECK((boost::math::isnan)(r1.avg(i)));
        BOOST_CHECK((boost::math::isnan)(r1.var(i)));
        BOOST_CHECK((boost::math::isnan)(r1.std(i)));
    }

    // Merge two accumulators and get the results for expected data
    for (std::size_t k = 0; k < M; ++k) {  // Split on [0, k), [k, M)
        r1.clear();
        r2.clear();
        for (std::size_t j = 0; j < k; ++j) r1(data[j]);  // First  portion
        BOOST_REQUIRE_EQUAL(k, r1.count());
        for (std::size_t j = k; j < M; ++j) r2(data[j]);  // Second portion
        BOOST_REQUIRE_EQUAL(M - k, r2.count());
        r1(r2);                                           // Merge the two

        BOOST_TEST_MESSAGE("Testing merge behavior at split " << k);
        const T close_enough = 25*(M + 1)*std::numeric_limits<T>::epsilon();
        for (std::size_t i = 0; i < N; ++i) {
            BOOST_CHECK_CLOSE(min[M-1][i], r1.min(i),         close_enough);
            BOOST_CHECK_CLOSE(max[M-1][i], r1.max(i),         close_enough);
            BOOST_CHECK_CLOSE(avg[M-1][i], r1.avg(i),         close_enough);
            BOOST_CHECK_CLOSE(var[M-1][i], r1.var(i),         close_enough);
            BOOST_CHECK_CLOSE(r1.std(i),   sqrt(var[M-1][i]), close_enough);
        }
    }
}
