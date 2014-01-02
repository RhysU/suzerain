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

#include <suzerain/state.hpp>

#define BOOST_TEST_MAIN
#include <boost/concept_check.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

#include "test_tools.hpp"

// Explicit template instantiation to flush out compilation errors
template class suzerain::contiguous_state<4,double>;
template class suzerain::interleaved_state<4,double>;

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

/** Helper for specifying 3D extents information */
static suzerain::array<int,3> size(int x, int y, int z)
{
    suzerain::array<int,3> a = {{ x, y, z }};
    return a;
}

/** Helper for specifying 4D extents information */
static suzerain::array<int,4> size(int w, int x, int y, int z)
{
    suzerain::array<int,4> a = {{ w, x, y, z }};
    return a;
}

/** Helper for specifying 3D extents information for load and verify */
static suzerain::array<int,3> size223() {
    return size(2,2,3);
}

/** Helper for specifying 4D extents information for load and verify */
static suzerain::array<int,4> size2234() {
    return size(2,2,3,4);
}

/** Load a 3D instance with real test data */
template<
    template <std::size_t,typename> class State,
    typename FPT,
    typename Scale
>
static void load(
    State<3,FPT> &state,
    const Scale scale,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3U);

    state[0][0][0] = FPT(  2) * scale;
    state[0][1][0] = FPT(  3) * scale;
    state[0][0][1] = FPT(  5) * scale;
    state[0][1][1] = FPT(  7) * scale;
    state[0][0][2] = FPT( 11) * scale;
    state[0][1][2] = FPT( 13) * scale;
    state[1][0][0] = FPT(102) * scale;
    state[1][1][0] = FPT(103) * scale;
    state[1][0][1] = FPT(105) * scale;
    state[1][1][1] = FPT(107) * scale;
    state[1][0][2] = FPT(111) * scale;
    state[1][1][2] = FPT(113) * scale;
}

/** Load a 4D instance with real test data */
template<
    template <std::size_t,typename> class State,
    typename FPT,
    typename Scale
>
static void load(
    State<4,FPT> &state,
    const Scale scale,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3U);
    BOOST_REQUIRE_EQUAL(state.shape()[3], 4U);

    state[0][0][0][0] = FPT(  2) * scale;
    state[0][1][0][0] = FPT(  3) * scale;
    state[0][0][1][0] = FPT(  5) * scale;
    state[0][1][1][0] = FPT(  7) * scale;
    state[0][0][2][0] = FPT( 11) * scale;
    state[0][1][2][0] = FPT( 13) * scale;

    state[0][0][0][1] = FPT( 17) * scale;
    state[0][1][0][1] = FPT( 19) * scale;
    state[0][0][1][1] = FPT( 23) * scale;
    state[0][1][1][1] = FPT( 29) * scale;
    state[0][0][2][1] = FPT( 31) * scale;
    state[0][1][2][1] = FPT( 37) * scale;

    state[0][0][0][2] = FPT( 41) * scale;
    state[0][1][0][2] = FPT( 43) * scale;
    state[0][0][1][2] = FPT( 47) * scale;
    state[0][1][1][2] = FPT( 53) * scale;
    state[0][0][2][2] = FPT( 59) * scale;
    state[0][1][2][2] = FPT( 61) * scale;

    state[0][0][0][3] = FPT( 67) * scale;
    state[0][1][0][3] = FPT( 71) * scale;
    state[0][0][1][3] = FPT( 73) * scale;
    state[0][1][1][3] = FPT( 79) * scale;
    state[0][0][2][3] = FPT( 83) * scale;
    state[0][1][2][3] = FPT( 89) * scale;

    state[1][0][0][0] = FPT(102) * scale;
    state[1][1][0][0] = FPT(103) * scale;
    state[1][0][1][0] = FPT(105) * scale;
    state[1][1][1][0] = FPT(107) * scale;
    state[1][0][2][0] = FPT(111) * scale;
    state[1][1][2][0] = FPT(113) * scale;

    state[1][0][0][1] = FPT(117) * scale;
    state[1][1][0][1] = FPT(119) * scale;
    state[1][0][1][1] = FPT(123) * scale;
    state[1][1][1][1] = FPT(129) * scale;
    state[1][0][2][1] = FPT(131) * scale;
    state[1][1][2][1] = FPT(137) * scale;

    state[1][0][0][2] = FPT(141) * scale;
    state[1][1][0][2] = FPT(143) * scale;
    state[1][0][1][2] = FPT(147) * scale;
    state[1][1][1][2] = FPT(153) * scale;
    state[1][0][2][2] = FPT(159) * scale;
    state[1][1][2][2] = FPT(161) * scale;

    state[1][0][0][3] = FPT(167) * scale;
    state[1][1][0][3] = FPT(171) * scale;
    state[1][0][1][3] = FPT(173) * scale;
    state[1][1][1][3] = FPT(179) * scale;
    state[1][0][2][3] = FPT(183) * scale;
    state[1][1][2][3] = FPT(189) * scale;
}

/** Verify a 3D instance against real test data */
template<
    template <std::size_t,typename> class State,
    typename FPT,
    typename Scale
>
static void verify(
    const State<3,FPT> &state,
    const Scale scale,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3U);

    BOOST_CHECK_EQUAL(state[0][0][0], FPT(  2) * scale);
    BOOST_CHECK_EQUAL(state[0][1][0], FPT(  3) * scale);
    BOOST_CHECK_EQUAL(state[0][0][1], FPT(  5) * scale);
    BOOST_CHECK_EQUAL(state[0][1][1], FPT(  7) * scale);
    BOOST_CHECK_EQUAL(state[0][0][2], FPT( 11) * scale);
    BOOST_CHECK_EQUAL(state[0][1][2], FPT( 13) * scale);
    BOOST_CHECK_EQUAL(state[1][0][0], FPT(102) * scale);
    BOOST_CHECK_EQUAL(state[1][1][0], FPT(103) * scale);
    BOOST_CHECK_EQUAL(state[1][0][1], FPT(105) * scale);
    BOOST_CHECK_EQUAL(state[1][1][1], FPT(107) * scale);
    BOOST_CHECK_EQUAL(state[1][0][2], FPT(111) * scale);
    BOOST_CHECK_EQUAL(state[1][1][2], FPT(113) * scale);
}

/** Verify a 4D instance against real test data */
template<
    template <std::size_t,typename> class State,
    typename FPT,
    typename Scale
>
static void verify(
    const State<4,FPT> &state,
    const Scale scale,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3U);
    BOOST_REQUIRE_EQUAL(state.shape()[3], 4U);

    BOOST_CHECK_EQUAL(state[0][0][0][0], FPT(  2) * scale);
    BOOST_CHECK_EQUAL(state[0][1][0][0], FPT(  3) * scale);
    BOOST_CHECK_EQUAL(state[0][0][1][0], FPT(  5) * scale);
    BOOST_CHECK_EQUAL(state[0][1][1][0], FPT(  7) * scale);
    BOOST_CHECK_EQUAL(state[0][0][2][0], FPT( 11) * scale);
    BOOST_CHECK_EQUAL(state[0][1][2][0], FPT( 13) * scale);

    BOOST_CHECK_EQUAL(state[0][0][0][1], FPT( 17) * scale);
    BOOST_CHECK_EQUAL(state[0][1][0][1], FPT( 19) * scale);
    BOOST_CHECK_EQUAL(state[0][0][1][1], FPT( 23) * scale);
    BOOST_CHECK_EQUAL(state[0][1][1][1], FPT( 29) * scale);
    BOOST_CHECK_EQUAL(state[0][0][2][1], FPT( 31) * scale);
    BOOST_CHECK_EQUAL(state[0][1][2][1], FPT( 37) * scale);

    BOOST_CHECK_EQUAL(state[0][0][0][2], FPT( 41) * scale);
    BOOST_CHECK_EQUAL(state[0][1][0][2], FPT( 43) * scale);
    BOOST_CHECK_EQUAL(state[0][0][1][2], FPT( 47) * scale);
    BOOST_CHECK_EQUAL(state[0][1][1][2], FPT( 53) * scale);
    BOOST_CHECK_EQUAL(state[0][0][2][2], FPT( 59) * scale);
    BOOST_CHECK_EQUAL(state[0][1][2][2], FPT( 61) * scale);

    BOOST_CHECK_EQUAL(state[0][0][0][3], FPT( 67) * scale);
    BOOST_CHECK_EQUAL(state[0][1][0][3], FPT( 71) * scale);
    BOOST_CHECK_EQUAL(state[0][0][1][3], FPT( 73) * scale);
    BOOST_CHECK_EQUAL(state[0][1][1][3], FPT( 79) * scale);
    BOOST_CHECK_EQUAL(state[0][0][2][3], FPT( 83) * scale);
    BOOST_CHECK_EQUAL(state[0][1][2][3], FPT( 89) * scale);

    BOOST_CHECK_EQUAL(state[1][0][0][0], FPT(102) * scale);
    BOOST_CHECK_EQUAL(state[1][1][0][0], FPT(103) * scale);
    BOOST_CHECK_EQUAL(state[1][0][1][0], FPT(105) * scale);
    BOOST_CHECK_EQUAL(state[1][1][1][0], FPT(107) * scale);
    BOOST_CHECK_EQUAL(state[1][0][2][0], FPT(111) * scale);
    BOOST_CHECK_EQUAL(state[1][1][2][0], FPT(113) * scale);

    BOOST_CHECK_EQUAL(state[1][0][0][1], FPT(117) * scale);
    BOOST_CHECK_EQUAL(state[1][1][0][1], FPT(119) * scale);
    BOOST_CHECK_EQUAL(state[1][0][1][1], FPT(123) * scale);
    BOOST_CHECK_EQUAL(state[1][1][1][1], FPT(129) * scale);
    BOOST_CHECK_EQUAL(state[1][0][2][1], FPT(131) * scale);
    BOOST_CHECK_EQUAL(state[1][1][2][1], FPT(137) * scale);

    BOOST_CHECK_EQUAL(state[1][0][0][2], FPT(141) * scale);
    BOOST_CHECK_EQUAL(state[1][1][0][2], FPT(143) * scale);
    BOOST_CHECK_EQUAL(state[1][0][1][2], FPT(147) * scale);
    BOOST_CHECK_EQUAL(state[1][1][1][2], FPT(153) * scale);
    BOOST_CHECK_EQUAL(state[1][0][2][2], FPT(159) * scale);
    BOOST_CHECK_EQUAL(state[1][1][2][2], FPT(161) * scale);

    BOOST_CHECK_EQUAL(state[1][0][0][3], FPT(167) * scale);
    BOOST_CHECK_EQUAL(state[1][1][0][3], FPT(171) * scale);
    BOOST_CHECK_EQUAL(state[1][0][1][3], FPT(173) * scale);
    BOOST_CHECK_EQUAL(state[1][1][1][3], FPT(179) * scale);
    BOOST_CHECK_EQUAL(state[1][0][2][3], FPT(183) * scale);
    BOOST_CHECK_EQUAL(state[1][1][2][3], FPT(189) * scale);
}

/** Load a 3D instance with complex test data */
template<
    template <std::size_t,typename> class State,
    typename FPT,
    typename Scale
>
static void load(
    State<3,std::complex<FPT> > &state,
    const Scale scale,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3U);

    typedef typename std::complex<FPT> complex;
    const FPT scaleFactor = boost::numeric_cast<FPT>(scale);
    state[0][0][0] = complex(  2, -  2) * scaleFactor;
    state[0][1][0] = complex(  3, -  3) * scaleFactor;
    state[0][0][1] = complex(  5, -  5) * scaleFactor;
    state[0][1][1] = complex(  7, -  7) * scaleFactor;
    state[0][0][2] = complex( 11, - 11) * scaleFactor;
    state[0][1][2] = complex( 13, - 13) * scaleFactor;
    state[1][0][0] = complex(102, -102) * scaleFactor;
    state[1][1][0] = complex(103, -103) * scaleFactor;
    state[1][0][1] = complex(105, -105) * scaleFactor;
    state[1][1][1] = complex(107, -107) * scaleFactor;
    state[1][0][2] = complex(111, -111) * scaleFactor;
    state[1][1][2] = complex(113, -113) * scaleFactor;
}

/** Load a 4D instance with complex test data */
template<
    template <std::size_t,typename> class State,
    typename FPT,
    typename Scale
>
static void load(
    State<4,std::complex<FPT> > &state,
    const Scale scale,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3U);
    BOOST_REQUIRE_EQUAL(state.shape()[3], 4U);

    typedef typename std::complex<FPT> complex;
    const FPT scaleFactor = boost::numeric_cast<FPT>(scale);
    state[0][0][0][0] = complex(  2, -  2) * scaleFactor;
    state[0][1][0][0] = complex(  3, -  3) * scaleFactor;
    state[0][0][1][0] = complex(  5, -  5) * scaleFactor;
    state[0][1][1][0] = complex(  7, -  7) * scaleFactor;
    state[0][0][2][0] = complex( 11, - 11) * scaleFactor;
    state[0][1][2][0] = complex( 13, - 13) * scaleFactor;

    state[0][0][0][1] = complex( 17, - 17) * scaleFactor;
    state[0][1][0][1] = complex( 19, - 19) * scaleFactor;
    state[0][0][1][1] = complex( 23, - 23) * scaleFactor;
    state[0][1][1][1] = complex( 29, - 29) * scaleFactor;
    state[0][0][2][1] = complex( 31, - 31) * scaleFactor;
    state[0][1][2][1] = complex( 37, - 37) * scaleFactor;

    state[0][0][0][2] = complex( 41, - 41) * scaleFactor;
    state[0][1][0][2] = complex( 43, - 43) * scaleFactor;
    state[0][0][1][2] = complex( 47, - 47) * scaleFactor;
    state[0][1][1][2] = complex( 53, - 53) * scaleFactor;
    state[0][0][2][2] = complex( 59, - 59) * scaleFactor;
    state[0][1][2][2] = complex( 61, - 61) * scaleFactor;

    state[0][0][0][3] = complex( 67, - 67) * scaleFactor;
    state[0][1][0][3] = complex( 71, - 71) * scaleFactor;
    state[0][0][1][3] = complex( 73, - 73) * scaleFactor;
    state[0][1][1][3] = complex( 79, - 79) * scaleFactor;
    state[0][0][2][3] = complex( 83, - 83) * scaleFactor;
    state[0][1][2][3] = complex( 89, - 89) * scaleFactor;

    state[1][0][0][0] = complex(102, -102) * scaleFactor;
    state[1][1][0][0] = complex(103, -103) * scaleFactor;
    state[1][0][1][0] = complex(105, -105) * scaleFactor;
    state[1][1][1][0] = complex(107, -107) * scaleFactor;
    state[1][0][2][0] = complex(111, -111) * scaleFactor;
    state[1][1][2][0] = complex(113, -113) * scaleFactor;

    state[1][0][0][1] = complex(117, -117) * scaleFactor;
    state[1][1][0][1] = complex(119, -119) * scaleFactor;
    state[1][0][1][1] = complex(123, -123) * scaleFactor;
    state[1][1][1][1] = complex(129, -129) * scaleFactor;
    state[1][0][2][1] = complex(131, -131) * scaleFactor;
    state[1][1][2][1] = complex(137, -137) * scaleFactor;

    state[1][0][0][2] = complex(141, -141) * scaleFactor;
    state[1][1][0][2] = complex(143, -143) * scaleFactor;
    state[1][0][1][2] = complex(147, -147) * scaleFactor;
    state[1][1][1][2] = complex(153, -153) * scaleFactor;
    state[1][0][2][2] = complex(159, -159) * scaleFactor;
    state[1][1][2][2] = complex(161, -161) * scaleFactor;

    state[1][0][0][3] = complex(167, -167) * scaleFactor;
    state[1][1][0][3] = complex(171, -171) * scaleFactor;
    state[1][0][1][3] = complex(173, -173) * scaleFactor;
    state[1][1][1][3] = complex(179, -179) * scaleFactor;
    state[1][0][2][3] = complex(183, -183) * scaleFactor;
    state[1][1][2][3] = complex(189, -189) * scaleFactor;
}

/** Verify a 3D instance against complex test data */
template<
    template <std::size_t,typename> class State,
    typename FPT,
    typename Scale
>
static void verify(
    const State<3,std::complex<FPT> > &state,
    const Scale scale,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3U);

    typedef typename std::complex<FPT> complex;
    const FPT scaleFactor = boost::numeric_cast<FPT>(scale);
    BOOST_CHECK_EQUAL(state[0][0][0], complex(  2, -  2) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][0], complex(  3, -  3) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][1], complex(  5, -  5) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][1], complex(  7, -  7) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][2], complex( 11, - 11) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][2], complex( 13, - 13) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][0], complex(102, -102) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][0], complex(103, -103) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][1], complex(105, -105) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][1], complex(107, -107) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][2], complex(111, -111) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][2], complex(113, -113) * scaleFactor);
}

/** Verify a 4D instance against complex test data */
template<
    template <std::size_t,typename> class State,
    typename FPT,
    typename Scale
>
static void verify(
    const State<4,std::complex<FPT> > &state,
    const Scale scale,
    typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2U);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3U);

    typedef typename std::complex<FPT> complex;
    const FPT scaleFactor = boost::numeric_cast<FPT>(scale);
    BOOST_CHECK_EQUAL(state[0][0][0][0], complex(  2, -  2) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][0][0], complex(  3, -  3) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][1][0], complex(  5, -  5) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][1][0], complex(  7, -  7) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][2][0], complex( 11, - 11) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][2][0], complex( 13, - 13) * scaleFactor);

    BOOST_CHECK_EQUAL(state[0][0][0][1], complex( 17, - 17) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][0][1], complex( 19, - 19) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][1][1], complex( 23, - 23) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][1][1], complex( 29, - 29) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][2][1], complex( 31, - 31) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][2][1], complex( 37, - 37) * scaleFactor);

    BOOST_CHECK_EQUAL(state[0][0][0][2], complex( 41, - 41) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][0][2], complex( 43, - 43) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][1][2], complex( 47, - 47) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][1][2], complex( 53, - 53) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][2][2], complex( 59, - 59) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][2][2], complex( 61, - 61) * scaleFactor);

    BOOST_CHECK_EQUAL(state[0][0][0][3], complex( 67, - 67) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][0][3], complex( 71, - 71) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][1][3], complex( 73, - 73) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][1][3], complex( 79, - 79) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][0][2][3], complex( 83, - 83) * scaleFactor);
    BOOST_CHECK_EQUAL(state[0][1][2][3], complex( 89, - 89) * scaleFactor);

    BOOST_CHECK_EQUAL(state[1][0][0][0], complex(102, -102) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][0][0], complex(103, -103) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][1][0], complex(105, -105) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][1][0], complex(107, -107) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][2][0], complex(111, -111) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][2][0], complex(113, -113) * scaleFactor);

    BOOST_CHECK_EQUAL(state[1][0][0][1], complex(117, -117) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][0][1], complex(119, -119) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][1][1], complex(123, -123) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][1][1], complex(129, -129) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][2][1], complex(131, -131) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][2][1], complex(137, -137) * scaleFactor);

    BOOST_CHECK_EQUAL(state[1][0][0][2], complex(141, -141) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][0][2], complex(143, -143) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][1][2], complex(147, -147) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][1][2], complex(153, -153) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][2][2], complex(159, -159) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][2][2], complex(161, -161) * scaleFactor);

    BOOST_CHECK_EQUAL(state[1][0][0][3], complex(167, -167) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][0][3], complex(171, -171) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][1][3], complex(173, -173) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][1][3], complex(179, -179) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][0][2][3], complex(183, -183) * scaleFactor);
    BOOST_CHECK_EQUAL(state[1][1][2][3], complex(189, -189) * scaleFactor);
}

template<
    std::size_t NumDims,
    template <std::size_t,typename> class State1,
    template <std::size_t,typename> class State2,
    typename Element
>
static void test_assignment_helper(
    State1<NumDims,Element> &foo,
    State2<NumDims,Element> &bar)
{
    load(foo, 1);

    foo.assign_from(foo); // Self
    verify(foo, 1);

    bar.assign_from(foo);
    verify(bar, 1);
}

template<
    std::size_t NumDims,
    template <std::size_t,typename> class State,
    typename Element
>
static void test_scale_helper(State<NumDims,Element> &foo)
{
    BOOST_TEST_PASSPOINT();
    load(foo, 1);
    BOOST_TEST_PASSPOINT();
    verify(foo, 1);

    BOOST_TEST_PASSPOINT();
    foo.scale(1);
    BOOST_TEST_PASSPOINT();
    verify(foo, 1);
    BOOST_TEST_PASSPOINT();

    foo.scale(2);
    BOOST_TEST_PASSPOINT();
    verify(foo, 2);
    BOOST_TEST_PASSPOINT();

    foo.scale(0);
    BOOST_TEST_PASSPOINT();
    verify(foo, 0);
    BOOST_TEST_PASSPOINT();
}

template<
    std::size_t NumDims,
    template <std::size_t,typename> class State1,
    template <std::size_t,typename> class State2,
    typename Element
>
static void test_add_scaled_helper(
    State1<NumDims,Element> &foo,
    State2<NumDims,Element> &bar)
{
    load(foo, 1);
    load(bar, 2);
    foo.add_scaled(3, bar);
    verify(foo, 7);
}

template<
    std::size_t NumDims,
    template <std::size_t,typename> class State1,
    template <std::size_t,typename> class State2,
    typename Element
>
static void test_exchange_helper(
    State1<NumDims,Element> &foo,
    State2<NumDims,Element> &bar)
{
    load(foo, 1);
    load(bar, 2);
    foo.exchange(bar);
    verify(foo, 2);
    verify(bar, 1);
    bar.exchange(foo);
    verify(foo, 1);
    verify(bar, 2);
}

BOOST_AUTO_TEST_SUITE( InterleavedState )

using suzerain::interleaved_state;

// Types to be tested with interleaved_state
typedef boost::mpl::list<
    double
    ,float
    ,std::complex<double>
    ,std::complex<float>
> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors, T, test_types )
{
    using boost::numeric_cast;

    // Regular constructor
    interleaved_state<4,T> foo(size2234());
    BOOST_CHECK_EQUAL(foo.range().size(), 2*2*3*4);
    load(foo, 1);
    verify(foo, 1);

    // Copy construct a second instance from the first
    interleaved_state<4,T> bar(foo);
    BOOST_CHECK_EQUAL(bar.range().size(), 2*2*3*4);
    verify(bar, 1);

    // Modify first instance's data
    for (int i = 0; i < numeric_cast<int>(foo.shape()[0]); ++i)
        for (int j = 0; j < numeric_cast<int>(foo.shape()[1]); ++j)
            for (int k = 0; k < numeric_cast<int>(foo.shape()[2]); ++k)
                for (int l = 0; l < numeric_cast<int>(foo.shape()[3]); ++l)
                    foo[i][j][k][l] += (i+1)*(j+1)*(k+1)*(l+1);

    // Ensure copy constructed data in second instance not modified
    verify(bar, 1);

    // Create padded instance and ensure content lies within padding
    interleaved_state<4,T> baz(size2234(), 2*2*3*4*7);
    BOOST_CHECK_EQUAL(baz.range().size(), 2*2*3*4*7);
    BOOST_CHECK_GE(baz.range().begin(), &(baz[0][0][0][0]));
    BOOST_CHECK_LT(&(baz[1][1][2][0]), baz.range().end());

    // Ensure padded information present propagated in copy operations
    interleaved_state<4,T> qux(baz);
    BOOST_CHECK_EQUAL(qux.range().size(), 2*2*3*4*7);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, T, test_types )
{
    BOOST_TEST_MESSAGE("Instances without padding");
    {
        interleaved_state<4,T> foo(size2234()), bar(size2234());
        test_assignment_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("First instance with padding");
    {
        interleaved_state<4,T> foo(size2234(), 87), bar(size2234());
        test_assignment_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Second instance with padding");
    {
        interleaved_state<4,T> foo(size2234()), bar(size2234(), 87);
        test_assignment_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        interleaved_state<4,T> foo(size2234(), 87), bar(size2234(), 78);
        test_assignment_helper(foo, bar);
    }

    // Operation between two nonconforming states throws
    interleaved_state<4,T> foo(size2234());
    interleaved_state<4,T> baz(size(2,2,2,2));
    BOOST_CHECK_THROW(baz.assign_from(foo), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( storage_order, T, test_types )
{
    interleaved_state<4,T> foo(size(2,2,2,2));

    BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +     1, &(foo[0][1][0][0]));
    BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +     2, &(foo[1][0][0][0]));
    BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +   2*2, &(foo[0][0][1][0]));
    BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) + 2*2*2, &(foo[0][0][0][1]));

    BOOST_CHECK_EQUAL(     1, foo.strides()[1] );
    BOOST_CHECK_EQUAL(     2, foo.strides()[0] );
    BOOST_CHECK_EQUAL(   2*2, foo.strides()[2] );
    BOOST_CHECK_EQUAL( 2*2*2, foo.strides()[3] );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( is_isomorphic, T, test_types )
{
    interleaved_state<4,T> foo(  size(2,2,2,2));
    interleaved_state<4,T> bar(  size(2,2,2,2));
    interleaved_state<4,T> baz(  size(1,2,2,2));
    interleaved_state<4,T> qux(  size(2,1,2,2));
    interleaved_state<4,T> quux( size(2,2,1,2));
    interleaved_state<4,T> quuux(size(2,2,2,1));

    BOOST_CHECK_EQUAL(true,  foo.is_isomorphic(foo));
    BOOST_CHECK_EQUAL(true,  foo.is_isomorphic(bar));
    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(baz));
    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(qux));
    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(quux));
    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(quuux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( scale, T, test_types )
{
    BOOST_TEST_MESSAGE("Instance without padding");
    {
        interleaved_state<4,T> foo(size2234());
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding");
    {
        interleaved_state<4,T> foo(size2234(), 87);
        test_scale_helper(foo);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( add_scaled, T, test_types )
{
    BOOST_TEST_MESSAGE("Instances without padding");
    {
        interleaved_state<4,T> foo(size2234()), bar(size2234());
        test_add_scaled_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("First instance with padding");
    {
        interleaved_state<4,T> foo(size2234(), 87), bar(size2234());
        test_add_scaled_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Second instance with padding");
    {
        interleaved_state<4,T> foo(size2234()), bar(size2234(), 87);
        test_add_scaled_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        interleaved_state<4,T> foo(size2234(), 87), bar(size2234(), 78);
        test_add_scaled_helper(foo, bar);
    }

    // Operation between two nonconforming states throws
    interleaved_state<4,T> foo(size2234());
    interleaved_state<4,T> baz(size(2,2,2,2));
    BOOST_CHECK_THROW(foo.add_scaled(3, baz), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange, T, test_types )
{
    BOOST_TEST_MESSAGE("Instances without padding");
    {
        interleaved_state<4,T> foo(size2234()), bar(size2234());
        test_exchange_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("First instance with padding");
    {
        interleaved_state<4,T> foo(size2234(), 87), bar(size2234());
        test_exchange_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Second instance with padding");
    {
        interleaved_state<4,T> foo(size2234()), bar(size2234(), 87);
        test_exchange_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        interleaved_state<4,T> foo(size2234(), 87), bar(size2234(), 78);
        test_exchange_helper(foo, bar);
    }

    // Operation between two nonconforming states throws
    interleaved_state<4,T> foo(size2234());
    interleaved_state<4,T> baz(size(2,2,2,2));
    BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( concept_check, T, test_types )
{
    using boost::detail::multi_array::MutableMultiArrayConcept;
    BOOST_CONCEPT_ASSERT((MutableMultiArrayConcept<interleaved_state<4,T>,4>));
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( ContiguousState )

using suzerain::contiguous_state;

// Types to be tested with contiguous_state
typedef boost::mpl::list<
    double
    ,float
    ,std::complex<double>
    ,std::complex<float>
> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors3, T, test_types )
{
    using boost::numeric_cast;

    // Regular constructor
    contiguous_state<3,T> foo(size223());
    BOOST_CHECK_EQUAL(foo.range().size(), 2*2*3);
    load(foo, 1);
    verify(foo, 1);

    // Copy construct a second instance from the first
    contiguous_state<3,T> bar(foo);
    BOOST_CHECK_EQUAL(bar.range().size(), 2*2*3);
    verify(bar, 1);

    // Modify first instance's data
    for (int i = 0; i < numeric_cast<int>(foo.shape()[0]); ++i)
        for (int j = 0; j < numeric_cast<int>(foo.shape()[1]); ++j)
            for (int k = 0; k < numeric_cast<int>(foo.shape()[2]); ++k)
                foo[i][j][k] += (i+1)*(j+1)*(k+1);

    // Ensure copy constructed data in second instance not modified
    verify(bar, 1);

    // Create padded instance and ensure content lies within padding
    contiguous_state<3,T> baz(size223(), size(7,1,1));
    BOOST_CHECK_EQUAL(baz.range().size(), 2*7);
    BOOST_CHECK_EQUAL(baz.strides()[0], 7);
    BOOST_CHECK_GE(baz.range().begin(), &(baz[0][0][0]));
    BOOST_CHECK_LT(&(baz[1][1][2]), baz.range().end());
    load(baz, 1);
    verify(baz, 1);

    // Ensure padded information present propagated in copy operations
    contiguous_state<3,T> qux(baz);
    BOOST_CHECK_EQUAL(qux.range().size(), 2*7);
    BOOST_CHECK_EQUAL(baz.strides()[0], qux.strides()[0]);
    BOOST_CHECK_EQUAL(baz.strides()[1], qux.strides()[1]);
    BOOST_CHECK_EQUAL(baz.strides()[2], qux.strides()[2]);
    verify(qux, 1);

    // Modify first padded instance's data
    for (int i = 0; i < numeric_cast<int>(foo.shape()[0]); ++i)
        for (int j = 0; j < numeric_cast<int>(foo.shape()[1]); ++j)
            for (int k = 0; k < numeric_cast<int>(foo.shape()[2]); ++k)
                baz[i][j][k] += (i+1)*(j+1)*(k+1);

    // Ensure copy constructed data in second padded instance not modified
    verify(qux, 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors4, T, test_types )
{
    using boost::numeric_cast;

    // Regular constructor
    contiguous_state<4,T> foo(size2234());
    BOOST_CHECK_EQUAL(foo.range().size(), 2*2*3*4);
    load(foo, 1);
    verify(foo, 1);

    // Copy construct a second instance from the first
    contiguous_state<4,T> bar(foo);
    BOOST_CHECK_EQUAL(bar.range().size(), 2*2*3*4);
    verify(bar, 1);

    // Modify first instance's data
    for (int i = 0; i < numeric_cast<int>(foo.shape()[0]); ++i)
        for (int j = 0; j < numeric_cast<int>(foo.shape()[1]); ++j)
            for (int k = 0; k < numeric_cast<int>(foo.shape()[2]); ++k)
                for (int l = 0; l < numeric_cast<int>(foo.shape()[3]); ++l)
                    foo[i][j][k][l] += (i+1)*(j+1)*(k+1)*(l+1);

    // Ensure copy constructed data in second instance not modified
    verify(bar, 1);

    // Create padded instance and ensure content lies within padding
    contiguous_state<4,T> baz(size2234(), size(27,1,1,1));
    BOOST_CHECK_EQUAL(baz.range().size(), 2*27);
    BOOST_CHECK_EQUAL(baz.strides()[0], 27);
    BOOST_CHECK_GE(baz.range().begin(), &(baz[0][0][0][0]));
    BOOST_CHECK_LT(&(baz[1][1][2][3]), baz.range().end());
    load(baz, 1);
    verify(baz, 1);

    // Ensure padded information present propagated in copy operations
    contiguous_state<4,T> qux(baz);
    BOOST_CHECK_EQUAL(qux.range().size(), 2*27);
    BOOST_CHECK_EQUAL(baz.strides()[0], qux.strides()[0]);
    BOOST_CHECK_EQUAL(baz.strides()[1], qux.strides()[1]);
    BOOST_CHECK_EQUAL(baz.strides()[2], qux.strides()[2]);
    BOOST_CHECK_EQUAL(baz.strides()[3], qux.strides()[3]);
    verify(qux, 1);

    // Modify first padded instance's data
    for (int i = 0; i < numeric_cast<int>(foo.shape()[0]); ++i)
        for (int j = 0; j < numeric_cast<int>(foo.shape()[1]); ++j)
            for (int k = 0; k < numeric_cast<int>(foo.shape()[2]); ++k)
                for (int l = 0; l < numeric_cast<int>(foo.shape()[2]); ++l)
                    baz[i][j][k][l] += (i+1)*(j+1)*(k+1)*(l+1);

    // Ensure copy constructed data in second padded instance not modified
    verify(qux, 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment3, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        contiguous_state<3,T> foo(size223()), bar(size223());
        test_assignment_helper(foo, bar);

        // Operation between two nonconforming states throws
        contiguous_state<3,T> baz(size(2,2,2));
        BOOST_CHECK_THROW(baz.assign_from(foo), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<3,T> foo(size223(),size(7,1,1));
        contiguous_state<3,T> bar(size223(),size(10,1,1));
        test_assignment_helper(foo, bar);

        contiguous_state<3,T> baz(size223(),size(1,2,1));
        contiguous_state<3,T> qux(size223(),size(1,1,4));
        test_assignment_helper(baz, qux);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        contiguous_state<3,T> foo(size223()), bar(size223(),size(7,1,1));
        test_assignment_helper(foo, bar);

        contiguous_state<3,T> baz(size223(),size(1,2,1));
        test_assignment_helper(foo, baz);

        contiguous_state<3,T> qux(size223(),size(1,2,3));
        test_assignment_helper(foo, qux);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        contiguous_state<3,T> foo(size223(),size(7,1,1)), bar(size223());
        test_assignment_helper(foo, bar);

        contiguous_state<3,T> baz(size223(),size(1,3,1));
        test_assignment_helper(baz, bar);

        contiguous_state<3,T> qux(size223(),size(4,4,4));
        test_assignment_helper(qux, bar);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment4, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        contiguous_state<4,T> foo(size2234()), bar(size2234());
        test_assignment_helper(foo, bar);

        // Operation between two nonconforming states throws
        contiguous_state<4,T> baz(size(2,2,2,2));
        BOOST_CHECK_THROW(baz.assign_from(foo), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<4,T> foo(size2234(),size(27,1,1,1));
        contiguous_state<4,T> bar(size2234(),size(29,1,1,1));
        test_assignment_helper(foo, bar);

        contiguous_state<4,T> baz(size2234(),size(1,2,1,1));
        contiguous_state<4,T> qux(size2234(),size(1,1,4,1));
        test_assignment_helper(baz, qux);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        contiguous_state<4,T> foo(size2234());
        contiguous_state<4,T> bar(size2234(),size(33,1,1,1));
        test_assignment_helper(foo, bar);

        contiguous_state<4,T> baz(size2234(),size(1,2,1,1));
        test_assignment_helper(foo, baz);

        contiguous_state<4,T> qux(size2234(),size(1,2,3,1));
        test_assignment_helper(foo, qux);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        contiguous_state<4,T> foo(size2234(),size(72,1,1,8));
        contiguous_state<4,T> bar(size2234());
        test_assignment_helper(foo, bar);

        contiguous_state<4,T> baz(size2234(),size(1,3,1,1));
        test_assignment_helper(baz, bar);

        contiguous_state<4,T> qux(size2234(),size(4,4,4,4));
        test_assignment_helper(qux, bar);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( storage_order3, T, test_types )
{
    BOOST_TEST_MESSAGE("Instance without padding");
    {
        contiguous_state<3,T> foo(size(2,2,2));

        BOOST_CHECK_EQUAL( &(foo[0][0][0]) + 2*2, &(foo[1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   1, &(foo[0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   2, &(foo[0][0][1]));

        BOOST_CHECK_EQUAL( 2*2, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   1, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   2, foo.strides()[2] );
    }

    BOOST_TEST_MESSAGE("Instance with padding 0");
    {
        contiguous_state<3,T> foo(size(2,2,2), size(10,1,1));

        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +  10, &(foo[1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   1, &(foo[0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   2, &(foo[0][0][1]));

        BOOST_CHECK_EQUAL(  10, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   1, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   2, foo.strides()[2] );
    }

    BOOST_TEST_MESSAGE("Instance with padding 1");
    {
        contiguous_state<3,T> foo(size(2,2,2), size(1,2,1));

        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   8, &(foo[1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   2, &(foo[0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   4, &(foo[0][0][1]));

        BOOST_CHECK_EQUAL(   8, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   2, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   4, foo.strides()[2] );
    }

    BOOST_TEST_MESSAGE("Instance with padding 2");
    {
        contiguous_state<3,T> foo(size(2,2,2), size(1,1,3));

        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   6, &(foo[1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   1, &(foo[0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   3, &(foo[0][0][1]));

        BOOST_CHECK_EQUAL(   6, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   1, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   3, foo.strides()[2] );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( storage_order4, T, test_types )
{
    BOOST_TEST_MESSAGE("Instance without padding");
    {
        contiguous_state<4,T> foo(size(2,2,2,2));

        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) + 8, &(foo[1][0][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) + 1, &(foo[0][1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) + 2, &(foo[0][0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) + 4, &(foo[0][0][0][1]));

        BOOST_CHECK_EQUAL( 8, foo.strides()[0] );
        BOOST_CHECK_EQUAL( 1, foo.strides()[1] );
        BOOST_CHECK_EQUAL( 2, foo.strides()[2] );
        BOOST_CHECK_EQUAL( 4, foo.strides()[3] );
    }

    BOOST_TEST_MESSAGE("Instance with padding 0");
    {
        contiguous_state<4,T> foo(size(2,2,2,2), size(10,1,1,1));

        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  10, &(foo[1][0][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +   1, &(foo[0][1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +   2, &(foo[0][0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +   4, &(foo[0][0][0][1]));

        BOOST_CHECK_EQUAL(  10, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   1, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   2, foo.strides()[2] );
        BOOST_CHECK_EQUAL(   4, foo.strides()[3] );
    }

    BOOST_TEST_MESSAGE("Instance with padding 1");
    {
        contiguous_state<4,T> foo(size(2,2,2,2), size(1,2,1,1));

        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) + 16, &(foo[1][0][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  2, &(foo[0][1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  4, &(foo[0][0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  8, &(foo[0][0][0][1]));

        BOOST_CHECK_EQUAL(  16, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   2, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   4, foo.strides()[2] );
        BOOST_CHECK_EQUAL(   8, foo.strides()[3] );
    }

    BOOST_TEST_MESSAGE("Instance with padding 2");
    {
        contiguous_state<4,T> foo(size(2,2,2,2), size(1,1,3,1));

        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) + 12, &(foo[1][0][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  1, &(foo[0][1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  3, &(foo[0][0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  6, &(foo[0][0][0][1]));

        BOOST_CHECK_EQUAL(  12, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   1, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   3, foo.strides()[2] );
        BOOST_CHECK_EQUAL(   6, foo.strides()[3] );
    }

    BOOST_TEST_MESSAGE("Instance with padding 2");
    {
        contiguous_state<4,T> foo(size(2,2,2,2), size(1,1,1,5));

        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) + 10, &(foo[1][0][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  1, &(foo[0][1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  2, &(foo[0][0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0][0]) +  5, &(foo[0][0][0][1]));

        BOOST_CHECK_EQUAL(  10, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   1, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   2, foo.strides()[2] );
        BOOST_CHECK_EQUAL(   5, foo.strides()[3] );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( is_isomorphic3, T, test_types )
{
    contiguous_state<3,T> foo(  size(2,2,2));
    contiguous_state<3,T> bar(  size(2,2,2), size(7,1,1));
    contiguous_state<3,T> baz(  size(2,2,2), size(10,1,1));
    contiguous_state<3,T> qux(  size(1,2,2));
    contiguous_state<3,T> quux( size(2,1,2));
    contiguous_state<3,T> quuux(size(2,2,1));

    BOOST_CHECK_EQUAL(true,  foo.is_isomorphic(foo));
    BOOST_CHECK_EQUAL(true,  foo.is_isomorphic(bar));
    BOOST_CHECK_EQUAL(true,  foo.is_isomorphic(baz));
    BOOST_CHECK_EQUAL(true,  bar.is_isomorphic(foo));
    BOOST_CHECK_EQUAL(true,  bar.is_isomorphic(bar));
    BOOST_CHECK_EQUAL(true,  bar.is_isomorphic(baz));
    BOOST_CHECK_EQUAL(true,  baz.is_isomorphic(foo));
    BOOST_CHECK_EQUAL(true,  baz.is_isomorphic(bar));
    BOOST_CHECK_EQUAL(true,  baz.is_isomorphic(baz));

    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(qux));
    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(quux));
    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(quuux));
    BOOST_CHECK_EQUAL(false, bar.is_isomorphic(qux));
    BOOST_CHECK_EQUAL(false, bar.is_isomorphic(quux));
    BOOST_CHECK_EQUAL(false, bar.is_isomorphic(quuux));
    BOOST_CHECK_EQUAL(false, baz.is_isomorphic(qux));
    BOOST_CHECK_EQUAL(false, baz.is_isomorphic(quux));
    BOOST_CHECK_EQUAL(false, baz.is_isomorphic(quuux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( is_isomorphic4, T, test_types )
{
    contiguous_state<4,T> foo(  size(2,2,2,2));
    contiguous_state<4,T> bar(  size(2,2,2,2), size(37,1,1,15));
    contiguous_state<4,T> baz(  size(2,2,2,2), size(10,1,1,2));
    contiguous_state<4,T> qux(  size(1,2,2,2));
    contiguous_state<4,T> quux( size(2,1,2,2));
    contiguous_state<4,T> quuux(size(2,2,1,2));

    BOOST_CHECK_EQUAL(true,  foo.is_isomorphic(foo));
    BOOST_CHECK_EQUAL(true,  foo.is_isomorphic(bar));
    BOOST_CHECK_EQUAL(true,  foo.is_isomorphic(baz));
    BOOST_CHECK_EQUAL(true,  bar.is_isomorphic(foo));
    BOOST_CHECK_EQUAL(true,  bar.is_isomorphic(bar));
    BOOST_CHECK_EQUAL(true,  bar.is_isomorphic(baz));
    BOOST_CHECK_EQUAL(true,  baz.is_isomorphic(foo));
    BOOST_CHECK_EQUAL(true,  baz.is_isomorphic(bar));
    BOOST_CHECK_EQUAL(true,  baz.is_isomorphic(baz));

    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(qux));
    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(quux));
    BOOST_CHECK_EQUAL(false, foo.is_isomorphic(quuux));
    BOOST_CHECK_EQUAL(false, bar.is_isomorphic(qux));
    BOOST_CHECK_EQUAL(false, bar.is_isomorphic(quux));
    BOOST_CHECK_EQUAL(false, bar.is_isomorphic(quuux));
    BOOST_CHECK_EQUAL(false, baz.is_isomorphic(qux));
    BOOST_CHECK_EQUAL(false, baz.is_isomorphic(quux));
    BOOST_CHECK_EQUAL(false, baz.is_isomorphic(quuux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( scale3, T, test_types )
{
    BOOST_TEST_MESSAGE("Instance without padding");
    {
        contiguous_state<3,T> foo(size223());
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding 0");
    {
        contiguous_state<3,T> foo(size223(), size(7,1,1));
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding 1");
    {
        contiguous_state<3,T> foo(size223(), size(1,2,1));
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding 2");
    {
        contiguous_state<3,T> foo(size223(), size(1,1,3));
        test_scale_helper(foo);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( scale4, T, test_types )
{
    BOOST_TEST_MESSAGE("Instance without padding");
    {
        contiguous_state<4,T> foo(size2234());
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding 0");
    {
        contiguous_state<4,T> foo(size2234(), size(26,1,1,1));
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding 1");
    {
        contiguous_state<4,T> foo(size2234(), size(1,2,1,1));
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding 2");
    {
        contiguous_state<4,T> foo(size2234(), size(1,1,7,1));
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding 3");
    {
        contiguous_state<4,T> foo(size2234(), size(1,1,1,20));
        test_scale_helper(foo);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( add_scaled3, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        contiguous_state<3,T> foo(size223()), bar(size223());
        test_add_scaled_helper(foo, bar);

        // Operation between two nonconforming states throws
        contiguous_state<3,T> baz(size(2,2,2));
        BOOST_CHECK_THROW(foo.add_scaled(3, baz), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<3,T> foo(size223(),size(7,1,1));
        contiguous_state<3,T> bar(size223(),size(10,1,1));
        test_add_scaled_helper(foo, bar);

        contiguous_state<3,T> baz(size223(),size(1,2,1));
        test_add_scaled_helper(foo, baz);

        contiguous_state<3,T> qux(size223(),size(1,1,7));
        test_add_scaled_helper(foo, qux);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        contiguous_state<3,T> foo(size223(),size(7,1,1)), bar(size223());
        test_add_scaled_helper(foo, bar);

        contiguous_state<3,T> baz(size223(),size(1,2,1));
        test_add_scaled_helper(baz, bar);

        contiguous_state<3,T> qux(size223(),size(1,1,7));
        test_add_scaled_helper(qux, bar);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        contiguous_state<3,T> foo(size223()), bar(size223(),size(10,1,1));
        test_add_scaled_helper(foo, bar);

        contiguous_state<3,T> qux(size223(),size(1,9,25));
        test_add_scaled_helper(foo, qux);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( add_scaled4, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        contiguous_state<4,T> foo(size2234()), bar(size2234());
        test_add_scaled_helper(foo, bar);

        // Operation between two nonconforming states throws
        contiguous_state<4,T> baz(size(2,2,2,2));
        BOOST_CHECK_THROW(foo.add_scaled(3, baz), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<4,T> foo(size2234(),size(35,1,1,1));
        contiguous_state<4,T> bar(size2234(),size(37,1,1,1));
        test_add_scaled_helper(foo, bar);

        contiguous_state<4,T> baz(size2234(),size(1,2,1,1));
        test_add_scaled_helper(foo, baz);

        contiguous_state<4,T> qux(size2234(),size(1,1,7,1));
        test_add_scaled_helper(foo, qux);

        contiguous_state<4,T> quux(size2234(),size(1,1,1,15));
        test_add_scaled_helper(foo, quux);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        contiguous_state<4,T> foo(size2234(),size(63,1,1,1));
        contiguous_state<4,T> bar(size2234());
        test_add_scaled_helper(foo, bar);

        contiguous_state<4,T> baz(size2234(),size(1,2,1,1));
        test_add_scaled_helper(baz, bar);

        contiguous_state<4,T> qux(size2234(),size(1,1,7,1));
        test_add_scaled_helper(qux, bar);

        contiguous_state<4,T> quux(size2234(),size(1,1,1,23));
        test_add_scaled_helper(quux, bar);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        contiguous_state<4,T> foo(size2234());
        contiguous_state<4,T> bar(size2234(),size(41,1,1,1));
        test_add_scaled_helper(foo, bar);

        contiguous_state<4,T> qux(size2234(),size(1,9,25,1));
        test_add_scaled_helper(foo, qux);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange3, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        contiguous_state<3,T> foo(size223()), bar(size223());
        test_exchange_helper(foo, bar);

        // Operation between two nonconforming states throws
        contiguous_state<3,T> baz(size(2,2,2));
        BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<3,T> foo(size223(),size(11,1,1));
        contiguous_state<3,T> bar(size223(),size(13,1,1));
        test_exchange_helper(foo, bar);

        contiguous_state<3,T> baz(size223(),size(1,4,1));
        test_exchange_helper(foo, baz);

        contiguous_state<3,T> qux(size223(),size(25,1,5));
        test_exchange_helper(foo, qux);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        contiguous_state<3,T> foo(size223(),size(9,1,1)), bar(size223());
        test_exchange_helper(foo, bar);

        contiguous_state<3,T> baz(size223(),size(1,4,10)), qux(size223());
        test_exchange_helper(baz, qux);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        contiguous_state<3,T> foo(size223()), bar(size223(),size(17,1,1));
        test_exchange_helper(foo, bar);

        contiguous_state<3,T> baz(size223()), qux(size223(),size(11,17,1));
        test_exchange_helper(baz, qux);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange4, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        contiguous_state<4,T> foo(size2234()), bar(size2234());
        test_exchange_helper(foo, bar);

        // Operation between two nonconforming states throws
        contiguous_state<4,T> baz(size(2,2,2,2));
        BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<4,T> foo(size2234(),size(53,1,1,1));
        contiguous_state<4,T> bar(size2234(),size(57,1,1,1));
        test_exchange_helper(foo, bar);

        contiguous_state<4,T> baz(size2234(),size(1,4,1,1));
        test_exchange_helper(foo, baz);

        contiguous_state<4,T> qux(size2234(),size(59,1,5,12));
        test_exchange_helper(foo, qux);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        contiguous_state<4,T> foo(size2234(),size(27,1,1,1));
        contiguous_state<4,T> bar(size2234());
        test_exchange_helper(foo, bar);

        contiguous_state<4,T> baz(size2234(),size(1,4,10,1));
        contiguous_state<4,T> qux(size2234());
        test_exchange_helper(baz, qux);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        contiguous_state<4,T> foo(size2234());
        contiguous_state<4,T> bar(size2234(),size(67,1,1,1));
        test_exchange_helper(foo, bar);

        contiguous_state<4,T> baz(size2234());
        contiguous_state<4,T> qux(size2234(),size(33,17,1,8));
        test_exchange_helper(baz, qux);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( concept_check3, T, test_types )
{
    using boost::detail::multi_array::MutableMultiArrayConcept;
    BOOST_CONCEPT_ASSERT((MutableMultiArrayConcept<contiguous_state<3,T>,3>));
    BOOST_CHECK(true); // Avoids "did not run any assertions" message
}

BOOST_AUTO_TEST_CASE_TEMPLATE( concept_check4, T, test_types )
{
    using boost::detail::multi_array::MutableMultiArrayConcept;
    BOOST_CONCEPT_ASSERT((MutableMultiArrayConcept<contiguous_state<4,T>,4>));
    BOOST_CHECK(true); // Avoids "did not run any assertions" message
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( Contiguous_interoperate_Interleaved )

using suzerain::contiguous_state;
using suzerain::interleaved_state;

// Types to be tested with both state_base subclasses
typedef boost::mpl::list<
    double
    ,float
    ,std::complex<double>
    ,std::complex<float>
> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, T, test_types )
{
    BOOST_TEST_MESSAGE("Instances without padding");
    {
        contiguous_state<4,T>  foo(size2234());
        interleaved_state<4,T> bar(size2234());
        test_assignment_helper(foo, bar);
        test_assignment_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("First instance with padding");
    {
        contiguous_state<4,T>  foo(size2234(), size(27,1,1,1));
        interleaved_state<4,T> bar(size2234());
        test_assignment_helper(foo, bar);
        test_assignment_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("Second instance with padding");
    {
        contiguous_state<4,T>  foo(size2234());
        interleaved_state<4,T> bar(size2234(), 87);
        test_assignment_helper(foo, bar);
        test_assignment_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<4,T>  foo(size2234(), size(27,1,1,1));
        interleaved_state<4,T> bar(size2234(), 87);
        test_assignment_helper(foo, bar);
        test_assignment_helper(bar, foo);
    }

    // Operation between two nonconforming states throws
    contiguous_state<4,T>  foo(size2234());
    interleaved_state<4,T> baz(size(2,2,2,2));
    BOOST_CHECK_THROW(baz.assign_from(foo), std::logic_error);
    BOOST_CHECK_THROW(foo.assign_from(baz), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( is_isomorphic, T, test_types )
{
    interleaved_state<4,T> i_foo(  size(2,2,2,2));
    interleaved_state<4,T> i_bar(  size(2,2,2,2));
    interleaved_state<4,T> i_baz(  size(1,2,2,2));
    interleaved_state<4,T> i_qux(  size(2,1,2,2));
    interleaved_state<4,T> i_quux( size(2,2,1,2));
    interleaved_state<4,T> i_quuux(size(2,2,2,1));
    contiguous_state<4,T>  c_foo(  size(2,2,2,2));
    contiguous_state<4,T>  c_bar(  size(2,2,2,2));
    contiguous_state<4,T>  c_baz(  size(1,2,2,2));
    contiguous_state<4,T>  c_qux(  size(2,1,2,2));
    contiguous_state<4,T>  c_quux( size(2,2,1,2));
    contiguous_state<4,T>  c_quuux(size(2,2,2,1));

    BOOST_CHECK_EQUAL(true,  i_foo.is_isomorphic(c_foo));
    BOOST_CHECK_EQUAL(true,  i_foo.is_isomorphic(c_bar));
    BOOST_CHECK_EQUAL(false, i_foo.is_isomorphic(c_baz));
    BOOST_CHECK_EQUAL(false, i_foo.is_isomorphic(c_qux));
    BOOST_CHECK_EQUAL(false, i_foo.is_isomorphic(c_quux));
    BOOST_CHECK_EQUAL(false, i_foo.is_isomorphic(c_quuux));
    BOOST_CHECK_EQUAL(true,  c_foo.is_isomorphic(i_foo));
    BOOST_CHECK_EQUAL(true,  c_foo.is_isomorphic(i_bar));
    BOOST_CHECK_EQUAL(false, c_foo.is_isomorphic(i_baz));
    BOOST_CHECK_EQUAL(false, c_foo.is_isomorphic(i_qux));
    BOOST_CHECK_EQUAL(false, c_foo.is_isomorphic(i_quux));
    BOOST_CHECK_EQUAL(false, c_foo.is_isomorphic(i_quuux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( add_scaled, T, test_types )
{
    BOOST_TEST_MESSAGE("Instances without padding");
    {
        contiguous_state<4,T>  foo(size2234());
        interleaved_state<4,T> bar(size2234());
        test_add_scaled_helper(foo, bar);
        test_add_scaled_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("First instance with padding");
    {
        contiguous_state<4,T>  foo(size2234(), size(27,1,1,1));
        interleaved_state<4,T> bar(size2234());
        test_add_scaled_helper(foo, bar);
        test_add_scaled_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("Second instance with padding");
    {
        contiguous_state<4,T>  foo(size2234());
        interleaved_state<4,T> bar(size2234(), 87);
        test_add_scaled_helper(foo, bar);
        test_add_scaled_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<4,T>  foo(size2234(), size(27,1,1,1));
        interleaved_state<4,T> bar(size2234(), 87);
        test_add_scaled_helper(foo, bar);
        test_add_scaled_helper(bar, foo);
    }

    // Operation between two nonconforming states throws
    contiguous_state<4,T>  foo(size2234());
    interleaved_state<4,T> baz(size(2,2,2,2));
    BOOST_CHECK_THROW(foo.add_scaled(3, baz), std::logic_error);
    BOOST_CHECK_THROW(baz.add_scaled(3, foo), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange, T, test_types )
{
    BOOST_TEST_MESSAGE("Instances without padding");
    {
        contiguous_state<4,T>  foo(size2234());
        interleaved_state<4,T> bar(size2234());
        test_exchange_helper(foo, bar);
        test_exchange_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("First instance with padding");
    {
        contiguous_state<4,T>  foo(size2234(), size(27,1,1,1));
        interleaved_state<4,T> bar(size2234());
        test_exchange_helper(foo, bar);
        test_exchange_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("Second instance with padding");
    {
        contiguous_state<4,T>  foo(size2234());
        interleaved_state<4,T> bar(size2234(), 87);
        test_exchange_helper(foo, bar);
        test_exchange_helper(bar, foo);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        contiguous_state<4,T>  foo(size2234(), size(27,1,1,1));
        interleaved_state<4,T> bar(size2234(), 87);
        test_exchange_helper(foo, bar);
        test_exchange_helper(bar, foo);
    }

    // Operation between two nonconforming states throws
    contiguous_state<4,T>  foo(size2234());
    interleaved_state<4,T> baz(size(2,2,2,2));
    BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);
    BOOST_CHECK_THROW(baz.exchange(foo), std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()
