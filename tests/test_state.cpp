#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/state_impl.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

#pragma warning(disable:383)

// Types to be tested with InterleavedState
typedef boost::mpl::list<
    double
   ,float
   ,std::complex<double>
   ,std::complex<float>
> test_types;

BOOST_AUTO_TEST_SUITE( InterleavedState )

static boost::array<std::size_t,3> size3(
      std::size_t x, std::size_t y, std::size_t z)
{
   boost::array<std::size_t,3> a = { x, y, z };
   return a;
}

template< typename FPT, typename Scale >
static void load223(
      typename suzerain::InterleavedState<3,FPT> &state,
      const Scale scale,
      typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3);

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

template< typename FPT, typename Scale >
static void verify223(
      const typename suzerain::InterleavedState<3,FPT> &state,
      const Scale scale,
      typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3);

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

template< typename FPT, typename Scale >
static void load223(
      typename suzerain::InterleavedState<3,std::complex<FPT> > &state,
      const Scale scale,
      typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3);

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

template< typename FPT, typename Scale >
static void verify223(
      const typename suzerain::InterleavedState<3,std::complex<FPT> > &state,
      const Scale scale,
      typename boost::enable_if<boost::is_floating_point<FPT> >::type *dummy = 0)
{
    SUZERAIN_UNUSED(dummy);
    BOOST_REQUIRE_EQUAL(state.shape()[0], 2);
    BOOST_REQUIRE_EQUAL(state.shape()[1], 2);
    BOOST_REQUIRE_EQUAL(state.shape()[2], 3);

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

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors, T, test_types )
{
    // Regular constructor
    suzerain::InterleavedState<3,T> foo(size3(2,2,3));
    BOOST_CHECK_EQUAL(foo.raw_memory_count(), 2*2*3);
    load223(foo, 1);
    verify223(foo, 1);

    // Copy construct a second instance from the first
    suzerain::InterleavedState<3,T> bar(foo);
    BOOST_CHECK_EQUAL(bar.raw_memory_count(), 2*2*3);
    verify223(bar, 1);

    // Modify first instance's data
    for (int i = 0; i < foo.shape()[0]; ++i)
        for (int j = 0; j < foo.shape()[1]; ++j)
            for (int k = 0; k < foo.shape()[2]; ++k)
                foo[i][j][k] += (i+1)*(j+1)*(k+1);

    // Ensure copy constructed data in second instance not modified
    verify223(bar, 1);

    // Create padded instance and ensure content lies within padding
    suzerain::InterleavedState<3,T> baz(size3(2,2,3), 2*2*3*7);
    BOOST_CHECK_EQUAL(baz.raw_memory_count(), 2*2*3*7);
    BOOST_CHECK_GE(baz.raw_memory(), &(baz[0][0][0]));
    BOOST_CHECK_LT(&(baz[1][1][2]), baz.raw_memory() + baz.raw_memory_count());

    // Ensure padded information present propagated in copy operations
    suzerain::InterleavedState<3,T> qux(baz);
    BOOST_CHECK_EQUAL(qux.raw_memory_count(), 2*2*3*7);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, T, test_types )
{
    suzerain::InterleavedState<3,T> foo(size3(2,2,3)), bar(size3(2,2,3));
    load223(foo, 1);

    foo.assign(foo); // Self
    verify223(foo, 1);

    bar.assign(foo);
    verify223(bar, 1);

    // Operation between two nonconforming states throws
    suzerain::InterleavedState<3,T> baz(size3(2,2,2));
    BOOST_CHECK_THROW(baz.assign(foo), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( storage_order, T, test_types )
{
    suzerain::InterleavedState<3,T> foo(size3(2,2,2));

    BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   1, &(foo[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   2, &(foo[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo[0][0][0]) + 2*2, &(foo[0][0][1]));

    BOOST_CHECK_EQUAL(   1, foo.strides()[0] );
    BOOST_CHECK_EQUAL(   2, foo.strides()[1] );
    BOOST_CHECK_EQUAL( 2*2, foo.strides()[2] );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( isConformant, T, test_types )
{
    suzerain::InterleavedState<3,T> foo( size3(2,2,2));
    suzerain::InterleavedState<3,T> bar( size3(2,2,2));
    suzerain::InterleavedState<3,T> baz( size3(1,2,2));
    suzerain::InterleavedState<3,T> qux( size3(2,1,2));
    suzerain::InterleavedState<3,T> quux(size3(2,2,1));

    BOOST_CHECK_EQUAL(true,  foo.isConformant(foo));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(bar));
    BOOST_CHECK_EQUAL(false, foo.isConformant(baz));
    BOOST_CHECK_EQUAL(false, foo.isConformant(qux));
    BOOST_CHECK_EQUAL(false, foo.isConformant(quux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( scale, T, test_types )
{
    suzerain::InterleavedState<3,T> foo(size3(2,2,3));
    load223(foo, 1);
    verify223(foo, 1);

    foo.scale(1);
    verify223(foo, 1);

    foo.scale(2);
    verify223(foo, 2);

    foo.scale(0);
    verify223(foo, 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( addScaled, T, test_types )
{
    suzerain::InterleavedState<3,T> foo(size3(2,2,3)), bar(size3(2,2,3));
    load223(foo, 1);
    load223(bar, 2);
    foo.addScaled(3, bar);
    verify223(foo, 7);

    // Operation between two nonconforming states throws
    suzerain::InterleavedState<3,T> baz(size3(2,2,2));
    BOOST_CHECK_THROW(foo.addScaled(3, baz), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange, T, test_types )
{
    suzerain::InterleavedState<3,T> foo(size3(2,2,3)), bar(size3(2,2,3));
    load223(foo, 1);
    load223(bar, 2);
    foo.exchange(bar);
    verify223(foo, 2);
    verify223(bar, 1);
    bar.exchange(foo);
    verify223(foo, 1);
    verify223(bar, 2);

    // Operation between two nonconforming states throws
    suzerain::InterleavedState<3,T> baz(size3(2,2,2));
    BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( comparison, T, test_types )
{
    suzerain::InterleavedState<3,T> foo(size3(2,2,3)), bar(size3(2,2,3));
    load223(foo, 1);
    load223(bar, 1);

    BOOST_CHECK(foo == bar);
    BOOST_CHECK(bar == foo);
    foo[0][0][0] += foo[0][0][0];
    BOOST_CHECK(foo != bar);
    BOOST_CHECK(bar != foo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( concept_check, T, test_types )
{
    // FIXME Add MultiArray Concept checking for several types...
}

BOOST_AUTO_TEST_SUITE_END()
