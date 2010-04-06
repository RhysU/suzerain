#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/state_impl.hpp>
#include <boost/concept/assert.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

#pragma warning(disable:383)

// Helper for specifying extents information
static boost::array<std::size_t,3> size3(
    std::size_t x, std::size_t y, std::size_t z)
{
    boost::array<std::size_t,3> a = { x, y, z };
    return a;
}

static boost::array<std::size_t,3> size223() {
    return size3(2,2,3);
}

template<
    template <std::size_t,typename,typename> class State,
    typename Allocator,
    typename FPT,
    typename Scale
>
static void load223(
    State<3,FPT,Allocator> &state,
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

template<
    template <std::size_t,typename,typename> class State,
    typename Allocator,
    typename FPT,
    typename Scale
>
static void verify223(
    const State<3,FPT,Allocator> &state,
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

template<
    template <std::size_t,typename,typename> class State,
    typename Allocator,
    typename FPT,
    typename Scale
>
static void load223(
    State<3,std::complex<FPT>,Allocator> &state,
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

template<
    template <std::size_t,typename,typename> class State,
    typename Allocator,
    typename FPT,
    typename Scale
>
static void verify223(
    const State<3,std::complex<FPT>,Allocator> &state,
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

template<
    template <std::size_t,typename,typename> class State1,
    template <std::size_t,typename,typename> class State2,
    typename Allocator1,
    typename Allocator2,
    typename Element
>
static void test_assignment_helper(
    State1<3,Element,Allocator1> &foo,
    State2<3,Element,Allocator2> &bar)
{
    load223(foo, 1);

    foo.assign(foo); // Self
    verify223(foo, 1);

    bar.assign(foo);
    verify223(bar, 1);
}

template<
    template <std::size_t,typename,typename> class State,
    typename Allocator,
    typename Element
>
static void test_scale_helper(State<3,Element,Allocator> &foo)
{
    load223(foo, 1);
    verify223(foo, 1);

    foo.scale(1);
    verify223(foo, 1);

    foo.scale(2);
    verify223(foo, 2);

    foo.scale(0);
    verify223(foo, 0);
}

template<
    template <std::size_t,typename,typename> class State1,
    template <std::size_t,typename,typename> class State2,
    typename Allocator1,
    typename Allocator2,
    typename Element
>
static void test_addScaled_helper(
    State1<3,Element,Allocator1> &foo,
    State2<3,Element,Allocator2> &bar)
{
    load223(foo, 1);
    load223(bar, 2);
    foo.addScaled(3, bar);
    verify223(foo, 7);
}

template<
    template <std::size_t,typename,typename> class State1,
    template <std::size_t,typename,typename> class State2,
    typename Allocator1,
    typename Allocator2,
    typename Element
>
static void test_exchange_helper(
    State1<3,Element,Allocator1> &foo,
    State2<3,Element,Allocator2> &bar)
{
    load223(foo, 1);
    load223(bar, 2);
    foo.exchange(bar);
    verify223(foo, 2);
    verify223(bar, 1);
    bar.exchange(foo);
    verify223(foo, 1);
    verify223(bar, 2);
}

// FIXME Correct InterleavedState implementation and enable tests
#ifdef FIXME_INTERLEAVED_STATE_TESTS_DISABLED
#warning 'FIXME: suzerain::InterleavedState tests disabled'
BOOST_AUTO_TEST_SUITE( InterleavedState )

using suzerain::InterleavedState;

// Types to be tested with InterleavedState
typedef boost::mpl::list<
    double
   ,float
   ,std::complex<double>
   ,std::complex<float>
> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors, T, test_types )
{
    // Regular constructor
    InterleavedState<3,T> foo(size223());
    BOOST_CHECK_EQUAL(
          std::distance(foo.memory_begin(),foo.memory_end()), 2*2*3);
    load223(foo, 1);
    verify223(foo, 1);

    // Copy construct a second instance from the first
    InterleavedState<3,T> bar(foo);
    BOOST_CHECK_EQUAL(
          std::distance(bar.memory_begin(),bar.memory_end()), 2*2*3);
    verify223(bar, 1);

    // Modify first instance's data
    for (int i = 0; i < foo.shape()[0]; ++i)
        for (int j = 0; j < foo.shape()[1]; ++j)
            for (int k = 0; k < foo.shape()[2]; ++k)
                foo[i][j][k] += (i+1)*(j+1)*(k+1);

    // Ensure copy constructed data in second instance not modified
    verify223(bar, 1);

    // Create padded instance and ensure content lies within padding
    InterleavedState<3,T> baz(size223(), 2*2*3*7);
    BOOST_CHECK_EQUAL(
          std::distance(baz.memory_begin(),baz.memory_end()), 2*2*3*7);
    BOOST_CHECK_GE(baz.memory_begin(), &(baz[0][0][0]));
    BOOST_CHECK_LT(&(baz[1][1][2]), baz.memory_end());

    // Ensure padded information present propagated in copy operations
    InterleavedState<3,T> qux(baz);
    BOOST_CHECK_EQUAL(
          std::distance(qux.memory_begin(),qux.memory_end()), 2*2*3*7);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, T, test_types )
{
    InterleavedState<3,T> foo(size223()), bar(size223());
    test_assignment_helper(foo, bar);

    // Operation between two nonconforming states throws
    InterleavedState<3,T> baz(size3(2,2,2));
    BOOST_CHECK_THROW(baz.assign(foo), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( storage_order, T, test_types )
{
    InterleavedState<3,T> foo(size3(2,2,2));

    BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   1, &(foo[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   2, &(foo[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo[0][0][0]) + 2*2, &(foo[0][0][1]));

    BOOST_CHECK_EQUAL(   1, foo.strides()[0] );
    BOOST_CHECK_EQUAL(   2, foo.strides()[1] );
    BOOST_CHECK_EQUAL( 2*2, foo.strides()[2] );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( isConformant, T, test_types )
{
    InterleavedState<3,T> foo( size3(2,2,2));
    InterleavedState<3,T> bar( size3(2,2,2));
    InterleavedState<3,T> baz( size3(1,2,2));
    InterleavedState<3,T> qux( size3(2,1,2));
    InterleavedState<3,T> quux(size3(2,2,1));

    BOOST_CHECK_EQUAL(true,  foo.isConformant(foo));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(bar));
    BOOST_CHECK_EQUAL(false, foo.isConformant(baz));
    BOOST_CHECK_EQUAL(false, foo.isConformant(qux));
    BOOST_CHECK_EQUAL(false, foo.isConformant(quux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( scale, T, test_types )
{
    InterleavedState<3,T> foo(size223());
    test_scale_helper(foo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( addScaled, T, test_types )
{
    InterleavedState<3,T> foo(size223()), bar(size223());
    test_addScaled_helper(foo, bar);

    // Operation between two nonconforming states throws
    InterleavedState<3,T> baz(size3(2,2,2));
    BOOST_CHECK_THROW(foo.addScaled(3, baz), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange, T, test_types )
{
    InterleavedState<3,T> foo(size223()), bar(size223());
    test_exchange_helper(foo, bar);

    // Operation between two nonconforming states throws
    InterleavedState<3,T> baz(size3(2,2,2));
    BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( concept_check, T, test_types )
{
    using boost::detail::multi_array::MutableMultiArrayConcept;
    BOOST_CONCEPT_ASSERT((MutableMultiArrayConcept<InterleavedState<3,T>,3>));
}

BOOST_AUTO_TEST_SUITE_END()
#endif /* FIXME_INTERLEAVED_STATE_TESTS_DISABLED */


BOOST_AUTO_TEST_SUITE( NoninterleavedState )

using suzerain::NoninterleavedState;

// Types to be tested with NoninterleavedState
typedef boost::mpl::list<
    double
   ,float
   ,std::complex<double>
   ,std::complex<float>
> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors, T, test_types )
{
    // Regular constructor
    NoninterleavedState<3,T> foo(size223());
    BOOST_CHECK_EQUAL(
          std::distance(foo.memory_begin(),foo.memory_end()), 2*2*3);
    load223(foo, 1);
    verify223(foo, 1);

    // Copy construct a second instance from the first
    NoninterleavedState<3,T> bar(foo);
    BOOST_CHECK_EQUAL(
          std::distance(bar.memory_begin(),bar.memory_end()), 2*2*3);
    verify223(bar, 1);

    // Modify first instance's data
    for (int i = 0; i < foo.shape()[0]; ++i)
        for (int j = 0; j < foo.shape()[1]; ++j)
            for (int k = 0; k < foo.shape()[2]; ++k)
                foo[i][j][k] += (i+1)*(j+1)*(k+1);

    // Ensure copy constructed data in second instance not modified
    verify223(bar, 1);

    // Create padded instance and ensure content lies within padding
    NoninterleavedState<3,T> baz(size223(), size3(7,1,1));
    BOOST_CHECK_EQUAL(
          std::distance(baz.memory_begin(),baz.memory_end()), 2*7);
    BOOST_CHECK_EQUAL(baz.strides()[0], 7);
    BOOST_CHECK_GE(baz.memory_begin(), &(baz[0][0][0]));
    BOOST_CHECK_LT(&(baz[1][1][2]), baz.memory_end());
    load223(baz, 1);
    verify223(baz, 1);

    // Ensure padded information present propagated in copy operations
    NoninterleavedState<3,T> qux(baz);
    BOOST_CHECK_EQUAL(
          std::distance(qux.memory_begin(),qux.memory_end()), 2*7);
    BOOST_CHECK_EQUAL(baz.strides()[0], qux.strides()[0]);
    BOOST_CHECK_EQUAL(baz.strides()[1], qux.strides()[1]);
    BOOST_CHECK_EQUAL(baz.strides()[2], qux.strides()[2]);
    verify223(qux, 1);

    // Modify first padded instance's data
    for (int i = 0; i < foo.shape()[0]; ++i)
        for (int j = 0; j < foo.shape()[1]; ++j)
            for (int k = 0; k < foo.shape()[2]; ++k)
                baz[i][j][k] += (i+1)*(j+1)*(k+1);

    // Ensure copy constructed data in second padded instance not modified
    verify223(qux, 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        NoninterleavedState<3,T> foo(size223()), bar(size223());
        test_assignment_helper(foo, bar);

        // Operation between two nonconforming states throws
        NoninterleavedState<3,T> baz(size3(2,2,2));
        BOOST_CHECK_THROW(baz.assign(foo), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        NoninterleavedState<3,T> foo(size223(),size3(7,1,1));
        NoninterleavedState<3,T> bar(size223(),size3(10,1,1));
        test_assignment_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        NoninterleavedState<3,T> foo(size223()), bar(size223(),size3(7,1,1));
        test_assignment_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        NoninterleavedState<3,T> foo(size223(),size3(7,1,1)), bar(size223());
        test_assignment_helper(foo, bar);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( storage_order, T, test_types )
{
    BOOST_TEST_MESSAGE("Instance without padding");
    {
        NoninterleavedState<3,T> foo(size3(2,2,2));

        BOOST_CHECK_EQUAL( &(foo[0][0][0]) + 2*2, &(foo[1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   1, &(foo[0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   2, &(foo[0][0][1]));

        BOOST_CHECK_EQUAL( 2*2, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   1, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   2, foo.strides()[2] );
    }

    BOOST_TEST_MESSAGE("Instance with padding");
    {
        NoninterleavedState<3,T> foo(size3(2,2,2), size3(10,1,1));

        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +  10, &(foo[1][0][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   1, &(foo[0][1][0]));
        BOOST_CHECK_EQUAL( &(foo[0][0][0]) +   2, &(foo[0][0][1]));

        BOOST_CHECK_EQUAL(  10, foo.strides()[0] );
        BOOST_CHECK_EQUAL(   1, foo.strides()[1] );
        BOOST_CHECK_EQUAL(   2, foo.strides()[2] );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( isConformant, T, test_types )
{
    NoninterleavedState<3,T> foo(  size3(2,2,2));
    NoninterleavedState<3,T> bar(  size3(2,2,2), size3(7,1,1));
    NoninterleavedState<3,T> baz(  size3(2,2,2), size3(10,1,1));
    NoninterleavedState<3,T> qux(  size3(1,2,2));
    NoninterleavedState<3,T> quux( size3(2,1,2));
    NoninterleavedState<3,T> quuux(size3(2,2,1));

    BOOST_CHECK_EQUAL(true,  foo.isConformant(foo));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(bar));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(baz));
    BOOST_CHECK_EQUAL(true,  bar.isConformant(foo));
    BOOST_CHECK_EQUAL(true,  bar.isConformant(bar));
    BOOST_CHECK_EQUAL(true,  bar.isConformant(baz));
    BOOST_CHECK_EQUAL(true,  baz.isConformant(foo));
    BOOST_CHECK_EQUAL(true,  baz.isConformant(bar));
    BOOST_CHECK_EQUAL(true,  baz.isConformant(baz));

    BOOST_CHECK_EQUAL(false, foo.isConformant(qux));
    BOOST_CHECK_EQUAL(false, foo.isConformant(quux));
    BOOST_CHECK_EQUAL(false, foo.isConformant(quuux));
    BOOST_CHECK_EQUAL(false, bar.isConformant(qux));
    BOOST_CHECK_EQUAL(false, bar.isConformant(quux));
    BOOST_CHECK_EQUAL(false, bar.isConformant(quuux));
    BOOST_CHECK_EQUAL(false, baz.isConformant(qux));
    BOOST_CHECK_EQUAL(false, baz.isConformant(quux));
    BOOST_CHECK_EQUAL(false, baz.isConformant(quuux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( scale, T, test_types )
{
    BOOST_TEST_MESSAGE("Instance without padding");
    {
        NoninterleavedState<3,T> foo(size223());
        test_scale_helper(foo);
    }

    BOOST_TEST_MESSAGE("Instance with padding");
    {
        NoninterleavedState<3,T> foo(size223(), size3(7,1,1));
        test_scale_helper(foo);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( addScaled, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        NoninterleavedState<3,T> foo(size223()), bar(size223());
        test_addScaled_helper(foo, bar);

        // Operation between two nonconforming states throws
        NoninterleavedState<3,T> baz(size3(2,2,2));
        BOOST_CHECK_THROW(foo.addScaled(3, baz), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        NoninterleavedState<3,T> foo(size223(),size3(7,1,1));
        NoninterleavedState<3,T> bar(size223(),size3(10,1,1));
        test_addScaled_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        NoninterleavedState<3,T> foo(size223(),size3(7,1,1)), bar(size223());
        test_addScaled_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        NoninterleavedState<3,T> foo(size223()), bar(size223(),size3(10,1,1));
        test_addScaled_helper(foo, bar);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange, T, test_types )
{
    BOOST_TEST_MESSAGE("Both instances without padding");
    {
        NoninterleavedState<3,T> foo(size223()), bar(size223());
        test_exchange_helper(foo, bar);

        // Operation between two nonconforming states throws
        NoninterleavedState<3,T> baz(size3(2,2,2));
        BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);
    }

    BOOST_TEST_MESSAGE("Both instances with padding");
    {
        NoninterleavedState<3,T> foo(size223(),size3(11,1,1));
        NoninterleavedState<3,T> bar(size223(),size3(13,1,1));
        test_exchange_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Target instance with padding");
    {
        NoninterleavedState<3,T> foo(size223(),size3(9,1,1)), bar(size223());
        test_exchange_helper(foo, bar);
    }

    BOOST_TEST_MESSAGE("Source instance with padding");
    {
        NoninterleavedState<3,T> foo(size223()), bar(size223(),size3(17,1,1));
        test_exchange_helper(foo, bar);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( concept_check, T, test_types )
{
    using boost::detail::multi_array::MutableMultiArrayConcept;
    BOOST_CONCEPT_ASSERT((MutableMultiArrayConcept<NoninterleavedState<3,T>,3>));
}

BOOST_AUTO_TEST_SUITE_END()
