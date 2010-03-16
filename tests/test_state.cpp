#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/state.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

#pragma warning(disable:383)

typedef boost::mpl::list_c<bool,true,false> bool_values;

BOOST_AUTO_TEST_SUITE( RealState )

BOOST_AUTO_TEST_CASE_TEMPLATE( declare_pointer, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::RealState<double,Interleaved::value> *state;
    (void)state;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    // Regular constructor
    suzerain::RealState<double,Interleaved::value> foo(2, 2, 3);
    BOOST_CHECK_EQUAL(foo.variable_count, 2);
    BOOST_CHECK_EQUAL(foo.vector_length, 2);
    BOOST_CHECK_EQUAL(foo.vector_count, 3);

    foo.data[0][0][0] =   2.0;
    foo.data[0][1][0] =   3.0;
    foo.data[0][0][1] =   5.0;
    foo.data[0][1][1] =   7.0;
    foo.data[0][0][2] =  11.0;
    foo.data[0][1][2] =  13.0;
    foo.data[1][0][0] = 102.0;
    foo.data[1][1][0] = 103.0;
    foo.data[1][0][1] = 105.0;
    foo.data[1][1][1] = 107.0;
    foo.data[1][0][2] = 111.0;
    foo.data[1][1][2] = 113.0;

    BOOST_CHECK_EQUAL(foo.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][2], 113.0);

    // Copy construct a second instance
    suzerain::RealState<double,Interleaved::value> bar(foo);
    BOOST_CHECK_EQUAL(bar.variable_count, 2);
    BOOST_CHECK_EQUAL(bar.vector_length, 2);
    BOOST_CHECK_EQUAL(bar.vector_count, 3);

    // Modify first instance's data
    foo.data[0][0][0] +=   2.0;
    foo.data[0][1][0] +=   3.0;
    foo.data[0][0][1] +=   5.0;
    foo.data[0][1][1] +=   7.0;
    foo.data[0][0][2] +=  11.0;
    foo.data[0][1][2] +=  13.0;
    foo.data[1][0][0] += 102.0;
    foo.data[1][1][0] += 103.0;
    foo.data[1][0][1] += 105.0;
    foo.data[1][1][1] += 107.0;
    foo.data[1][0][2] += 111.0;
    foo.data[1][1][2] += 113.0;

    // Ensure copy constructed data in second instance not modified
    BOOST_CHECK_EQUAL(bar.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][2], 113.0);

    // Copy construct a third instance with a different interleaving
    suzerain::RealState<double,!Interleaved::value> baz(bar);
    BOOST_CHECK_NE(bar.interleaved(), baz.interleaved());
    BOOST_CHECK_EQUAL(baz.variable_count, 2);
    BOOST_CHECK_EQUAL(baz.vector_length, 2);
    BOOST_CHECK_EQUAL(baz.vector_count, 3);

    // Modify second instance's data
    bar.data[0][0][0] +=   2.0;
    bar.data[0][1][0] +=   3.0;
    bar.data[0][0][1] +=   5.0;
    bar.data[0][1][1] +=   7.0;
    bar.data[0][0][2] +=  11.0;
    bar.data[0][1][2] +=  13.0;
    bar.data[1][0][0] += 102.0;
    bar.data[1][1][0] += 103.0;
    bar.data[1][0][1] += 105.0;
    bar.data[1][1][1] += 107.0;
    bar.data[1][0][2] += 111.0;
    bar.data[1][1][2] += 113.0;

    // Ensure copy constructed data in third instance not modified
    BOOST_CHECK_EQUAL(baz.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(baz.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(baz.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(baz.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(baz.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(baz.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(baz.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(baz.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(baz.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(baz.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(baz.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(baz.data[1][1][2], 113.0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::RealState<double,Interleaved::value> foo(2, 2, 3), bar(2,2,3);

    foo.data[0][0][0] =   2.0;
    foo.data[0][1][0] =   3.0;
    foo.data[0][0][1] =   5.0;
    foo.data[0][1][1] =   7.0;
    foo.data[0][0][2] =  11.0;
    foo.data[0][1][2] =  13.0;
    foo.data[1][0][0] = 102.0;
    foo.data[1][1][0] = 103.0;
    foo.data[1][0][1] = 105.0;
    foo.data[1][1][1] = 107.0;
    foo.data[1][0][2] = 111.0;
    foo.data[1][1][2] = 113.0;

    // Self assignment
    foo = foo;

    BOOST_CHECK_EQUAL(foo.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][2], 113.0);

    // Normal assignment
    bar = foo;

    BOOST_CHECK_EQUAL(bar.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][2], 113.0);

    // Ensure we catch an operation between two nonconforming states
    suzerain::RealState<double,Interleaved::value> baz(2, 2, 2);
    BOOST_CHECK_THROW(baz = foo, std::logic_error);

    // Ensure we catch an operation between two different subclasses
    BOOST_CHECK_THROW(baz = suzerain::ComplexState<double>(1,2,3),
                      std::bad_cast);

    // Assignment from different interleaving
    suzerain::RealState<double,!Interleaved::value> qux(2, 2, 3);
    BOOST_CHECK_NE(qux.interleaved(), foo.interleaved());
    qux = foo;
    BOOST_CHECK_EQUAL(qux.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(qux.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(qux.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(qux.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(qux.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(qux.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(qux.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(qux.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(qux.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(qux.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(qux.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(qux.data[1][1][2], 113.0);
}

BOOST_AUTO_TEST_CASE( interleaved_storage_order )
{
    suzerain::RealState<double,true> foo(2, 2, 2);

    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   1, &(foo.data[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   2, &(foo.data[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) + 2*2, &(foo.data[0][0][1]));

    BOOST_CHECK_EQUAL(   1, foo.data.strides()[0] );
    BOOST_CHECK_EQUAL(   2, foo.data.strides()[1] );
    BOOST_CHECK_EQUAL( 2*2, foo.data.strides()[2] );
}

BOOST_AUTO_TEST_CASE( non_interleaved_storage_order )
{
    suzerain::RealState<double,false> foo(2, 2, 2);

    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) + 2*2, &(foo.data[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   1, &(foo.data[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   2, &(foo.data[0][0][1]));

    BOOST_CHECK_EQUAL( 2*2, foo.data.strides()[0] );
    BOOST_CHECK_EQUAL(   1, foo.data.strides()[1] );
    BOOST_CHECK_EQUAL(   2, foo.data.strides()[2] );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( isConformant, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::RealState<double,Interleaved::value>  foo(2, 2, 2);
    suzerain::RealState<double,Interleaved::value>  bar(2, 2, 2);
    suzerain::RealState<double,Interleaved::value>  baz(1, 2, 2);
    suzerain::RealState<double,Interleaved::value>  qux(2, 1, 2);
    suzerain::RealState<double,Interleaved::value>  quux(2, 2, 1);
    suzerain::RealState<double,!Interleaved::value> quuux(2, 2, 2);

    BOOST_CHECK_EQUAL(true,  foo.isConformant(foo));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(bar));
    BOOST_CHECK_EQUAL(false, foo.isConformant(baz));
    BOOST_CHECK_EQUAL(false, foo.isConformant(qux));
    BOOST_CHECK_EQUAL(false, foo.isConformant(quux));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(quuux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( scale, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::RealState<double,Interleaved::value> foo(1, 2, 3);
    foo.data[0][0][0] =  2.0;
    foo.data[0][1][0] =  3.0;
    foo.data[0][0][1] =  5.0;
    foo.data[0][1][1] =  7.0;
    foo.data[0][0][2] = 11.0;
    foo.data[0][1][2] = 13.0;

    foo.scale(1.0);

    // Ensure foo.data is unchanged
    BOOST_CHECK_EQUAL(foo.data[0][0][0],  2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],  3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],  5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],  7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2], 11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2], 13.0);

    foo.scale(2.0);

    // Ensure foo.data was modified appropriately
    BOOST_CHECK_EQUAL(foo.data[0][0][0], 2.0* 2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0], 2.0* 3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1], 2.0* 5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1], 2.0* 7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2], 2.0*11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2], 2.0*13.0);

    foo.scale(0.0);

    // Ensure foo.data was zeroed
    BOOST_CHECK_EQUAL(foo.data[0][0][0], 0.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0], 0.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1], 0.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1], 0.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2], 0.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2], 0.0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( addScaled, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::RealState<double,Interleaved::value> foo(2, 2, 3);
    foo.data[0][0][0] =   2.0;
    foo.data[0][1][0] =   3.0;
    foo.data[0][0][1] =   5.0;
    foo.data[0][1][1] =   7.0;
    foo.data[0][0][2] =  11.0;
    foo.data[0][1][2] =  13.0;
    foo.data[1][0][0] = 102.0;
    foo.data[1][1][0] = 103.0;
    foo.data[1][0][1] = 105.0;
    foo.data[1][1][1] = 107.0;
    foo.data[1][0][2] = 111.0;
    foo.data[1][1][2] = 113.0;

    suzerain::RealState<double,Interleaved::value> bar(2, 2, 3);
    bar.data[0][0][0] =  17.0;
    bar.data[0][1][0] =  19.0;
    bar.data[0][0][1] =  23.0;
    bar.data[0][1][1] =  29.0;
    bar.data[0][0][2] =  31.0;
    bar.data[0][1][2] =  37.0;
    bar.data[1][0][0] = 117.0;
    bar.data[1][1][0] = 119.0;
    bar.data[1][0][1] = 123.0;
    bar.data[1][1][1] = 129.0;
    bar.data[1][0][2] = 131.0;
    bar.data[1][1][2] = 137.0;

    foo.addScaled(3.0, bar);

    // Ensure foo.data was modified
    BOOST_CHECK_EQUAL(foo.data[0][0][0],   2.0 + 3.0* 17.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],   3.0 + 3.0* 19.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],   5.0 + 3.0* 23.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],   7.0 + 3.0* 29.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2],  11.0 + 3.0* 31.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2],  13.0 + 3.0* 37.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][0], 102.0 + 3.0*117.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][0], 103.0 + 3.0*119.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][1], 105.0 + 3.0*123.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][1], 107.0 + 3.0*129.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][2], 111.0 + 3.0*131.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][2], 113.0 + 3.0*137.0);

    // Ensure we catch an operation between two nonconforming states
    suzerain::RealState<double,Interleaved::value> baz(2, 2, 2);
    BOOST_CHECK_THROW(foo.addScaled(3.0, baz), std::logic_error);

    // Check functionality against a different interleaving
    suzerain::RealState<double,!Interleaved::value> qux(2, 2, 3);
    BOOST_CHECK_NE(qux.interleaved(), bar.interleaved());
    qux.data[0][0][0] =   2.0;
    qux.data[0][1][0] =   3.0;
    qux.data[0][0][1] =   5.0;
    qux.data[0][1][1] =   7.0;
    qux.data[0][0][2] =  11.0;
    qux.data[0][1][2] =  13.0;
    qux.data[1][0][0] = 102.0;
    qux.data[1][1][0] = 103.0;
    qux.data[1][0][1] = 105.0;
    qux.data[1][1][1] = 107.0;
    qux.data[1][0][2] = 111.0;
    qux.data[1][1][2] = 113.0;

    bar.addScaled(4.0, qux);

    // Ensure bar.data was modified
    BOOST_CHECK_EQUAL(bar.data[0][0][0],  17.0 + 4.0*  2.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][0],  19.0 + 4.0*  3.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][1],  23.0 + 4.0*  5.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][1],  29.0 + 4.0*  7.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][2],  31.0 + 4.0* 11.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][2],  37.0 + 4.0* 13.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][0], 117.0 + 4.0*102.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][0], 119.0 + 4.0*103.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][1], 123.0 + 4.0*105.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][1], 129.0 + 4.0*107.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][2], 131.0 + 4.0*111.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][2], 137.0 + 4.0*113.0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(comparison_and_assignment,
                              Interleaved, bool_values)
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::RealState<double,Interleaved::value> foo(2, 2, 3);

    foo.data[0][0][0] =   2.0;
    foo.data[0][1][0] =   3.0;
    foo.data[0][0][1] =   5.0;
    foo.data[0][1][1] =   7.0;
    foo.data[0][0][2] =  11.0;
    foo.data[0][1][2] =  13.0;
    foo.data[1][0][0] = 102.0;
    foo.data[1][1][0] = 103.0;
    foo.data[1][0][1] = 105.0;
    foo.data[1][1][1] = 107.0;
    foo.data[1][0][2] = 111.0;
    foo.data[1][1][2] = 113.0;

    suzerain::RealState<double,Interleaved::value> bar(2, 2, 3);
    bar.data = foo.data;
    BOOST_CHECK_EQUAL(true,  bar.data == foo.data);
    BOOST_CHECK_EQUAL(true,  foo.data == bar.data);
    bar.data[0][1][0] = 555.0;
    BOOST_CHECK_EQUAL(false, bar.data == foo.data);
    BOOST_CHECK_EQUAL(false, foo.data == bar.data);

    BOOST_CHECK_EQUAL(foo.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][2], 113.0);

    BOOST_CHECK_EQUAL(bar.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][0], 555.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][2], 113.0);

    suzerain::RealState<double,!Interleaved::value> baz(2, 2, 3);
    BOOST_CHECK_NE(baz.interleaved(), bar.interleaved());
    baz.data = bar.data;
    BOOST_CHECK_EQUAL(true,  baz.data == bar.data);
    BOOST_CHECK_EQUAL(true,  bar.data == baz.data);
    bar.data[0][1][0] = 777.0;
    BOOST_CHECK_EQUAL(false, baz.data == bar.data);
    BOOST_CHECK_EQUAL(false, bar.data == baz.data);

    BOOST_CHECK_EQUAL(baz.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(baz.data[0][1][0], 555.0);
    BOOST_CHECK_EQUAL(baz.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(baz.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(baz.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(baz.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(baz.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(baz.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(baz.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(baz.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(baz.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(baz.data[1][1][2], 113.0);

    BOOST_CHECK_EQUAL(bar.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][0], 777.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][2], 113.0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::RealState<double,Interleaved::value> foo(2, 2, 3);

    foo.data[0][0][0] =   2.0;
    foo.data[0][1][0] =   3.0;
    foo.data[0][0][1] =   5.0;
    foo.data[0][1][1] =   7.0;
    foo.data[0][0][2] =  11.0;
    foo.data[0][1][2] =  13.0;
    foo.data[1][0][0] = 102.0;
    foo.data[1][1][0] = 103.0;
    foo.data[1][0][1] = 105.0;
    foo.data[1][1][1] = 107.0;
    foo.data[1][0][2] = 111.0;
    foo.data[1][1][2] = 113.0;

    suzerain::RealState<double,Interleaved::value> bar(2, 2, 3);

    bar.data[0][0][0] =  17.0;
    bar.data[0][1][0] =  19.0;
    bar.data[0][0][1] =  23.0;
    bar.data[0][1][1] =  29.0;
    bar.data[0][0][2] =  31.0;
    bar.data[0][1][2] =  37.0;
    bar.data[1][0][0] = 117.0;
    bar.data[1][1][0] = 119.0;
    bar.data[1][0][1] = 123.0;
    bar.data[1][1][1] = 129.0;
    bar.data[1][0][2] = 131.0;
    bar.data[1][1][2] = 137.0;

    foo.exchange(bar);

    BOOST_CHECK_EQUAL(foo.data[0][0][0],  17.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],  19.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],  23.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],  29.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2],  31.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2],  37.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][0], 117.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][0], 119.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][1], 123.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][1], 129.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][2], 131.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][2], 137.0);

    BOOST_CHECK_EQUAL(bar.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(bar.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(bar.data[1][1][2], 113.0);

    // Ensure we catch an operation between two nonconforming states
    suzerain::RealState<double,Interleaved::value> baz(2, 2, 2);
    BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);

    // Ensure we catch an operation between two different subclasses
    suzerain::ComplexState<double,Interleaved::value> qux(2,2,3);
    BOOST_CHECK_THROW(foo.exchange(qux), std::bad_cast);

    // Test exchange against different interleaving
    suzerain::RealState<double,!Interleaved::value> quux(2, 2, 3);
    BOOST_CHECK_NE(foo.interleaved(), quux.interleaved());

    quux.data[0][0][0] =   2.0;
    quux.data[0][1][0] =   3.0;
    quux.data[0][0][1] =   5.0;
    quux.data[0][1][1] =   7.0;
    quux.data[0][0][2] =  11.0;
    quux.data[0][1][2] =  13.0;
    quux.data[1][0][0] = 102.0;
    quux.data[1][1][0] = 103.0;
    quux.data[1][0][1] = 105.0;
    quux.data[1][1][1] = 107.0;
    quux.data[1][0][2] = 111.0;
    quux.data[1][1][2] = 113.0;

    foo.exchange(quux);

    BOOST_CHECK_EQUAL(quux.data[0][0][0],  17.0);
    BOOST_CHECK_EQUAL(quux.data[0][1][0],  19.0);
    BOOST_CHECK_EQUAL(quux.data[0][0][1],  23.0);
    BOOST_CHECK_EQUAL(quux.data[0][1][1],  29.0);
    BOOST_CHECK_EQUAL(quux.data[0][0][2],  31.0);
    BOOST_CHECK_EQUAL(quux.data[0][1][2],  37.0);
    BOOST_CHECK_EQUAL(quux.data[1][0][0], 117.0);
    BOOST_CHECK_EQUAL(quux.data[1][1][0], 119.0);
    BOOST_CHECK_EQUAL(quux.data[1][0][1], 123.0);
    BOOST_CHECK_EQUAL(quux.data[1][1][1], 129.0);
    BOOST_CHECK_EQUAL(quux.data[1][0][2], 131.0);
    BOOST_CHECK_EQUAL(quux.data[1][1][2], 137.0);

    BOOST_CHECK_EQUAL(foo.data[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][0], 102.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][0], 103.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][1], 105.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][1], 107.0);
    BOOST_CHECK_EQUAL(foo.data[1][0][2], 111.0);
    BOOST_CHECK_EQUAL(foo.data[1][1][2], 113.0);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( ComplexState )

// Ensures data, real, imag, and components views are self-consistent
template<typename FPT, bool Interleaved>
static
void checkSelfConsistent(const suzerain::ComplexState<FPT,Interleaved>& s) {
    for (int i = 0; i < s.variable_count; ++i) {
        for (int j = 0; j < s.vector_length; ++j) {
            for (int k = 0; k < s.vector_count; ++k) {
                BOOST_CHECK_EQUAL(s.data[i][j][k].real(), s.real[i][j][k]);
                BOOST_CHECK_EQUAL(s.real[i][j][k], s.components[0][i][j][k]);

                BOOST_CHECK_EQUAL(s.data[i][j][k].imag(), s.imag[i][j][k]);
                BOOST_CHECK_EQUAL(s.imag[i][j][k], s.components[1][i][j][k]);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( declare_pointer, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::ComplexState<double,Interleaved::value> *state;
    (void)state;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( constructor, Interleaved, bool_values )
{
    typedef std::complex<double> complex;
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    // Regular constructor
    suzerain::ComplexState<double,Interleaved::value> foo(2, 2, 3);
    BOOST_CHECK_EQUAL(foo.variable_count, 2);
    BOOST_CHECK_EQUAL(foo.vector_length, 2);
    BOOST_CHECK_EQUAL(foo.vector_count, 3);

    for (int i = 0; i < foo.variable_count; ++i) {
        foo.data[i][0][0] = complex(100*i + 2.0, -100*i - 2.0);
        foo.data[i][1][0] = complex(100*i + 3.0, -100*i - 3.0);
        foo.data[i][0][1] = complex(100*i + 5.0, -100*i - 5.0);
        foo.data[i][1][1] = complex(100*i + 7.0, -100*i - 7.0);
        foo.data[i][0][2] = complex(100*i +11.0, -100*i -11.0);
        foo.data[i][1][2] = complex(100*i +13.0, -100*i -13.0);
    }

    for (int i = 0; i < foo.variable_count; ++i) {
        BOOST_CHECK_EQUAL(foo.data[i][0][0].real(),  100*i +  2.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][0].imag(), -100*i -  2.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][0].real(),  100*i +  3.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][0].imag(), -100*i -  3.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][1].real(),  100*i +  5.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][1].imag(), -100*i -  5.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][1].real(),  100*i +  7.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][1].imag(), -100*i -  7.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][2].real(),  100*i + 11.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][2].imag(), -100*i - 11.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][2].real(),  100*i + 13.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][2].imag(), -100*i - 13.0);
    }

    checkSelfConsistent(foo);

    // Copy construct a second instance
    suzerain::ComplexState<double,Interleaved::value> bar(foo);
    BOOST_CHECK_EQUAL(bar.variable_count, 2);
    BOOST_CHECK_EQUAL(bar.vector_length, 2);
    BOOST_CHECK_EQUAL(bar.vector_count, 3);

    // Modify first instance's data
    for (int i = 0; i < foo.variable_count; ++i) {
        foo.data[i][0][0] += complex( 2.0, 555.0);
        foo.data[i][1][0] += complex( 3.0, 555.0);
        foo.data[i][0][1] += complex( 5.0, 555.0);
        foo.data[i][1][1] += complex( 7.0, 555.0);
        foo.data[i][0][2] += complex(11.0, 555.0);
        foo.data[i][1][2] += complex(13.0, 555.0);
    }

    // Ensure copy constructed data in second instance not modified
    for (int i = 0; i < bar.variable_count; ++i) {
        BOOST_CHECK_EQUAL(bar.data[i][0][0].real(),  100*i +  2.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][0].imag(), -100*i -  2.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][0].real(),  100*i +  3.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][0].imag(), -100*i -  3.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][1].real(),  100*i +  5.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][1].imag(), -100*i -  5.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][1].real(),  100*i +  7.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][1].imag(), -100*i -  7.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][2].real(),  100*i + 11.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][2].imag(), -100*i - 11.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][2].real(),  100*i + 13.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][2].imag(), -100*i - 13.0);
    }

    checkSelfConsistent(bar);

    // Copy construct a third instance with a different interleaving
    suzerain::ComplexState<double,!Interleaved::value> baz(bar);
    BOOST_CHECK_NE(bar.interleaved(), baz.interleaved());
    BOOST_CHECK_EQUAL(baz.variable_count, 2);
    BOOST_CHECK_EQUAL(baz.vector_length, 2);
    BOOST_CHECK_EQUAL(baz.vector_count, 3);

    // Modify second instance's data
    for (int i = 0; i < bar.variable_count; ++i) {
        for (int j = 0; j < bar.vector_length; ++j) {
            for (int k = 0; k < bar.vector_count; ++k) {
                bar.data[i][j][k] += complex((i+1)*(j+1), -1*(j+1)*(k+1));
            }
        }
    }

    // Ensure copy constructed data in third instance not modified
    for (int i = 0; i < baz.variable_count; ++i) {
        BOOST_CHECK_EQUAL(baz.data[i][0][0].real(),  100*i +  2.0);
        BOOST_CHECK_EQUAL(baz.data[i][0][0].imag(), -100*i -  2.0);
        BOOST_CHECK_EQUAL(baz.data[i][1][0].real(),  100*i +  3.0);
        BOOST_CHECK_EQUAL(baz.data[i][1][0].imag(), -100*i -  3.0);
        BOOST_CHECK_EQUAL(baz.data[i][0][1].real(),  100*i +  5.0);
        BOOST_CHECK_EQUAL(baz.data[i][0][1].imag(), -100*i -  5.0);
        BOOST_CHECK_EQUAL(baz.data[i][1][1].real(),  100*i +  7.0);
        BOOST_CHECK_EQUAL(baz.data[i][1][1].imag(), -100*i -  7.0);
        BOOST_CHECK_EQUAL(baz.data[i][0][2].real(),  100*i + 11.0);
        BOOST_CHECK_EQUAL(baz.data[i][0][2].imag(), -100*i - 11.0);
        BOOST_CHECK_EQUAL(baz.data[i][1][2].real(),  100*i + 13.0);
        BOOST_CHECK_EQUAL(baz.data[i][1][2].imag(), -100*i - 13.0);
    }

    checkSelfConsistent(baz);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, Interleaved, bool_values )
{
    typedef std::complex<double> complex;
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::ComplexState<double,Interleaved::value> foo(2, 2, 3), bar(2,2,3);
    for (int i = 0; i < foo.variable_count; ++i) {
        foo.data[i][0][0] = complex( 100*i + 2.0, -100*i - 2.0);
        foo.data[i][1][0] = complex( 100*i + 3.0, -100*i - 3.0);
        foo.data[i][0][1] = complex( 100*i + 5.0, -100*i - 5.0);
        foo.data[i][1][1] = complex( 100*i + 7.0, -100*i - 7.0);
        foo.data[i][0][2] = complex( 100*i +11.0, -100*i -11.0);
        foo.data[i][1][2] = complex( 100*i +13.0, -100*i -13.0);
    }

    // Self assignment
    foo = foo;

    for (int i = 0; i < foo.variable_count; ++i) {
        BOOST_CHECK_EQUAL(foo.data[i][0][0].real(),  100*i + 2.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][0].imag(), -100*i - 2.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][0].real(),  100*i + 3.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][0].imag(), -100*i - 3.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][1].real(),  100*i + 5.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][1].imag(), -100*i - 5.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][1].real(),  100*i + 7.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][1].imag(), -100*i - 7.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][2].real(),  100*i +11.0);
        BOOST_CHECK_EQUAL(foo.data[i][0][2].imag(), -100*i -11.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][2].real(),  100*i +13.0);
        BOOST_CHECK_EQUAL(foo.data[i][1][2].imag(), -100*i -13.0);
    }

    checkSelfConsistent(foo);

    // Normal assignment
    bar = foo;

    for (int i = 0; i < bar.variable_count; ++i) {
        BOOST_CHECK_EQUAL(bar.data[i][0][0].real(),  100*i + 2.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][0].imag(), -100*i - 2.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][0].real(),  100*i + 3.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][0].imag(), -100*i - 3.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][1].real(),  100*i + 5.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][1].imag(), -100*i - 5.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][1].real(),  100*i + 7.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][1].imag(), -100*i - 7.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][2].real(),  100*i +11.0);
        BOOST_CHECK_EQUAL(bar.data[i][0][2].imag(), -100*i -11.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][2].real(),  100*i +13.0);
        BOOST_CHECK_EQUAL(bar.data[i][1][2].imag(), -100*i -13.0);
    }

    checkSelfConsistent(bar);

    // Ensure we catch an operation between two nonconforming states
    suzerain::ComplexState<double,Interleaved::value> baz(2, 2, 2);
    BOOST_CHECK_THROW(baz = foo, std::logic_error);

    // Ensure we catch an operation between two different subclasses
    suzerain::RealState<double,Interleaved::value> qux(2,2,3);
    BOOST_CHECK_THROW(baz = qux, std::bad_cast);

    // Assignment from different interleaving
    suzerain::ComplexState<double,!Interleaved::value> quux(2, 2, 3);
    BOOST_CHECK_NE(quux.interleaved(), foo.interleaved());
    quux = bar;

    for (int i = 0; i < quux.variable_count; ++i) {
        BOOST_CHECK_EQUAL(quux.data[i][0][0].real(),  100*i + 2.0);
        BOOST_CHECK_EQUAL(quux.data[i][0][0].imag(), -100*i - 2.0);
        BOOST_CHECK_EQUAL(quux.data[i][1][0].real(),  100*i + 3.0);
        BOOST_CHECK_EQUAL(quux.data[i][1][0].imag(), -100*i - 3.0);
        BOOST_CHECK_EQUAL(quux.data[i][0][1].real(),  100*i + 5.0);
        BOOST_CHECK_EQUAL(quux.data[i][0][1].imag(), -100*i - 5.0);
        BOOST_CHECK_EQUAL(quux.data[i][1][1].real(),  100*i + 7.0);
        BOOST_CHECK_EQUAL(quux.data[i][1][1].imag(), -100*i - 7.0);
        BOOST_CHECK_EQUAL(quux.data[i][0][2].real(),  100*i +11.0);
        BOOST_CHECK_EQUAL(quux.data[i][0][2].imag(), -100*i -11.0);
        BOOST_CHECK_EQUAL(quux.data[i][1][2].real(),  100*i +13.0);
        BOOST_CHECK_EQUAL(quux.data[i][1][2].imag(), -100*i -13.0);
    }

    checkSelfConsistent(quux);
}

BOOST_AUTO_TEST_CASE( interleaved_storage_order )
{
    suzerain::ComplexState<double,true> foo(2, 2, 2);

    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   1, &(foo.data[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   2, &(foo.data[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) + 2*2, &(foo.data[0][0][1]));

    BOOST_CHECK_EQUAL(   1, foo.data.strides()[0] );
    BOOST_CHECK_EQUAL(   2, foo.data.strides()[1] );
    BOOST_CHECK_EQUAL( 2*2, foo.data.strides()[2] );

    BOOST_CHECK_EQUAL( &(foo.components[0][0][0][0]) +     1,
                       &(foo.components[1][0][0][0]));
    BOOST_CHECK_EQUAL( &(foo.components[0][0][0][0]) +     2,
                       &(foo.components[0][1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.components[0][0][0][0]) +   2*2,
                       &(foo.components[0][0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.components[0][0][0][0]) + 2*2*2,
                       &(foo.components[0][0][0][1]));

    BOOST_CHECK_EQUAL(     1, foo.components.strides()[0] );
    BOOST_CHECK_EQUAL(     2, foo.components.strides()[1] );
    BOOST_CHECK_EQUAL(   2*2, foo.components.strides()[2] );
    BOOST_CHECK_EQUAL( 2*2*2, foo.components.strides()[3] );

    BOOST_CHECK_EQUAL( &(foo.real[0][0][0]) +   1*2, &(foo.real[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.real[0][0][0]) +   2*2, &(foo.real[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.real[0][0][0]) + 2*2*2, &(foo.real[0][0][1]));

    BOOST_CHECK_EQUAL(   1*2, foo.real.strides()[0] );
    BOOST_CHECK_EQUAL(   2*2, foo.real.strides()[1] );
    BOOST_CHECK_EQUAL( 2*2*2, foo.real.strides()[2] );

    BOOST_CHECK_EQUAL( &(foo.imag[0][0][0]) +   1*2, &(foo.imag[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.imag[0][0][0]) +   2*2, &(foo.imag[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.imag[0][0][0]) + 2*2*2, &(foo.imag[0][0][1]));

    BOOST_CHECK_EQUAL(   1*2, foo.imag.strides()[0] );
    BOOST_CHECK_EQUAL(   2*2, foo.imag.strides()[1] );
    BOOST_CHECK_EQUAL( 2*2*2, foo.imag.strides()[2] );
}

BOOST_AUTO_TEST_CASE( non_interleaved_storage_order )
{
    suzerain::ComplexState<double,false> foo(2, 2, 2);

    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) + 2*2, &(foo.data[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   1, &(foo.data[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   2, &(foo.data[0][0][1]));

    BOOST_CHECK_EQUAL( 2*2, foo.data.strides()[0] );
    BOOST_CHECK_EQUAL(   1, foo.data.strides()[1] );
    BOOST_CHECK_EQUAL(   2, foo.data.strides()[2] );

    BOOST_CHECK_EQUAL( &(foo.components[0][0][0][0]) +     1,
                       &(foo.components[1][0][0][0]));
    BOOST_CHECK_EQUAL( &(foo.components[0][0][0][0]) + 2*2*2,
                       &(foo.components[0][1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.components[0][0][0][0]) +     2,
                       &(foo.components[0][0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.components[0][0][0][0]) +   2*2,
                       &(foo.components[0][0][0][1]));

    BOOST_CHECK_EQUAL(     1, foo.components.strides()[0] );
    BOOST_CHECK_EQUAL( 2*2*2, foo.components.strides()[1] );
    BOOST_CHECK_EQUAL(     2, foo.components.strides()[2] );
    BOOST_CHECK_EQUAL(   2*2, foo.components.strides()[3] );

    BOOST_CHECK_EQUAL( &(foo.real[0][0][0]) + 2*2*2, &(foo.real[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.real[0][0][0]) +   1*2, &(foo.real[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.real[0][0][0]) +   2*2, &(foo.real[0][0][1]));

    BOOST_CHECK_EQUAL( 2*2*2, foo.real.strides()[0] );
    BOOST_CHECK_EQUAL(   1*2, foo.real.strides()[1] );
    BOOST_CHECK_EQUAL(   2*2, foo.real.strides()[2] );

    BOOST_CHECK_EQUAL( &(foo.imag[0][0][0]) + 2*2*2, &(foo.imag[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.imag[0][0][0]) +   1*2, &(foo.imag[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.imag[0][0][0]) +   2*2, &(foo.imag[0][0][1]));

    BOOST_CHECK_EQUAL( 2*2*2, foo.imag.strides()[0] );
    BOOST_CHECK_EQUAL(   1*2, foo.imag.strides()[1] );
    BOOST_CHECK_EQUAL(   2*2, foo.imag.strides()[2] );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( isConformant, Interleaved, bool_values )
{
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::ComplexState<double,Interleaved::value>  foo(2, 2, 2);
    suzerain::ComplexState<double,Interleaved::value>  bar(2, 2, 2);
    suzerain::ComplexState<double,Interleaved::value>  baz(1, 2, 2);
    suzerain::ComplexState<double,Interleaved::value>  qux(2, 1, 2);
    suzerain::ComplexState<double,Interleaved::value>  quux(2, 2, 1);
    suzerain::ComplexState<double,!Interleaved::value> quuux(2, 2, 2);

    BOOST_CHECK_EQUAL(true,  foo.isConformant(foo));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(bar));
    BOOST_CHECK_EQUAL(false, foo.isConformant(baz));
    BOOST_CHECK_EQUAL(false, foo.isConformant(qux));
    BOOST_CHECK_EQUAL(false, foo.isConformant(quux));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(quuux));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( scale, Interleaved, bool_values )
{
    typedef std::complex<double> complex;
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::ComplexState<double,Interleaved::value> foo(1, 2, 3);
    foo.data[0][0][0] = complex( 2.0, - 2.0);
    foo.data[0][1][0] = complex( 3.0, - 3.0);
    foo.data[0][0][1] = complex( 5.0, - 5.0);
    foo.data[0][1][1] = complex( 7.0, - 7.0);
    foo.data[0][0][2] = complex(11.0, -11.0);
    foo.data[0][1][2] = complex(13.0, -13.0);

    foo.scale(1.0);

    // Ensure foo.data is unchanged
    BOOST_CHECK_EQUAL(foo.real[0][0][0],   2.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][0],   3.0);
    BOOST_CHECK_EQUAL(foo.real[0][0][1],   5.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][1],   7.0);
    BOOST_CHECK_EQUAL(foo.real[0][0][2],  11.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][2],  13.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][0], - 2.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][0], - 3.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][1], - 5.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][1], - 7.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][2], -11.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][2], -13.0);

    checkSelfConsistent(foo);

    foo.scale(2.0);

    // Ensure foo.data was modified appropriately
    BOOST_CHECK_EQUAL(foo.real[0][0][0],  2.0* 2.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][0],  2.0* 3.0);
    BOOST_CHECK_EQUAL(foo.real[0][0][1],  2.0* 5.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][1],  2.0* 7.0);
    BOOST_CHECK_EQUAL(foo.real[0][0][2],  2.0*11.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][2],  2.0*13.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][0], -2.0* 2.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][0], -2.0* 3.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][1], -2.0* 5.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][1], -2.0* 7.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][2], -2.0*11.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][2], -2.0*13.0);

    checkSelfConsistent(foo);

    foo.scale(0.0);

    // Ensure foo.data was zeroed
    BOOST_CHECK_EQUAL(foo.real[0][0][0], 0.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][0], 0.0);
    BOOST_CHECK_EQUAL(foo.real[0][0][1], 0.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][1], 0.0);
    BOOST_CHECK_EQUAL(foo.real[0][0][2], 0.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][2], 0.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][0], 0.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][0], 0.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][1], 0.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][1], 0.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][2], 0.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][2], 0.0);

    checkSelfConsistent(foo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( addScaled, Interleaved, bool_values )
{
    typedef std::complex<double> complex;
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    suzerain::ComplexState<double,Interleaved::value> foo(2, 2, 3);
    for (int i = 0; i < foo.variable_count; ++i) {
        foo.data[i][0][0] = complex( 100*i + 2.0, -100*i - 2.0);
        foo.data[i][1][0] = complex( 100*i + 3.0, -100*i - 3.0);
        foo.data[i][0][1] = complex( 100*i + 5.0, -100*i - 5.0);
        foo.data[i][1][1] = complex( 100*i + 7.0, -100*i - 7.0);
        foo.data[i][0][2] = complex( 100*i +11.0, -100*i -11.0);
        foo.data[i][1][2] = complex( 100*i +13.0, -100*i -13.0);
    }

    suzerain::ComplexState<double,Interleaved::value> bar(2, 2, 3);
    for (int i = 0; i < foo.variable_count; ++i) {
        bar.data[i][0][0] = complex( 100*i +17.0, -100*i -17.0);
        bar.data[i][1][0] = complex( 100*i +19.0, -100*i -19.0);
        bar.data[i][0][1] = complex( 100*i +23.0, -100*i -23.0);
        bar.data[i][1][1] = complex( 100*i +29.0, -100*i -29.0);
        bar.data[i][0][2] = complex( 100*i +31.0, -100*i -31.0);
        bar.data[i][1][2] = complex( 100*i +37.0, -100*i -37.0);
    }

    foo.addScaled(3.0, bar);

    // Ensure foo.data was modified
    for (int i = 0; i < foo.variable_count; ++i) {
        BOOST_CHECK_EQUAL(foo.real[i][0][0], (100*i+ 2.0) + 3.0*(100*i+17.0));
        BOOST_CHECK_EQUAL(foo.real[i][1][0], (100*i+ 3.0) + 3.0*(100*i+19.0));
        BOOST_CHECK_EQUAL(foo.real[i][0][1], (100*i+ 5.0) + 3.0*(100*i+23.0));
        BOOST_CHECK_EQUAL(foo.real[i][1][1], (100*i+ 7.0) + 3.0*(100*i+29.0));
        BOOST_CHECK_EQUAL(foo.real[i][0][2], (100*i+11.0) + 3.0*(100*i+31.0));
        BOOST_CHECK_EQUAL(foo.real[i][1][2], (100*i+13.0) + 3.0*(100*i+37.0));

        BOOST_CHECK_EQUAL(foo.imag[i][0][0], (-100*i- 2.0) + 3.0*(-100*i-17.0));
        BOOST_CHECK_EQUAL(foo.imag[i][1][0], (-100*i- 3.0) + 3.0*(-100*i-19.0));
        BOOST_CHECK_EQUAL(foo.imag[i][0][1], (-100*i- 5.0) + 3.0*(-100*i-23.0));
        BOOST_CHECK_EQUAL(foo.imag[i][1][1], (-100*i- 7.0) + 3.0*(-100*i-29.0));
        BOOST_CHECK_EQUAL(foo.imag[i][0][2], (-100*i-11.0) + 3.0*(-100*i-31.0));
        BOOST_CHECK_EQUAL(foo.imag[i][1][2], (-100*i-13.0) + 3.0*(-100*i-37.0));
    }

    checkSelfConsistent(foo);

    // Ensure we catch an operation between two nonconforming states
    suzerain::ComplexState<double,Interleaved::value> baz(2, 2, 2);
    BOOST_CHECK_THROW(foo.addScaled(3.0, baz), std::logic_error);

    // Ensure correct behavior on a different interleaving
    suzerain::ComplexState<double,!Interleaved::value> qux(2, 2, 3);
    BOOST_CHECK_NE(qux.interleaved(), bar.interleaved());
    for (int i = 0; i < foo.variable_count; ++i) {
        qux.data[i][0][0] = complex( 100*i + 2.0, -100*i - 2.0);
        qux.data[i][1][0] = complex( 100*i + 3.0, -100*i - 3.0);
        qux.data[i][0][1] = complex( 100*i + 5.0, -100*i - 5.0);
        qux.data[i][1][1] = complex( 100*i + 7.0, -100*i - 7.0);
        qux.data[i][0][2] = complex( 100*i +11.0, -100*i -11.0);
        qux.data[i][1][2] = complex( 100*i +13.0, -100*i -13.0);
    }

    qux.addScaled(3.0, bar);

    for (int i = 0; i < qux.variable_count; ++i) {
        BOOST_CHECK_EQUAL(qux.real[i][0][0], (100*i+ 2.0) + 3.0*(100*i+17.0));
        BOOST_CHECK_EQUAL(qux.real[i][1][0], (100*i+ 3.0) + 3.0*(100*i+19.0));
        BOOST_CHECK_EQUAL(qux.real[i][0][1], (100*i+ 5.0) + 3.0*(100*i+23.0));
        BOOST_CHECK_EQUAL(qux.real[i][1][1], (100*i+ 7.0) + 3.0*(100*i+29.0));
        BOOST_CHECK_EQUAL(qux.real[i][0][2], (100*i+11.0) + 3.0*(100*i+31.0));
        BOOST_CHECK_EQUAL(qux.real[i][1][2], (100*i+13.0) + 3.0*(100*i+37.0));

        BOOST_CHECK_EQUAL(qux.imag[i][0][0], (-100*i- 2.0) + 3.0*(-100*i-17.0));
        BOOST_CHECK_EQUAL(qux.imag[i][1][0], (-100*i- 3.0) + 3.0*(-100*i-19.0));
        BOOST_CHECK_EQUAL(qux.imag[i][0][1], (-100*i- 5.0) + 3.0*(-100*i-23.0));
        BOOST_CHECK_EQUAL(qux.imag[i][1][1], (-100*i- 7.0) + 3.0*(-100*i-29.0));
        BOOST_CHECK_EQUAL(qux.imag[i][0][2], (-100*i-11.0) + 3.0*(-100*i-31.0));
        BOOST_CHECK_EQUAL(qux.imag[i][1][2], (-100*i-13.0) + 3.0*(-100*i-37.0));
    }

    checkSelfConsistent(qux);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(comparison_and_assignment,
                              Interleaved, bool_values )
{
    typedef std::complex<double> complex;
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    // Check comparison and assignment to a same interleaving
    suzerain::ComplexState<double,Interleaved::value> foo(2, 2, 3);
    for (int i = 0; i < foo.variable_count; ++i) {
        foo.data[i][0][0] = complex( 100*i + 2.0, -100*i - 2.0);
        foo.data[i][1][0] = complex( 100*i + 3.0, -100*i - 3.0);
        foo.data[i][0][1] = complex( 100*i + 5.0, -100*i - 5.0);
        foo.data[i][1][1] = complex( 100*i + 7.0, -100*i - 7.0);
        foo.data[i][0][2] = complex( 100*i +11.0, -100*i -11.0);
        foo.data[i][1][2] = complex( 100*i +13.0, -100*i -13.0);
    }

    suzerain::ComplexState<double,Interleaved::value> bar(2, 2, 3);
    BOOST_CHECK_EQUAL(foo.interleaved(), bar.interleaved());
    bar.data = foo.data;
    BOOST_CHECK_EQUAL(true, bar.data == foo.data);
    BOOST_CHECK_EQUAL(true, foo.data == bar.data);
    for (int i = 0; i < bar.variable_count; ++i) {
        bar.data[i][1][0] = complex(100*i + 555.0, -100*i - 555.0);
        bar.data[i][1][1] = 100*i + 777.0;
    }
    BOOST_CHECK_EQUAL(false, bar.data == foo.data);
    BOOST_CHECK_EQUAL(false, foo.data == bar.data);

    for (int i = 0; i < foo.variable_count; ++i) {
        BOOST_CHECK_EQUAL(foo.data[i][0][0],complex(100*i+ 2.0,-100*i- 2.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][0],complex(100*i+ 3.0,-100*i- 3.0));
        BOOST_CHECK_EQUAL(foo.data[i][0][1],complex(100*i+ 5.0,-100*i- 5.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][1],complex(100*i+ 7.0,-100*i- 7.0));
        BOOST_CHECK_EQUAL(foo.data[i][0][2],complex(100*i+11.0,-100*i-11.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][2],complex(100*i+13.0,-100*i-13.0));
    }

    checkSelfConsistent(foo);

    for (int i = 0; i < bar.variable_count; ++i) {
        BOOST_CHECK_EQUAL(bar.data[i][0][0],complex(100*i+  2.0,-100*i-  2.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][0],complex(100*i+555.0,-100*i-555.0));
        BOOST_CHECK_EQUAL(bar.data[i][0][1],complex(100*i+  5.0,-100*i-  5.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][1],complex(100*i+777.0,         0.0));
        BOOST_CHECK_EQUAL(bar.data[i][0][2],complex(100*i+ 11.0,-100*i- 11.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][2],complex(100*i+ 13.0,-100*i- 13.0));
    }

    checkSelfConsistent(bar);

    // Check comparison and assignment to a different interleaving
    suzerain::ComplexState<double,!Interleaved::value> qux(2, 2, 3);
    BOOST_CHECK_NE(qux.interleaved(), bar.interleaved());
    qux.data = bar.data;
    BOOST_CHECK_EQUAL(true, qux.data == bar.data);
    BOOST_CHECK_EQUAL(true, bar.data == qux.data);
    for (int i = 0; i < bar.variable_count; ++i) {
        qux.data[i][1][0] = complex( 100*i + 3.0, -100*i - 3.0);
        qux.data[i][1][1] = complex( 100*i + 7.0, -100*i - 7.0);
    }
    BOOST_CHECK_EQUAL(false, qux.data == bar.data);
    BOOST_CHECK_EQUAL(false, bar.data == qux.data);

    for (int i = 0; i < qux.variable_count; ++i) {
        BOOST_CHECK_EQUAL(qux.data[i][0][0],complex(100*i+ 2.0,-100*i- 2.0));
        BOOST_CHECK_EQUAL(qux.data[i][1][0],complex(100*i+ 3.0,-100*i- 3.0));
        BOOST_CHECK_EQUAL(qux.data[i][0][1],complex(100*i+ 5.0,-100*i- 5.0));
        BOOST_CHECK_EQUAL(qux.data[i][1][1],complex(100*i+ 7.0,-100*i- 7.0));
        BOOST_CHECK_EQUAL(qux.data[i][0][2],complex(100*i+11.0,-100*i-11.0));
        BOOST_CHECK_EQUAL(qux.data[i][1][2],complex(100*i+13.0,-100*i-13.0));
    }

    checkSelfConsistent(qux);

    for (int i = 0; i < bar.variable_count; ++i) {
        BOOST_CHECK_EQUAL(bar.data[i][0][0],complex(100*i+  2.0,-100*i-  2.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][0],complex(100*i+555.0,-100*i-555.0));
        BOOST_CHECK_EQUAL(bar.data[i][0][1],complex(100*i+  5.0,-100*i-  5.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][1],complex(100*i+777.0,         0.0));
        BOOST_CHECK_EQUAL(bar.data[i][0][2],complex(100*i+ 11.0,-100*i- 11.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][2],complex(100*i+ 13.0,-100*i- 13.0));
    }

    checkSelfConsistent(bar);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( exchange, Interleaved, bool_values )
{
    typedef std::complex<double> complex;
    BOOST_TEST_MESSAGE("Interleaved = " << !!Interleaved::value);

    // Exchange with same interleaving
    suzerain::ComplexState<double,Interleaved::value> foo(2, 2, 3);
    for (int i = 0; i < foo.variable_count; ++i) {
        foo.data[i][0][0] = complex(100*i + 2.0, -100*i -  2.0);
        foo.data[i][1][0] = complex(100*i + 3.0, -100*i -  3.0);
        foo.data[i][0][1] = complex(100*i + 5.0, -100*i -  5.0);
        foo.data[i][1][1] = complex(100*i + 7.0, -100*i -  7.0);
        foo.data[i][0][2] = complex(100*i +11.0, -100*i - 11.0);
        foo.data[i][1][2] = complex(100*i +13.0, -100*i - 13.0);
    }

    suzerain::ComplexState<double,Interleaved::value> bar(2, 2, 3);
    for (int i = 0; i < bar.variable_count; ++i) {
        bar.data[i][0][0] = complex(100*i + 17.0, -100*i - 17.0);
        bar.data[i][1][0] = complex(100*i + 19.0, -100*i - 19.0);
        bar.data[i][0][1] = complex(100*i + 23.0, -100*i - 23.0);
        bar.data[i][1][1] = complex(100*i + 29.0, -100*i - 29.0);
        bar.data[i][0][2] = complex(100*i + 31.0, -100*i - 31.0);
        bar.data[i][1][2] = complex(100*i + 37.0, -100*i - 37.0);
    }

    foo.exchange(bar);

    for (int i = 0; i < foo.variable_count; ++i) {
        BOOST_CHECK_EQUAL(foo.data[i][0][0], complex(100*i+17.0, -100*i-17.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][0], complex(100*i+19.0, -100*i-19.0));
        BOOST_CHECK_EQUAL(foo.data[i][0][1], complex(100*i+23.0, -100*i-23.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][1], complex(100*i+29.0, -100*i-29.0));
        BOOST_CHECK_EQUAL(foo.data[i][0][2], complex(100*i+31.0, -100*i-31.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][2], complex(100*i+37.0, -100*i-37.0));
    }

    checkSelfConsistent(foo);

    for (int i = 0; i < bar.variable_count; ++i) {
        BOOST_CHECK_EQUAL(bar.data[i][0][0], complex(100*i+ 2.0, -100*i- 2.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][0], complex(100*i+ 3.0, -100*i- 3.0));
        BOOST_CHECK_EQUAL(bar.data[i][0][1], complex(100*i+ 5.0, -100*i- 5.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][1], complex(100*i+ 7.0, -100*i- 7.0));
        BOOST_CHECK_EQUAL(bar.data[i][0][2], complex(100*i+11.0, -100*i-11.0));
        BOOST_CHECK_EQUAL(bar.data[i][1][2], complex(100*i+13.0, -100*i-13.0));
    }

    checkSelfConsistent(bar);

    // Ensure we catch an operation between two nonconforming states
    suzerain::ComplexState<double,Interleaved::value> baz(2, 2, 2);
    BOOST_CHECK_THROW(foo.exchange(baz), std::logic_error);

    // Ensure we catch an operation between two different subclasses
    suzerain::RealState<double,Interleaved::value> qux(2,2,3);
    BOOST_CHECK_THROW(foo.exchange(qux), std::bad_cast);

    // Exchange with different interleaving
    suzerain::ComplexState<double,!Interleaved::value> quux(2, 2, 3);
    BOOST_CHECK_NE(quux.interleaved(), foo.interleaved());
    for (int i = 0; i < quux.variable_count; ++i) {
        quux.data[i][0][0] = complex(100*i + 2.0, -100*i -  2.0);
        quux.data[i][1][0] = complex(100*i + 3.0, -100*i -  3.0);
        quux.data[i][0][1] = complex(100*i + 5.0, -100*i -  5.0);
        quux.data[i][1][1] = complex(100*i + 7.0, -100*i -  7.0);
        quux.data[i][0][2] = complex(100*i +11.0, -100*i - 11.0);
        quux.data[i][1][2] = complex(100*i +13.0, -100*i - 13.0);
    }

    foo.exchange(quux);

    for (int i = 0; i < foo.variable_count; ++i) {
        BOOST_CHECK_EQUAL(foo.data[i][0][0], complex(100*i+ 2.0, -100*i- 2.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][0], complex(100*i+ 3.0, -100*i- 3.0));
        BOOST_CHECK_EQUAL(foo.data[i][0][1], complex(100*i+ 5.0, -100*i- 5.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][1], complex(100*i+ 7.0, -100*i- 7.0));
        BOOST_CHECK_EQUAL(foo.data[i][0][2], complex(100*i+11.0, -100*i-11.0));
        BOOST_CHECK_EQUAL(foo.data[i][1][2], complex(100*i+13.0, -100*i-13.0));
    }

    checkSelfConsistent(foo);

    for (int i = 0; i < quux.variable_count; ++i) {
        BOOST_CHECK_EQUAL(quux.data[i][0][0],complex(100*i+17.0, -100*i-17.0));
        BOOST_CHECK_EQUAL(quux.data[i][1][0],complex(100*i+19.0, -100*i-19.0));
        BOOST_CHECK_EQUAL(quux.data[i][0][1],complex(100*i+23.0, -100*i-23.0));
        BOOST_CHECK_EQUAL(quux.data[i][1][1],complex(100*i+29.0, -100*i-29.0));
        BOOST_CHECK_EQUAL(quux.data[i][0][2],complex(100*i+31.0, -100*i-31.0));
        BOOST_CHECK_EQUAL(quux.data[i][1][2],complex(100*i+37.0, -100*i-37.0));
    }

    checkSelfConsistent(quux);
}

BOOST_AUTO_TEST_SUITE_END()
