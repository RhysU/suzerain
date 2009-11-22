#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/state.hpp>
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( RealState )

BOOST_AUTO_TEST_CASE( declare_pointer )
{
    suzerain::RealState<double> *state = NULL;
}

BOOST_AUTO_TEST_CASE( constructor )
{
    suzerain::RealState<double> foo(1, 2, 3);
    BOOST_CHECK_EQUAL(foo.variable_count, 1);
    BOOST_CHECK_EQUAL(foo.vector_length, 2);
    BOOST_CHECK_EQUAL(foo.vector_count, 3);

    foo.data[0][0][0] =  2.0;
    foo.data[0][1][0] =  3.0;
    foo.data[0][0][1] =  5.0;
    foo.data[0][1][1] =  7.0;
    foo.data[0][0][2] = 11.0;
    foo.data[0][1][2] = 13.0;

    BOOST_CHECK_EQUAL(foo.data[0][0][0],  2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],  3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],  5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],  7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2], 11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2], 13.0);
}

BOOST_AUTO_TEST_CASE( fortran_storage_order )
{
    suzerain::RealState<double> foo(2, 2, 2);

    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   1, &(foo.data[1][0][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) +   2, &(foo.data[0][1][0]));
    BOOST_CHECK_EQUAL( &(foo.data[0][0][0]) + 2*2, &(foo.data[0][0][1]));

    BOOST_CHECK_EQUAL(   1, foo.data.strides()[0] );
    BOOST_CHECK_EQUAL(   2, foo.data.strides()[1] );
    BOOST_CHECK_EQUAL( 2*2, foo.data.strides()[2] );
}

BOOST_AUTO_TEST_CASE( isConformant )
{
    suzerain::RealState<double> foo(2, 2, 2);
    suzerain::RealState<double> bar(2, 2, 2);
    suzerain::RealState<double> baz(1, 2, 2);
    suzerain::RealState<double> qux(2, 1, 2);
    suzerain::RealState<double> quux(2, 2, 1);

    BOOST_CHECK_EQUAL(true,  foo.isConformant(&foo));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(&bar));
    BOOST_CHECK_EQUAL(false, foo.isConformant(&baz));
    BOOST_CHECK_EQUAL(false, foo.isConformant(&qux));
    BOOST_CHECK_EQUAL(false, foo.isConformant(&quux));
}

BOOST_AUTO_TEST_CASE( scaleAddScaled )
{
    suzerain::RealState<double> foo(1, 2, 3);
    foo.data[0][0][0] =  2.0;
    foo.data[0][1][0] =  3.0;
    foo.data[0][0][1] =  5.0;
    foo.data[0][1][1] =  7.0;
    foo.data[0][0][2] = 11.0;
    foo.data[0][1][2] = 13.0;

    suzerain::RealState<double> bar(1, 2, 3);
    bar.data[0][0][0] = 17.0;
    bar.data[0][1][0] = 19.0;
    bar.data[0][0][1] = 23.0;
    bar.data[0][1][1] = 29.0;
    bar.data[0][0][2] = 31.0;
    bar.data[0][1][2] = 37.0;

    foo.scaleAddScaled(2.0, 3.0, &bar);

    // Ensure foo.data was modified
    BOOST_CHECK_EQUAL(foo.data[0][0][0], 2.0* 2.0 + 3.0*17.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0], 2.0* 3.0 + 3.0*19.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1], 2.0* 5.0 + 3.0*23.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1], 2.0* 7.0 + 3.0*29.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2], 2.0*11.0 + 3.0*31.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2], 2.0*13.0 + 3.0*37.0);

    // Ensure bar.data was not modified
    BOOST_CHECK_EQUAL(bar.data[0][0][0], 17.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][0], 19.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][1], 23.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][1], 29.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][2], 31.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][2], 37.0);

    // Ensure we catch an operation between two nonconforming states
    suzerain::RealState<double> baz(2, 2, 2);
    BOOST_CHECK_THROW(foo.scaleAddScaled(2.0, 3.0, &baz),
                      suzerain::logic_error);
}

BOOST_AUTO_TEST_CASE( comparison_and_assignment )
{
    suzerain::RealState<double> foo(1, 2, 3);

    foo.data[0][0][0] =  2.0;
    foo.data[0][1][0] =  3.0;
    foo.data[0][0][1] =  5.0;
    foo.data[0][1][1] =  7.0;
    foo.data[0][0][2] = 11.0;
    foo.data[0][1][2] = 13.0;

    suzerain::RealState<double> bar(1, 2, 3);
    bar.data = foo.data;

    BOOST_CHECK_EQUAL(true, bar.data == foo.data);

    bar.data[0][1][0] = 555.0;

    BOOST_CHECK_EQUAL(false, bar.data == foo.data);

    BOOST_CHECK_EQUAL(foo.data[0][0][0],  2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0],  3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1],  5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1],  7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2], 11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2], 13.0);

    BOOST_CHECK_EQUAL(bar.data[0][0][0],  2.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][0],  555.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][1],  5.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][1],  7.0);
    BOOST_CHECK_EQUAL(bar.data[0][0][2], 11.0);
    BOOST_CHECK_EQUAL(bar.data[0][1][2], 13.0);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( ComplexState )

BOOST_AUTO_TEST_CASE( declare_pointer )
{
    suzerain::ComplexState<double> *state = NULL;
}

BOOST_AUTO_TEST_CASE( constructor )
{
    suzerain::ComplexState<double> foo(1, 2, 3);

    BOOST_CHECK_EQUAL(foo.variable_count, 1);
    BOOST_CHECK_EQUAL(foo.vector_length, 2);
    BOOST_CHECK_EQUAL(foo.vector_count, 3);

    typedef std::complex<double> complex;
    foo.data[0][0][0] = complex( 2.0, - 2.0);
    foo.data[0][1][0] = complex( 3.0, - 3.0);
    foo.data[0][0][1] = complex( 5.0, - 5.0);
    foo.data[0][1][1] = complex( 7.0, - 7.0);
    foo.data[0][0][2] = complex(11.0, -11.0);
    foo.data[0][1][2] = complex(13.0, -13.0);

    BOOST_CHECK_EQUAL(foo.data[0][0][0].real(),   2.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][0].imag(), - 2.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0].real(),   3.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][0].imag(), - 3.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1].real(),   5.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][1].imag(), - 5.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1].real(),   7.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][1].imag(), - 7.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2].real(),  11.0);
    BOOST_CHECK_EQUAL(foo.data[0][0][2].imag(), -11.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2].real(),  13.0);
    BOOST_CHECK_EQUAL(foo.data[0][1][2].imag(), -13.0);

    BOOST_CHECK_EQUAL(foo.components[0][0][0][0],   2.0);
    BOOST_CHECK_EQUAL(foo.components[1][0][0][0], - 2.0);
    BOOST_CHECK_EQUAL(foo.components[0][0][1][0],   3.0);
    BOOST_CHECK_EQUAL(foo.components[1][0][1][0], - 3.0);
    BOOST_CHECK_EQUAL(foo.components[0][0][0][1],   5.0);
    BOOST_CHECK_EQUAL(foo.components[1][0][0][1], - 5.0);
    BOOST_CHECK_EQUAL(foo.components[0][0][1][1],   7.0);
    BOOST_CHECK_EQUAL(foo.components[1][0][1][1], - 7.0);
    BOOST_CHECK_EQUAL(foo.components[0][0][0][2],  11.0);
    BOOST_CHECK_EQUAL(foo.components[1][0][0][2], -11.0);
    BOOST_CHECK_EQUAL(foo.components[0][0][1][2],  13.0);
    BOOST_CHECK_EQUAL(foo.components[1][0][1][2], -13.0);

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
}

BOOST_AUTO_TEST_CASE( fortran_storage_order )
{
    suzerain::ComplexState<double> foo(2, 2, 2);

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

BOOST_AUTO_TEST_CASE( isConformant )
{
    suzerain::ComplexState<double> foo(2, 2, 2);
    suzerain::ComplexState<double> bar(2, 2, 2);
    suzerain::ComplexState<double> baz(1, 2, 2);
    suzerain::ComplexState<double> qux(2, 1, 2);
    suzerain::ComplexState<double> quux(2, 2, 1);

    BOOST_CHECK_EQUAL(true,  foo.isConformant(&foo));
    BOOST_CHECK_EQUAL(true,  foo.isConformant(&bar));
    BOOST_CHECK_EQUAL(false, foo.isConformant(&baz));
    BOOST_CHECK_EQUAL(false, foo.isConformant(&qux));
    BOOST_CHECK_EQUAL(false, foo.isConformant(&quux));
}

BOOST_AUTO_TEST_CASE( scaleAddScaled )
{
    suzerain::ComplexState<double> foo(1, 2, 3);
    typedef std::complex<double> complex;
    foo.data[0][0][0] = complex( 2.0, - 2.0);
    foo.data[0][1][0] = complex( 3.0, - 3.0);
    foo.data[0][0][1] = complex( 5.0, - 5.0);
    foo.data[0][1][1] = complex( 7.0, - 7.0);
    foo.data[0][0][2] = complex(11.0, -11.0);
    foo.data[0][1][2] = complex(13.0, -13.0);

    suzerain::ComplexState<double> bar(1, 2, 3);
    bar.data[0][0][0] = complex(17.0, -17.0);
    bar.data[0][1][0] = complex(19.0, -19.0);
    bar.data[0][0][1] = complex(23.0, -23.0);
    bar.data[0][1][1] = complex(29.0, -29.0);
    bar.data[0][0][2] = complex(31.0, -31.0);
    bar.data[0][1][2] = complex(37.0, -37.0);

    foo.scaleAddScaled(2.0, 3.0, &bar);

    // Ensure foo.data was modified
    BOOST_CHECK_EQUAL(foo.real[0][0][0],  2.0* 2.0 + 3.0*17.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][0],  2.0* 3.0 + 3.0*19.0);
    BOOST_CHECK_EQUAL(foo.real[0][0][1],  2.0* 5.0 + 3.0*23.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][1],  2.0* 7.0 + 3.0*29.0);
    BOOST_CHECK_EQUAL(foo.real[0][0][2],  2.0*11.0 + 3.0*31.0);
    BOOST_CHECK_EQUAL(foo.real[0][1][2],  2.0*13.0 + 3.0*37.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][0], -2.0* 2.0 - 3.0*17.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][0], -2.0* 3.0 - 3.0*19.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][1], -2.0* 5.0 - 3.0*23.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][1], -2.0* 7.0 - 3.0*29.0);
    BOOST_CHECK_EQUAL(foo.imag[0][0][2], -2.0*11.0 - 3.0*31.0);
    BOOST_CHECK_EQUAL(foo.imag[0][1][2], -2.0*13.0 - 3.0*37.0);

    // Ensure bar.data was not modified
    BOOST_CHECK_EQUAL(bar.data[0][0][0], complex(17.0, -17.0));
    BOOST_CHECK_EQUAL(bar.data[0][1][0], complex(19.0, -19.0));
    BOOST_CHECK_EQUAL(bar.data[0][0][1], complex(23.0, -23.0));
    BOOST_CHECK_EQUAL(bar.data[0][1][1], complex(29.0, -29.0));
    BOOST_CHECK_EQUAL(bar.data[0][0][2], complex(31.0, -31.0));
    BOOST_CHECK_EQUAL(bar.data[0][1][2], complex(37.0, -37.0));

    // Ensure we catch an operation between two nonconforming states
    suzerain::ComplexState<double> baz(2, 2, 2);
    BOOST_CHECK_THROW(foo.scaleAddScaled(2.0, 3.0, &baz),
                      suzerain::logic_error);
}

BOOST_AUTO_TEST_CASE( comparison_and_assignment )
{
    suzerain::ComplexState<double> foo(1, 2, 3);

    typedef std::complex<double> complex;
    foo.data[0][0][0] = complex( 2.0, - 2.0);
    foo.data[0][1][0] = complex( 3.0, - 3.0);
    foo.data[0][0][1] = complex( 5.0, - 5.0);
    foo.data[0][1][1] = complex( 7.0, - 7.0);
    foo.data[0][0][2] = complex(11.0, -11.0);
    foo.data[0][1][2] = complex(13.0, -13.0);

    suzerain::ComplexState<double> bar(1, 2, 3);
    bar.data = foo.data;

    BOOST_CHECK_EQUAL(true, bar.data == foo.data);

    bar.data[0][1][0] = complex(555.0, -555.0);
    bar.data[0][1][1] = 777.0;

    BOOST_CHECK_EQUAL(false, bar.data == foo.data);

    BOOST_CHECK_EQUAL(foo.data[0][0][0], complex( 2.0, - 2.0));
    BOOST_CHECK_EQUAL(foo.data[0][1][0], complex( 3.0, - 3.0));
    BOOST_CHECK_EQUAL(foo.data[0][0][1], complex( 5.0, - 5.0));
    BOOST_CHECK_EQUAL(foo.data[0][1][1], complex( 7.0, - 7.0));
    BOOST_CHECK_EQUAL(foo.data[0][0][2], complex(11.0, -11.0));
    BOOST_CHECK_EQUAL(foo.data[0][1][2], complex(13.0, -13.0));

    BOOST_CHECK_EQUAL(bar.data[0][0][0], complex(  2.0, -  2.0));
    BOOST_CHECK_EQUAL(bar.data[0][1][0], complex(555.0, -555.0));
    BOOST_CHECK_EQUAL(bar.data[0][0][1], complex(  5.0, -  5.0));
    BOOST_CHECK_EQUAL(bar.data[0][1][1], complex(777.0,    0.0));
    BOOST_CHECK_EQUAL(bar.data[0][0][2], complex( 11.0, - 11.0));
    BOOST_CHECK_EQUAL(bar.data[0][1][2], complex( 13.0, - 13.0));
}

BOOST_AUTO_TEST_SUITE_END()
