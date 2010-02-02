#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/complex.hpp>

BOOST_AUTO_TEST_SUITE( fill_multi_array )

using boost::array_view_gen;
using boost::extents;
using boost::indices;
using boost::multi_array;
using boost::multi_array_ref;
using boost::multi_array_types::index;
using boost::multi_array_types::index_range;
using boost::scoped_array;
using suzerain::multi_array::fill;

const std::size_t NX = 3, NY = 4, NZ = 5, NZZ = 6;

typedef boost::mpl::list<int, float, double> element_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_multi_array_1D, T, element_test_types)
{
    // Fill boost::multi_array
    multi_array<T,1> foo(extents[NX]);
    fill(foo,123);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        BOOST_CHECK_EQUAL(123, foo[i]); // Ensure multi_array filled
    }

    // Fill view of boost::multi_array
    typename array_view_gen<multi_array<T,1>,1>::type bar
        = foo[ indices[index_range(0,foo.shape()[0]-1)] ];
    fill(bar,456);
    for (index i = 0; i < bar.shape()[0]; ++i) {
        BOOST_CHECK_EQUAL(456, foo[i]); // Ensure multi_array filled
        BOOST_CHECK_EQUAL(456, bar[i]); // Ensure view filled
    }
    for (index i = bar.shape()[0]; i < foo.shape()[0]; ++i) {
        BOOST_CHECK_EQUAL(123, foo[i]); // Ensure non-view untouched
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_multi_array_ref_1D, T, element_test_types)
{
    // Fill boost::multi_array_ref
    scoped_array<T> raw(new T[NX]);
    multi_array_ref<T,1> foo(raw.get(),extents[NX]);
    fill(foo,123);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        BOOST_CHECK_EQUAL(123, foo[i]); // Ensure multi_array_ref filled
    }

    // Fill view of boost::multi_array_ref
    typename array_view_gen<multi_array_ref<T,1>,1>::type bar = foo[
        indices[index_range(0,foo.shape()[0]-1)]
    ];
    fill(bar,456);
    for (index i = 0; i < bar.shape()[0]; ++i) {
        BOOST_CHECK_EQUAL(456, foo[i]); // Ensure multi_array_ref filled
        BOOST_CHECK_EQUAL(456, bar[i]); // Ensure view filled
    }
    for (index i = bar.shape()[0]; i < foo.shape()[0]; ++i) {
        BOOST_CHECK_EQUAL(123, foo[i]); // Ensure non-view untouched
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_multi_array_2D, T, element_test_types)
{
    // Fill boost::multi_array
    multi_array<T,2> foo(extents[NX][NY]);
    fill(foo,123);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            BOOST_CHECK_EQUAL(123, foo[i][j]); // Ensure multi_array filled
        }
    }

    // Fill view of boost::multi_array
    typename array_view_gen<multi_array<T,2>,2>::type bar = foo[
        indices[index_range(0,foo.shape()[0]-1)]
               [index_range(0,foo.shape()[1]-1)]
    ];
    fill(bar,456);
    for (index i = 0; i < bar.shape()[0]; ++i) {
        for (index j = 0; j < bar.shape()[1]; ++j) {
            BOOST_CHECK_EQUAL(456, foo[i][j]); // Ensure multi_array filled
            BOOST_CHECK_EQUAL(456, bar[i][j]); // Ensure view filled
        }
        for (index j = bar.shape()[1]; j < foo.shape()[1]; ++j) {
            BOOST_CHECK_EQUAL(123, foo[i][j]); // Ensure non-view untouched
        }
    }
    for (index i = bar.shape()[0]; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            BOOST_CHECK_EQUAL(123, foo[i][j]); // Ensure non-view untouched
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_multi_array_ref_2D, T, element_test_types)
{
    // Fill boost::multi_array_ref
    scoped_array<T> raw(new T[NX*NY]);
    multi_array_ref<T,2> foo(raw.get(),extents[NX][NY]);
    fill(foo,123);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            BOOST_CHECK_EQUAL(123, foo[i][j]); // Ensure multi_array_ref filled
        }
    }

    // Fill view of boost::multi_array_ref
    typename array_view_gen<multi_array<T,2>,2>::type bar = foo[
        indices[index_range(0,foo.shape()[0]-1)]
               [index_range(0,foo.shape()[1]-1)]
    ];
    fill(bar,456);
    for (index i = 0; i < bar.shape()[0]; ++i) {
        for (index j = 0; j < bar.shape()[1]; ++j) {
            BOOST_CHECK_EQUAL(456, foo[i][j]); // Ensure multi_array_ref filled
            BOOST_CHECK_EQUAL(456, bar[i][j]); // Ensure view filled
        }
        for (index j = bar.shape()[1]; j < foo.shape()[1]; ++j) {
            BOOST_CHECK_EQUAL(123, foo[i][j]); // Ensure non-view untouched
        }
    }
    for (index i = bar.shape()[0]; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            BOOST_CHECK_EQUAL(123, foo[i][j]); // Ensure non-view untouched
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_multi_array_3D, T, element_test_types)
{
    // Fill boost::multi_array
    multi_array<T,3> foo(extents[NX][NY][NZ]);
    fill(foo,123);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            for (index k = 0; k < foo.shape()[2]; ++k) {
                BOOST_CHECK_EQUAL(123, foo[i][j][k]); // filled?
            }
        }
    }

    // Fill view of boost::multi_array
    typename array_view_gen<multi_array<T,3>,3>::type bar = foo[
        indices[index_range(0,foo.shape()[0]-1)]
               [index_range(0,foo.shape()[1]-1)]
               [index_range(0,foo.shape()[2]-1)]
    ];
    fill(bar,456);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            for (index k = 0; k < foo.shape()[2]; ++k) {

                if (   i < bar.shape()[0]
                    && j < bar.shape()[1]
                    && k < bar.shape()[2]) {

                    BOOST_CHECK_EQUAL(456, foo[i][j][k]); // orig filled?
                    BOOST_CHECK_EQUAL(456, bar[i][j][k]); // view filled?
                } else {
                    BOOST_CHECK_EQUAL(123, foo[i][j][k]); // orig untouched?
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_multi_array_ref_3D, T, element_test_types)
{
    // Fill boost::multi_array_ref
    scoped_array<T> raw(new T[NX*NY*NZ]);
    multi_array_ref<T,3> foo(raw.get(),extents[NX][NY][NZ]);
    fill(foo,123);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            for (index k = 0; k < foo.shape()[2]; ++k) {
                BOOST_CHECK_EQUAL(123, foo[i][j][k]); // filled?
            }
        }
    }

    // Fill view of boost::multi_array_ref
    typename array_view_gen<multi_array_ref<T,3>,3>::type bar = foo[
        indices[index_range(0,foo.shape()[0]-1)]
               [index_range(0,foo.shape()[1]-1)]
               [index_range(0,foo.shape()[2]-1)]
    ];
    fill(bar,456);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            for (index k = 0; k < foo.shape()[2]; ++k) {

                if (   i < bar.shape()[0]
                    && j < bar.shape()[1]
                    && k < bar.shape()[2]) {

                    BOOST_CHECK_EQUAL(456, foo[i][j][k]); // orig filled?
                    BOOST_CHECK_EQUAL(456, bar[i][j][k]); // view filled?
                } else {
                    BOOST_CHECK_EQUAL(123, foo[i][j][k]); // orig untouched?
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_multi_array_4D, T, element_test_types)
{
    // Fill boost::multi_array
    multi_array<T,4> foo(extents[NX][NY][NZ][NZZ]);
    fill(foo,123);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            for (index k = 0; k < foo.shape()[2]; ++k) {
                for (index l = 0; l < foo.shape()[3]; ++l) {
                    BOOST_CHECK_EQUAL(123, foo[i][j][k][l]); // filled?
                }
            }
        }
    }

    // Fill view of boost::multi_array
    typename array_view_gen<multi_array<T,4>,4>::type bar = foo[
        indices[index_range(0,foo.shape()[0]-1)]
               [index_range(0,foo.shape()[1]-1)]
               [index_range(0,foo.shape()[2]-1)]
               [index_range(0,foo.shape()[3]-1)]
    ];
    fill(bar,456);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            for (index k = 0; k < foo.shape()[2]; ++k) {
                for (index l = 0; l < foo.shape()[3]; ++l) {

                    if (   i < bar.shape()[0]
                        && j < bar.shape()[1]
                        && k < bar.shape()[2]
                        && l < bar.shape()[3]) {

                        BOOST_CHECK_EQUAL(456, foo[i][j][k][l]); // orig ?
                        BOOST_CHECK_EQUAL(456, bar[i][j][k][l]); // view ?
                    } else {
                        BOOST_CHECK_EQUAL(123, foo[i][j][k][l]); // untouched?
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fill_multi_array_ref_4D, T, element_test_types)
{
    // Fill boost::multi_array_ref
    scoped_array<T> raw(new T[NX*NY*NZ*NZZ]);
    multi_array_ref<T,4> foo(raw.get(),extents[NX][NY][NZ][NZZ]);
    fill(foo,123);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            for (index k = 0; k < foo.shape()[2]; ++k) {
                for (index l = 0; l < foo.shape()[3]; ++l) {
                    BOOST_CHECK_EQUAL(123, foo[i][j][k][l]); // filled?
                }
            }
        }
    }

    // Fill view of boost::multi_array_ref
    typename array_view_gen<multi_array_ref<T,4>,4>::type bar = foo[
        indices[index_range(0,foo.shape()[0]-1)]
               [index_range(0,foo.shape()[1]-1)]
               [index_range(0,foo.shape()[2]-1)]
               [index_range(0,foo.shape()[3]-1)]
    ];
    fill(bar,456);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        for (index j = 0; j < foo.shape()[1]; ++j) {
            for (index k = 0; k < foo.shape()[2]; ++k) {
                for (index l = 0; l < foo.shape()[3]; ++l) {

                    if (   i < bar.shape()[0]
                        && j < bar.shape()[1]
                        && k < bar.shape()[2]
                        && l < bar.shape()[3]) {

                        BOOST_CHECK_EQUAL(456, foo[i][j][k][l]); // orig ?
                        BOOST_CHECK_EQUAL(456, bar[i][j][k][l]); // view ?
                    } else {
                        BOOST_CHECK_EQUAL(123, foo[i][j][k][l]); // untouched?
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

// Checks that the suzerain::functional::assign specialization
// in <suzerain/complex.hpp> is picked up within
// suzerain::multi_array::fill as expected.
BOOST_AUTO_TEST_SUITE( fill_complex_multiarray )

using boost::extents;
using boost::multi_array;
using boost::multi_array_types::index;
using suzerain::multi_array::fill;

const std::size_t NX = 3;

BOOST_AUTO_TEST_CASE(fill_complex_multiarray)
{
    multi_array<std::complex<double>,1> foo(extents[NX]);

    // Fill std::complex elements with FFTW-like complex array
    const double z[2] = { 1, -1 };
    fill(foo, z);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        BOOST_CHECK_EQUAL(foo[i].real(),  1);
        BOOST_CHECK_EQUAL(foo[i].imag(), -1);
    }

    // Fill std::complex elements with real value
    // Imaginary parts should be set to zero as well
    fill(foo, 0);
    for (index i = 0; i < foo.shape()[0]; ++i) {
        BOOST_CHECK_EQUAL(foo[i].real(), 0);
        BOOST_CHECK_EQUAL(foo[i].imag(), 0);
    }
}

BOOST_AUTO_TEST_SUITE_END()
