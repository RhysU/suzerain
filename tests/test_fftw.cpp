#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/fftw.hpp>

BOOST_AUTO_TEST_SUITE( rigor_suite )

using namespace suzerain::fftw;

BOOST_AUTO_TEST_CASE( compare_against_flags )
{
    BOOST_CHECK_EQUAL(FFTW_ESTIMATE,   estimate);
    BOOST_CHECK_EQUAL(FFTW_MEASURE,    measure);
    BOOST_CHECK_EQUAL(FFTW_PATIENT,    patient);
    BOOST_CHECK_EQUAL(FFTW_EXHAUSTIVE, exhaustive);
}

BOOST_AUTO_TEST_CASE( from_flags )
{
    BOOST_CHECK_EQUAL(estimate,   rigor_from(FFTW_ESTIMATE));
    BOOST_CHECK_EQUAL(measure,    rigor_from(FFTW_MEASURE));
    BOOST_CHECK_EQUAL(patient,    rigor_from(FFTW_PATIENT));
    BOOST_CHECK_EQUAL(exhaustive, rigor_from(FFTW_EXHAUSTIVE));
}

BOOST_AUTO_TEST_CASE( from_constchar )
{
    BOOST_CHECK_EQUAL(estimate,   rigor_from("estimate"));
    BOOST_CHECK_EQUAL(estimate,   rigor_from("ESTIMATE"));
    BOOST_CHECK_EQUAL(estimate,   rigor_from("es"));
    BOOST_CHECK_EQUAL(exhaustive, rigor_from("exhaustive"));
    BOOST_CHECK_EQUAL(exhaustive, rigor_from("EXHAUSTIVE"));
    BOOST_CHECK_EQUAL(exhaustive, rigor_from("EX"));
    BOOST_CHECK_EQUAL(measure,    rigor_from("measure"));
    BOOST_CHECK_EQUAL(measure,    rigor_from("MEAsURE"));
    BOOST_CHECK_EQUAL(measure,    rigor_from("meaSure"));
    BOOST_CHECK_EQUAL(measure,    rigor_from("m"));
    BOOST_CHECK_EQUAL(patient,    rigor_from("patient"));
    BOOST_CHECK_EQUAL(patient,    rigor_from("PATIENT"));
    BOOST_CHECK_EQUAL(patient,    rigor_from("P"));

    BOOST_CHECK_EQUAL(measure, rigor_from("e"));  // Ambiguous
    BOOST_CHECK_EQUAL(measure, rigor_from("E"));  // Ambiguous

    BOOST_CHECK_EQUAL(measure, rigor_from(""));                 // Unknown
    BOOST_CHECK_EQUAL(measure, rigor_from((const char *)NULL)); // Unknown
}

BOOST_AUTO_TEST_CASE( const_char )
{
    BOOST_CHECK( 0 == strcmp("estimate",   c_str(estimate)) );
    BOOST_CHECK( 0 == strcmp("measure",    c_str(measure)) );
    BOOST_CHECK( 0 == strcmp("patient",    c_str(patient)) );
    BOOST_CHECK( 0 == strcmp("exhaustive", c_str(exhaustive)) );
}

BOOST_AUTO_TEST_CASE( ostream )
{
    // No need to check all values, just one will do
    std::ostringstream oss;
    oss << estimate;
    BOOST_CHECK_EQUAL(oss.str(), "estimate");
}

BOOST_AUTO_TEST_CASE( round_trip )
{
    const rigor rigors[] = { estimate, patient, measure, exhaustive };
    for (int i = 0; i < sizeof(rigors)/sizeof(rigors[0]); ++i) {
        BOOST_CHECK_EQUAL(rigors[i], rigor_from(c_str(rigors[i])));
    }
}

BOOST_AUTO_TEST_SUITE_END()
