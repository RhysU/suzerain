/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

#include <suzerain/fftw.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

BOOST_AUTO_TEST_SUITE( rigor_suite )

using namespace suzerain::fftw;

BOOST_AUTO_TEST_CASE( compare_against_flags )
{
    BOOST_CHECK_EQUAL(FFTW_ESTIMATE,    estimate);
    BOOST_CHECK_EQUAL(FFTW_MEASURE,     measure);
    BOOST_CHECK_EQUAL(FFTW_PATIENT,     patient);
    BOOST_CHECK_EQUAL(FFTW_EXHAUSTIVE,  exhaustive);
    BOOST_CHECK_EQUAL(FFTW_WISDOM_ONLY, wisdom_only);
}

BOOST_AUTO_TEST_CASE( from_flags )
{
    BOOST_CHECK_EQUAL(estimate,    rigor_from(FFTW_ESTIMATE));
    BOOST_CHECK_EQUAL(measure,     rigor_from(FFTW_MEASURE));
    BOOST_CHECK_EQUAL(patient,     rigor_from(FFTW_PATIENT));
    BOOST_CHECK_EQUAL(exhaustive,  rigor_from(FFTW_EXHAUSTIVE));
    BOOST_CHECK_EQUAL(wisdom_only, rigor_from(FFTW_WISDOM_ONLY));
}

BOOST_AUTO_TEST_CASE( from_constchar )
{
    BOOST_CHECK_EQUAL(estimate,    rigor_from("estimate"));
    BOOST_CHECK_EQUAL(estimate,    rigor_from("ESTIMATE"));
    BOOST_CHECK_EQUAL(estimate,    rigor_from("es"));
    BOOST_CHECK_EQUAL(exhaustive,  rigor_from("exhaustive"));
    BOOST_CHECK_EQUAL(exhaustive,  rigor_from("EXHAUSTIVE"));
    BOOST_CHECK_EQUAL(exhaustive,  rigor_from("EX"));
    BOOST_CHECK_EQUAL(measure,     rigor_from("measure"));
    BOOST_CHECK_EQUAL(measure,     rigor_from("MEAsURE"));
    BOOST_CHECK_EQUAL(measure,     rigor_from("meaSure"));
    BOOST_CHECK_EQUAL(measure,     rigor_from("m"));
    BOOST_CHECK_EQUAL(patient,     rigor_from("patient"));
    BOOST_CHECK_EQUAL(patient,     rigor_from("PATIENT"));
    BOOST_CHECK_EQUAL(patient,     rigor_from("P"));
    BOOST_CHECK_EQUAL(estimate,    rigor_from("e"));
    BOOST_CHECK_EQUAL(estimate,    rigor_from("E"));
    BOOST_CHECK_EQUAL(wisdom_only, rigor_from("wisdom_only"));
    BOOST_CHECK_EQUAL(wisdom_only, rigor_from("WISDOM_ONLY"));
    BOOST_CHECK_EQUAL(wisdom_only, rigor_from("w"));
    BOOST_CHECK_EQUAL(wisdom_only, rigor_from("W"));

    BOOST_CHECK_EQUAL(measure, rigor_from(""));                 // Unknown
    BOOST_CHECK_EQUAL(measure, rigor_from((const char *)NULL)); // Unknown
}

BOOST_AUTO_TEST_CASE( const_char )
{
    BOOST_CHECK( 0 == strcmp("estimate",    c_str(estimate))    );
    BOOST_CHECK( 0 == strcmp("measure",     c_str(measure))     );
    BOOST_CHECK( 0 == strcmp("patient",     c_str(patient))     );
    BOOST_CHECK( 0 == strcmp("exhaustive",  c_str(exhaustive))  );
    BOOST_CHECK( 0 == strcmp("wisdom_only", c_str(wisdom_only)) );
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
    const rigor rigors[] = {
            estimate, patient, measure, exhaustive, wisdom_only
    };
    for (std::size_t i = 0; i < sizeof(rigors)/sizeof(rigors[0]); ++i) {
        BOOST_CHECK_EQUAL(rigors[i], rigor_from(c_str(rigors[i])));
    }
}

BOOST_AUTO_TEST_SUITE_END()
