//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
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

#include <suzerain/specification_zgbsv.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

using suzerain::zgbsv_specification;

BOOST_AUTO_TEST_CASE( default_constructor )
{
    zgbsv_specification s;
    BOOST_CHECK_NE("UNKNOWN", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_empty )
{
    zgbsv_specification s(" ");  // Whitespace
    BOOST_CHECK(s.method());
    BOOST_CHECK_NE("",        (std::string) s);
    BOOST_CHECK_NE("UNKNOWN", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_zgbsv )
{
    zgbsv_specification s("zgbsv");
    BOOST_CHECK(s.method() == zgbsv_specification::zgbsv);
    BOOST_CHECK_EQUAL("zgbsv", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_zgbsvx )
{
    {
        zgbsv_specification s("zgbsvx");  // Uppercase
        BOOST_CHECK(s.method() == zgbsv_specification::zgbsvx);
        BOOST_CHECK_EQUAL("zgbsvx,equil=false", (std::string) s);
    }

    {
        zgbsv_specification s("ZGBSVX,EQUIL=TRUE");  // Uppercase
        BOOST_CHECK(s.method() == zgbsv_specification::zgbsvx);
        BOOST_CHECK_EQUAL("zgbsvx,equil=true", (std::string) s);
    }
}

BOOST_AUTO_TEST_CASE( construct_zcgbsvx )
{
    {
        zgbsv_specification s("zcgbsvx");
        BOOST_CHECK(s.method() == zgbsv_specification::zcgbsvx);
        BOOST_CHECK_NE("zcgbsvx", (std::string) s);   // Much more required
    }

    {
        std::string str("zcgbsvx,reuse=true,aiter=2,siter=3,diter=4,tolsc=5");
        zgbsv_specification s(str);
        BOOST_CHECK(s.method() == zgbsv_specification::zcgbsvx);  // As constructed...?
        BOOST_CHECK_EQUAL(s.reuse(), true);
        BOOST_CHECK_EQUAL(s.aiter(), 2);
        BOOST_CHECK_EQUAL(s.siter(), 3);
        BOOST_CHECK_EQUAL(s.diter(), 4);
        BOOST_CHECK_EQUAL(s.tolsc(), 5);
        BOOST_CHECK_EQUAL(str, (std::string) s);       // Roundtripped?
    }

    {
        // Out of order
        zgbsv_specification a("zcgbsvx , reuse = true\t\n,\t\nAITER=\t \t2");
        zgbsv_specification b("zcgbsvx,reuse=true,aiter=2");
        zgbsv_specification c("zcgbsvx,aiter=2,reuse=true");
        BOOST_CHECK(a == b);
        BOOST_CHECK(b == c);
        BOOST_CHECK(a == c);
    }
}
