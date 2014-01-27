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

using suzerain::specification_zgbsv;

BOOST_AUTO_TEST_CASE( default_constructor )
{
    specification_zgbsv s;
    BOOST_CHECK_NE("UNKNOWN", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_empty )
{
    specification_zgbsv s(" ");  // Whitespace
    BOOST_CHECK(s.method());
    BOOST_CHECK_NE("",        (std::string) s);
    BOOST_CHECK_NE("UNKNOWN", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_zgbsv )
{
    specification_zgbsv s("zgbsv");
    BOOST_CHECK(s.method() == specification_zgbsv::zgbsv);
    BOOST_CHECK_EQUAL("zgbsv", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_zgbsvx )
{
    {
        specification_zgbsv s("zgbsvx");  // Uppercase
        BOOST_CHECK(s.method() == specification_zgbsv::zgbsvx);
        BOOST_CHECK_EQUAL("zgbsvx,equil=false", (std::string) s);
    }

    {
        specification_zgbsv s("ZGBSVX,EQUIL=TRUE");  // Uppercase
        BOOST_CHECK(s.method() == specification_zgbsv::zgbsvx);
        BOOST_CHECK_EQUAL("zgbsvx,equil=true", (std::string) s);
    }
}

BOOST_AUTO_TEST_CASE( construct_zcgbsvx )
{
    {
        specification_zgbsv s("zcgbsvx");
        BOOST_CHECK(s.method() == specification_zgbsv::zcgbsvx);
        BOOST_CHECK_NE("zcgbsvx", (std::string) s);   // Much more required
    }

    {
        std::string str("zcgbsvx,reuse=true,aiter=2,siter=3,diter=4,tolsc=5");
        specification_zgbsv s(str);
        BOOST_CHECK(s.method() == specification_zgbsv::zcgbsvx);  // As constructed...?
        BOOST_CHECK_EQUAL(s.reuse(), true);
        BOOST_CHECK_EQUAL(s.aiter(), 2);
        BOOST_CHECK_EQUAL(s.siter(), 3);
        BOOST_CHECK_EQUAL(s.diter(), 4);
        BOOST_CHECK_EQUAL(s.tolsc(), 5);
        BOOST_CHECK_EQUAL(str, (std::string) s);       // Roundtripped?
    }

    {
        // Out of order
        specification_zgbsv a("zcgbsvx , reuse = true\t\n,\t\nAITER=\t \t2");
        specification_zgbsv b("zcgbsvx,reuse=true,aiter=2");
        specification_zgbsv c("zcgbsvx,aiter=2,reuse=true");
        BOOST_CHECK(a == b);
        BOOST_CHECK(b == c);
        BOOST_CHECK(a == c);
    }
}
