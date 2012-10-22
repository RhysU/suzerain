//--------------------------------------------------------------------------
//
// Copyright (C) 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// spec_zgbsv.cpp: Encapsulates complex-valued banded solve options
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/spec_zgbsv.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

using suzerain::spec_zgbsv;

BOOST_AUTO_TEST_CASE( default_constructor )
{
    spec_zgbsv s;
    BOOST_CHECK_NE("UNKNOWN", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_empty )
{
    spec_zgbsv s(" ");  // Whitespace
    BOOST_CHECK(s.method);
    BOOST_CHECK_NE("",        (std::string) s);
    BOOST_CHECK_NE("UNKNOWN", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_zgbsv )
{
    spec_zgbsv s("zgbsv");
    BOOST_CHECK(s.method == spec_zgbsv::zgbsv);
    BOOST_CHECK_EQUAL("zgbsv", (std::string) s);
}

BOOST_AUTO_TEST_CASE( construct_zgbsvx )
{
    {
        spec_zgbsv s("zgbsvx");  // Uppercase
        BOOST_CHECK(s.method == spec_zgbsv::zgbsvx);
        BOOST_CHECK_EQUAL("zgbsvx,equil=false", (std::string) s);
    }

    {
        spec_zgbsv s("ZGBSVX,EQUIL=TRUE");  // Uppercase
        BOOST_CHECK(s.method == spec_zgbsv::zgbsvx);
        BOOST_CHECK_EQUAL("zgbsvx,equil=true", (std::string) s);
    }
}

BOOST_AUTO_TEST_CASE( construct_zcgbsvx )
{
    {
        spec_zgbsv s("zcgbsvx");
        BOOST_CHECK(s.method == spec_zgbsv::zcgbsvx);
        BOOST_CHECK_NE("zcgbsvx", (std::string) s);   // Much more required
    }

    {
        std::string str("zcgbsvx,reuse=true,aiter=2,siter=3,diter=4,tolsc=5");
        spec_zgbsv s(str);
        BOOST_CHECK(s.method == spec_zgbsv::zcgbsvx);  // As constructed...?
        BOOST_CHECK_EQUAL(s.reuse, true);
        BOOST_CHECK_EQUAL(s.aiter, 2);
        BOOST_CHECK_EQUAL(s.siter, 3);
        BOOST_CHECK_EQUAL(s.diter, 4);
        BOOST_CHECK_EQUAL(s.tolsc, 5);
        BOOST_CHECK_EQUAL(str, (std::string) s);       // Roundtripped?
    }

    {
        // Out of order
        spec_zgbsv a("zcgbsvx , reuse = true\t\n,\t\nAITER=\t \t2");
        spec_zgbsv b("zcgbsvx,reuse=true,aiter=2");
        spec_zgbsv c("zcgbsvx,aiter=2,reuse=true");
        BOOST_CHECK(a == b);
        BOOST_CHECK(b == c);
        BOOST_CHECK(a == c);
    }
}
