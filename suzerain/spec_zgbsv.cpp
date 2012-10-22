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

#include <suzerain/spec_zgbsv.hpp>

#include <algorithm>
#include <cctype>
#include <ios>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include <boost/spirit/version.hpp>
#if !defined(SPIRIT_VERSION) || SPIRIT_VERSION < 0x2010
#error "At least Spirit version 2.1 required"
#endif
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>

namespace suzerain {

spec_zgbsv::spec_zgbsv()
    : method(zgbsvx),
      reuse(true),
      aiter(10),
      siter(30),
      diter(100),
      tolsc(1)
{}

spec_zgbsv::spec_zgbsv(const std::string& spec)
{
    using namespace std;

    // Default construct instance
    spec_zgbsv s;

    // Advance past any leading whitespace
    std::string::const_iterator iter = spec.begin(), end = spec.end();
    while (iter != end && isspace(*iter)) ++iter;

    // If nontrivial content, parse with a Spirit grammar modifying defaults
    if (iter != end) {
        using boost::phoenix::ref;
        using boost::spirit::qi::_1;
        using boost::spirit::qi::bool_;
        using boost::spirit::qi::char_;
        using boost::spirit::qi::double_;
        using boost::spirit::qi::int_;
        using boost::spirit::qi::no_case;
        using boost::spirit::qi::no_case;

        const bool r = boost::spirit::qi::phrase_parse(iter, end, (
            // Grammar Begin
                  ( no_case["zgbsvx" ] [ ref(s.method) = spec_zgbsv::zgbsvx  ])
                | ( no_case["zgbsv"  ] [ ref(s.method) = spec_zgbsv::zgbsv   ])
                | ( no_case["zcgbsvx"] [ ref(s.method) = spec_zgbsv::zcgbsvx ]
                    >> -(   ( char_(',') >> no_case["reuse"] >> '='
                                         >> bool_  [ref(s.reuse) = _1])
                          ^ ( char_(',') >> no_case["aiter"] >> '='
                                         >> int_   [ref(s.aiter) = _1])
                          ^ ( char_(',') >> no_case["siter"] >> '='
                                         >> int_   [ref(s.siter) = _1])
                          ^ ( char_(',') >> no_case["diter"] >> '='
                                         >> int_   [ref(s.diter) = _1])
                          ^ ( char_(',') >> no_case["tolsc"] >> '='
                                         >> double_[ref(s.tolsc) = _1])
                        )
                  )
            // Grammar End
            ), boost::spirit::ascii::space);

        // Ensure the grammar matched and the entire specification was consumed
        if (!r || iter != end) {
            ostringstream os;
            os << "spec_zgbsv specification '" << spec
               << "' invalid beginning with '";
            copy(iter, end, ostream_iterator<string::value_type>(os));
            os << '\'';
            throw invalid_argument(os.str());
        }
    }

    // Assign the mutated defaults to *this
    *this = s;
}

spec_zgbsv::operator std::string () const
{
    switch (this->method) {
        default:                   return "UNKNOWN"  ;
        case spec_zgbsv::zgbsv:    return "zgbsv"    ;
        case spec_zgbsv::zgbsvx:   return "zgbsvx"   ;
        case spec_zgbsv::zcgbsvx:  /* Fall through */;
    }

    std::ostringstream os;
    os << std::boolalpha
       << "zcgbsvx"
       << ",reuse=" << this->reuse
       << ",aiter=" << this->aiter
       << ",siter=" << this->siter
       << ",diter=" << this->diter
       << ",tolsc=" << this->tolsc;

    return os.str();
}

bool spec_zgbsv::operator==(const spec_zgbsv &that) const
{
    switch (this->method) {
        default:                   return this->method == that.method;
        case spec_zgbsv::zcgbsvx:  /* Fall through */                ;
    }

    return this->method == that.method
        && this->reuse  == that.reuse
        && this->aiter  == that.aiter
        && this->diter  == that.diter
        && this->tolsc  == that.tolsc;
}

} // end namespace suzerain
