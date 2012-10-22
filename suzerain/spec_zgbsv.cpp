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
      equil(false),
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
                  ( no_case["zgbsvx" ] [ ref(s.method) = spec_zgbsv::zgbsvx  ]
                    >> -(   ( char_(',') >> no_case["equil"] >> '='
                                         >> no_case[bool_]  [ref(s.equil) = _1])
                        )
                  )
                | ( no_case["zgbsv"  ] [ ref(s.method) = spec_zgbsv::zgbsv   ])
                | ( no_case["zcgbsvx"] [ ref(s.method) = spec_zgbsv::zcgbsvx ]
                    >> -(   ( char_(',') >> no_case["reuse"] >> '='
                                         >> no_case[bool_]  [ref(s.reuse) = _1])
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
    std::ostringstream os;
    os << std::boolalpha;

    if        (this->method == spec_zgbsv::zgbsv  ) {
        os << "zgbsv";
    } else if (this->method == spec_zgbsv::zgbsvx ) {
        os << "zgbsvx"
           << ",equil=" << this->equil;
    } else if (this->method == spec_zgbsv::zcgbsvx) {
        os << "zcgbsvx"
           << ",reuse=" << this->reuse
           << ",aiter=" << this->aiter
           << ",siter=" << this->siter
           << ",diter=" << this->diter
           << ",tolsc=" << this->tolsc;
    } else {
        os << "UNKNOWN";
    }

    return os.str();
}

bool spec_zgbsv::operator==(const spec_zgbsv &that) const
{
    if (this->method != that.method)
        return false;

    switch (this->method) {  // Additional stipulations per method

    default:
        return true;

    case spec_zgbsv::zgbsvx:
        return this->equil == that.equil;

    case spec_zgbsv::zcgbsvx:
        return this->reuse  == that.reuse
            && this->aiter  == that.aiter
            && this->siter  == that.siter
            && this->diter  == that.diter
            && this->tolsc  == that.tolsc;

    }
}

} // end namespace suzerain
