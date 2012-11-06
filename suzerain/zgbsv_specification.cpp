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
// zgbsv_specification.cpp: Encapsulates complex-valued banded solve options
// $Id$

#include <suzerain/zgbsv_specification.hpp>

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

zgbsv_specification::zgbsv_specification()  // Defaults behave much like zgbsvx
    : method_(zcgbsvx),
      equil_(false),
      reuse_(false),
      aiter_(1),
      siter_(-1),
      diter_(5),
      tolsc_(0)
{}

zgbsv_specification::zgbsv_specification(const std::string& spec)
{
    using namespace std;

    // Default construct instance
    zgbsv_specification s;

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

        const bool r = boost::spirit::qi::phrase_parse(iter, end, (
            // Grammar Begin
                  ( no_case["zgbsvx" ] [ ref(s.method_) = zgbsv_specification::zgbsvx  ]
                    >> -(   ( char_(',') >> no_case["equil"] >> '='
                                         >> no_case[bool_]  [ref(s.equil_) = _1])
                        )
                  )
                | ( no_case["zgbsv"  ] [ ref(s.method_) = zgbsv_specification::zgbsv   ])
                | ( no_case["zcgbsvx"] [ ref(s.method_) = zgbsv_specification::zcgbsvx ]
                    >> -(   ( char_(',') >> no_case["reuse"] >> '='
                                         >> no_case[bool_]  [ref(s.reuse_) = _1])
                          ^ ( char_(',') >> no_case["aiter"] >> '='
                                         >> int_   [ref(s.aiter_) = _1])
                          ^ ( char_(',') >> no_case["siter"] >> '='
                                         >> int_   [ref(s.siter_) = _1])
                          ^ ( char_(',') >> no_case["diter"] >> '='
                                         >> int_   [ref(s.diter_) = _1])
                          ^ ( char_(',') >> no_case["tolsc"] >> '='
                                         >> double_[ref(s.tolsc_) = _1])
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

zgbsv_specification::operator std::string () const
{
    std::ostringstream os;
    os << std::boolalpha;

    if        (method_ == zgbsv_specification::zgbsv  ) {
        os << "zgbsv";
    } else if (method_ == zgbsv_specification::zgbsvx ) {
        os << "zgbsvx"
           << ",equil=" << equil_;
    } else if (method_ == zgbsv_specification::zcgbsvx) {
        os << "zcgbsvx"
           << ",reuse=" << reuse_
           << ",aiter=" << aiter_
           << ",siter=" << siter_
           << ",diter=" << diter_
           << ",tolsc=" << tolsc_;
    } else {
        os << "UNKNOWN";
    }

    return os.str();
}

bool zgbsv_specification::operator==(const zgbsv_specification &that) const
{
    if (method_ != that.method_)
        return false;

    switch (this->method_) {  // Additional stipulations by method

    default:
        return true;

    case zgbsv_specification::zgbsvx:
        return equil_ == that.equil_;

    case zgbsv_specification::zcgbsvx:
        return reuse_  == that.reuse_
            && aiter_  == that.aiter_
            && siter_  == that.siter_
            && diter_  == that.diter_
#pragma warning(push,disable:1572)
            && tolsc_  == that.tolsc_;
#pragma warning(pop)

    }
}

} // end namespace suzerain
