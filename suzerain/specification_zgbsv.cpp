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

/** @file
 * @copydoc specification_zgbsv.hpp
 */

#include <suzerain/specification_zgbsv.hpp>

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

bool zgbsv_specification::in_place() const
{
    switch (method_) {
    case zgbsv_specification::zgbsv:   return true;
    case zgbsv_specification::zgbsvx:  return false;
    case zgbsv_specification::zcgbsvx: return false;
    default:                           return false;
    }
}

const char * zgbsv_specification::mname() const
{
    switch (method_) {
    case zgbsv_specification::zgbsv:   return "suzerain_lapackext_zgbsv";
    case zgbsv_specification::zgbsvx:  return "suzerain_lapack_zgbsvx";
    case zgbsv_specification::zcgbsvx: return "suzerain_lapackext_zcgbsvx";
    default:                           return "UNKNOWN";
    }
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
