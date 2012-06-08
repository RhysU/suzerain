//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// exprparse_impl.hpp: constant arithmetic expression evaluation (private)
// $Id$

#ifndef SUZERAIN_EXPRPARSE_IMPL_HPP
#define SUZERAIN_EXPRPARSE_IMPL_HPP

#include <algorithm>
#include <cstring>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <suzerain/exprgrammar.hpp>

namespace suzerain {

namespace detail {

// Forward declaration to silence Intel remark #1418
template<typename FPT>
FPT exprparse_impl(const char *s, const char *name = NULL);

template<typename FPT>
FPT exprparse_impl(const char *s, const char *name)
{
    using namespace std;
    const char *       iter = s;
    const char * const end  = iter + strlen(s);
    FPT result;

    if (!exprgrammar::parse(iter, end, result)) {
        ostringstream m;
        if (name) {
            m << "exprparse error in " << name << " at '";
        } else {
            m << "exprparse error at '";
        }
        copy(iter, end, ostream_iterator<const char>(m));
        m << '\'';
        throw invalid_argument(m.str());
    }

    if (iter != end) {
        ostringstream m;
        m << "exprparse halted at position " << distance(s, iter) << " in ";
        if (name) {
            m << name << ' ';
        }
        m << '\'' << s << '\'';
        throw invalid_argument(m.str());
    }

    return result;
}

// Forward declaration to silence Intel remark #1418
template<typename FPT>
FPT exprparse_impl(const std::string &s, const char *name = NULL);

template<typename FPT>
FPT exprparse_impl(const std::string &s, const char *name)
{
    using namespace std;
    string::const_iterator       iter = s.begin();
    string::const_iterator const end  = s.end();
    FPT result;

    if (!exprgrammar::parse(iter, end, result)) {
        ostringstream m;
        if (name) {
            m << "exprparse error in " << name << " at '";
        } else {
            m << "exprparse error at '";
        }
        copy(iter, end, ostream_iterator<string::value_type>(m));
        m << '\'';
        throw invalid_argument(m.str());
    }

    if (iter != end) {
        ostringstream m;
        m << "exprparse halted at position "
          << distance(s.begin(), iter) << " in ";
        if (name) {
            m << name << ' ';
        }
        m << '\'' << s << '\'';
        throw invalid_argument(m.str());
    }

    return result;
}

} // end namespace detail

} // end namespace suzerain

#endif
