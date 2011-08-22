/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
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
 *
 * exprparse_impl.hpp: constant arithmetic expression evaluation (private)
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_EXPRPARSE_IMPL_HPP
#define __SUZERAIN_EXPRPARSE_IMPL_HPP

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
