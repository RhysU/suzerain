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
 * exprparse.hpp: constant arithmetic expression evaluation
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <algorithm>
#include <cstring>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/version.hpp>
#include <suzerain/exprgrammar.hpp>
#include <suzerain/exprparse.hpp>

namespace suzerain {


namespace detail {

#if BOOST_VERSION < 104100 // Spirit.Qi 2.1+ unavailable so use lexical_cast

template<typename FPT, typename StringType>
FPT exprparse_impl(const StringType& s)
{
    FPT result;
    try {
        result = boost::lexical_cast<FPT>(s);
    } catch (boost::bad_lexical_cast &) {
        ostringstream m;
        m << "exprparse lexical_cast error on input '";
        copy(iter, end, ostream_iterator<const char>(m));
        m << '\'';
        throw invalid_argument(m.str());
    }
    return result;
}

#else  // BOOST_VERSION >= 104100 so use suzerain::exprgrammar

template<typename FPT>
FPT exprparse_impl(const char *s)
{
    using namespace std;
    const char *       iter = s;
    const char * const end  = iter + strlen(s);
    FPT result;

    if (!exprgrammar::parse(iter, end, result)) {
        ostringstream m;
        m << "exprparse error on input '";
        copy(iter, end, ostream_iterator<const char>(m));
        m << '\'';
        throw invalid_argument(m.str());
    }

    if (iter != end) {
        ostringstream m;
        m << "exprparse halted at position ";
        m << distance(s, iter);
        m << " on input '";
        copy(iter, end, ostream_iterator<const char>(m));
        m << '\'';
        throw invalid_argument(m.str());
    }

    return result;
}

template<typename FPT>
FPT exprparse_impl(const std::string &s)
{
    using namespace std;
    string::const_iterator       iter = s.begin();
    string::const_iterator const end  = s.end();
    FPT result;

    if (!exprgrammar::parse(iter, end, result)) {
        ostringstream m;
        m << "exprparse error on input '";
        copy(iter, end, ostream_iterator<string::value_type>(m));
        m << '\'';
        throw invalid_argument(m.str());
    }

    if (iter != end) {
        ostringstream m;
        m << "exprparse halted at position ";
        m << distance(s.begin(), iter);
        m << " on input '";
        copy(iter, end, ostream_iterator<string::value_type>(m));
        m << '\'';
        throw invalid_argument(m.str());
    }

    return result;
}

#endif // BOOST_VERSION

} // end namespace detail

// To mitigate long compilation times associated with Boost Spirit.Qi parsers,
// only a small number of template instantiations are provided.  These are
// compiled once and may be linked quickly.

void exprparse(const char *s, float& v)
{
    v = detail::exprparse_impl<float>(s);
}

void exprparse(const char *s, double& v)
{
    v = detail::exprparse_impl<double>(s);
}

void exprparse(const char *s, long double& v)
{
    v = detail::exprparse_impl<long double>(s);
}

void exprparse(const std::string& s, float& v)
{
    v = detail::exprparse_impl<float>(s);
}

void exprparse(const std::string& s, double& v)
{
    v = detail::exprparse_impl<double>(s);
}

void exprparse(const std::string& s, long double& v)
{
    v = detail::exprparse_impl<long double>(s);
}

} // end namespace suzerain
