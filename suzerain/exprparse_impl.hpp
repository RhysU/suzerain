//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_EXPRPARSE_IMPL_HPP
#define SUZERAIN_EXPRPARSE_IMPL_HPP

/** @file
 * Implementation logic built using the grammar in exprgrammar.hpp.
 */

#include <algorithm>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>

#include "exprgrammar.hpp"

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

// Prior forward declaration suppresses Intel warnings
template<typename T>
void exprparse_range_impl(const std::string& s,
                          T *valmin, T *valmax,
                          const T defaultmin, const T defaultmax,
                          const T absmin, const T absmax,
                          const char *name = NULL);

// Prior forward declaration suppresses Intel warnings
template<typename T>
void exprparse_range_impl(const std::string& s,
                          T *valmin, T *valmax,
                          const T defaultmin, const T defaultmax,
                          const T absmin,     const T absmax,
                          const char *name)
{
    // Sanity check incoming arguments and provide default name behavior
    using namespace std;
    if (!(absmin     <= defaultmin))
        throw invalid_argument("!(absmin <= defaultmin)");
    if (!(defaultmin <= defaultmax))
        throw invalid_argument("!(defaultmin <= defaultmax)");
    if (!(defaultmax <= absmax    ))
        throw invalid_argument("!(defaultmax <= absmax)");
    if (!name) name = "exprparse_range input";

    // Split s on a mandatory colon into whitespace-trimmed s_{min,max}
    const size_t colonpos = s.find_first_of(':');
    if (colonpos == string::npos) {
        throw invalid_argument(string(name)
            + " not in format \"low:high\", \"[low]:high\", or low:[high].");
    }
    string s_min(s, 0, colonpos);
    string s_max(s, colonpos + 1);
    boost::algorithm::trim(s_min);
    boost::algorithm::trim(s_max);

    // Parse recognized formats into valmin and valmax
    if (s_min.length() == 0 && s_max.length() == 0) {
        throw invalid_argument(string(name)
            + " not in format \"low:high\", \"[low]:high\", or low:[high].");
    } else if (s_min.length() == 0) {
        *valmin = defaultmin;
        *valmax = exprparse_impl<T>(s_max, name);
    } else if (s_max.length() == 0) {
        *valmin = exprparse_impl<T>(s_min, name);
        *valmax = defaultmax;
    } else {
        *valmin = exprparse_impl<T>(s_min, name);
        *valmax = exprparse_impl<T>(s_max, name);
    }

    // Ensure valmin <= valmax
    if (*valmin > *valmax) swap(*valmin, *valmax);

    // Validate range is within [absmin, absmax]
    if (*valmin < absmin || absmax < *valmax) {
        ostringstream m;
        m << name << " value [" << *valmin << ":" << *valmax
          << "] is outside valid range [" << absmin << ":" << absmax <<  "]";
        throw invalid_argument(m.str());
    }
}

} // end namespace detail

} // end namespace suzerain

#endif
