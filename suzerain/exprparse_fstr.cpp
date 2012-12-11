//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// exprparse_fstr.cpp: arithmetic expression evaluation for float, string
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/exprparse.hpp>

#include "exprparse_impl.hpp"

// To mitigate long compilation times associated with Boost Spirit.Qi parsers,
// exprparse is broken into many separate compilation modules.

namespace suzerain {

void exprparse(const std::string& s, float& v, const char *name)
{
    v = detail::exprparse_impl<float>(s, name);
}

void exprparse_range(const std::string& s,
                     float *valmin, float *valmax,
                     const float defaultmin, const float defaultmax,
                     const float absmin, const float absmax,
                     const char *name)
{
    return detail::exprparse_range_impl(
            s, valmin, valmax, defaultmin, defaultmax, absmin, absmax, name);
}

} // namespace suzerain
