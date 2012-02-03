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
// exprparse_dchar.cpp: arithmetic expression evaluation for double, char*
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/exprparse.hpp>
#include "exprparse_impl.hpp"

// To mitigate long compilation times associated with Boost Spirit.Qi parsers,
// exprparse is broken into many separate compilation modules.

void suzerain::exprparse(const char *s, double& v, const char *name)
{
    v = suzerain::detail::exprparse_impl<double>(s, name);
}
