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
// ndx.cpp: Symbolic indices, names, and descriptions for Navier--Stokes
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/ndx.hpp>

namespace suzerain {

namespace ndx {

const boost::array<const char *, rho + 1u> identifier = {{
    "rho_E", "rho_u", "rho_v", "rho_w", "rho"
}};

const boost::array<const char *, rho + 1u> description = {{
    "total energy",
    "streamwise momentum",
    "wall-normal momentum",
    "spanwise momentum",
    "density"
}};

} // namespace ndx

} // namespace suzerain
