//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// noise_specification.cpp: classes handling specifying additive noise
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/noise_specification.hpp>

namespace suzerain {

noise_specification::noise_specification(
        real_t fluct_percent,
        unsigned long fluct_seed)
    : percent(percent)
    , kxfrac_min(0)
    , kxfrac_max(1)
    , kzfrac_min(0)
    , kzfrac_max(1)
    , seed(seed)
{
}

} // namespace suzerain
