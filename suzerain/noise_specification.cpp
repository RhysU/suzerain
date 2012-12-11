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
        real_t        fluct_percent,
        unsigned long fluct_seed,
        real_t        kxfrac_min,
        real_t        kxfrac_max,
        real_t        kzfrac_min,
        real_t        kzfrac_max)
    : percent   (fluct_percent)
    , seed      (fluct_seed)
    , kxfrac_min(kxfrac_min)
    , kxfrac_max(kxfrac_max)
    , kzfrac_min(kzfrac_min)
    , kzfrac_max(kzfrac_max)
{
}

} // namespace suzerain
