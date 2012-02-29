//--------------------------------------------------------------------------
//
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// timers.hpp: Conditional wrappers for performance timer routines
// $Id$

#ifndef TIMERS_HPP
#define TIMERS_HPP

#include <suzerain/suzerain-config.h>

#ifdef SUZERAIN_HAVE_GRVY
#include <grvy.h>
#define GRVY_TIMER_BEGIN(id)   grvy_timer_begin(id)
#define GRVY_TIMER_END(id)     grvy_timer_end(id)
#else
#define GRVY_TIMER_BEGIN(id)
#define GRVY_TIMER_END(id)
#endif

#endif /* TIMERS_HPP */
