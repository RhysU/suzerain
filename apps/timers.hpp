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
#define GRVY_TIMER_FINALIZE()  grvy_timer_finalize()
#define GRVY_TIMER_INIT(id)    grvy_timer_init(id)
#define GRVY_TIMER_RESET()     grvy_timer_reset()
#define GRVY_TIMER_SUMMARIZE() grvy_timer_summarize()
#else
#define GRVY_TIMER_BEGIN(id)
#define GRVY_TIMER_END(id)
#define GRVY_TIMER_FINALIZE()
#define GRVY_TIMER_INIT(id)
#define GRVY_TIMER_RESET()
#define GRVY_TIMER_SUMMARIZE()
#endif

#endif /* TIMERS_HPP */
