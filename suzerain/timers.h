//--------------------------------------------------------------------------
//
// Copyright (C) 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
// timers.hpp: Conditional wrappers for performance timer routines
// $Id$

#ifndef SUZERAIN_TIMERS_H
#define SUZERAIN_TIMERS_H

#include <suzerain/suzerain-config.h>

/** @file
 * Wraps GRVY-based performance timing logic in a way permitting its consistent
 * use without the need to worry if GRVY is installed.
 */

#ifdef SUZERAIN_HAVE_GRVY
#include <grvy.h>
#define SUZERAIN_TIMER_BEGIN(id)   grvy_timer_begin(id)
#define SUZERAIN_TIMER_END(id)     grvy_timer_end(id)
#else
#define SUZERAIN_TIMER_BEGIN(id)
#define SUZERAIN_TIMER_END(id)
#endif

#endif /* SUZERAIN_TIMERS_H */
