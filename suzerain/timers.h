//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_TIMERS_H
#define SUZERAIN_TIMERS_H

/** @file
 * Wraps GRVY-based performance timers to permit their consistent
 * use without the need to worry if GRVY is available.
 */

#include <suzerain/suzerain-config.h>

#ifdef SUZERAIN_HAVE_GRVY
#include <grvy.h>
#define SUZERAIN_TIMER_BEGIN(id)   grvy_timer_begin(id)
#define SUZERAIN_TIMER_END(id)     grvy_timer_end(id)
#else
#define SUZERAIN_TIMER_BEGIN(id)
#define SUZERAIN_TIMER_END(id)
#endif

#ifdef __cplusplus

/**
 * Invoke <tt>SUZERAIN_TIMER_BEGIN(id)</tt> and
 * <tt>SUZERAIN_TIMER_BEGIN(id)</tt> using a scope guard.  The C-style string
 * identifier <tt>id</tt> must have a lifetime longer than the scope guard,
 * otherwise the behavior is undefined.
 */
#define SUZERAIN_TIMER_SCOPED(id)                                   \
    SUZERAIN_TIMER_SCOPED_HELPER1(id,__LINE__)

#define SUZERAIN_TIMER_SCOPED_HELPER1(id,line)                      \
    SUZERAIN_TIMER_SCOPED_HELPER2(id,line)

#define SUZERAIN_TIMER_SCOPED_HELPER2(id,line)             \
class SuzerainTimerScoped_##line {                         \
    const char *id_;                                       \
public:                                                    \
    SuzerainTimerScoped_##line(const char *id_) : id_(id_) \
        { SUZERAIN_TIMER_BEGIN(id_); }                     \
    ~SuzerainTimerScoped_##line()                          \
        { SUZERAIN_TIMER_END  (id_); }                     \
} suzeraintimerscoped_##line(id)

#endif  // #ifdef __cplusplus

#endif /* SUZERAIN_TIMERS_H */
