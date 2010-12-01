/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * timecontroller.hpp: higher-level time advance management
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_TIMECONTROLLER_HPP
#define __SUZERAIN_TIMECONTROLLER_HPP

#include <suzerain/common.hpp>
#include <suzerain/timestepper.hpp>
#include <suzerain/traits.hpp>

/** @file
 * Provides time integration schemes.
 */

namespace suzerain
{

namespace timestepper
{

/**
 * Provides an abstract base class for higher-level time advance
 * management including periodic restart and statistics writing events.
 */
template<
    typename FPT     = double,
    typename Integer = unsigned long
>
class AbstractTimeController
{
private:
    struct Entry {
        FPT every_t;
        FPT every_nt;
    };

public:
    virtual ~AbstractTimeController() {}

protected:
    virtual FPT step(const FPT max_delta_t) const = 0;
};

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMECONTROLLER_HPP
