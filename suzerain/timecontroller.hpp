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
#include <suzerain/functional.hpp>
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

public:
    AbstractTimeController(FPT intial_t = 0,
                           FPT min_dt = 1e-8,
                           FPT max_dt = std::numeric_limits<FPT>::max());

    virtual ~AbstractTimeController();

    template<typename Callback>
    void addCallback(FPT every_dt,
                     Integer every_nt,
                     const Callback callback);

    FPT advanceTime(FPT final_t,
                    Integer final_nt = std::numeric_limits<Integer>::max());

protected:

    virtual FPT stepTime(FPT max_dt) const = 0;

private:

    struct Entry {
        FPT     every_dt;
        Integer every_nt;
        FPT     last_t;
        Integer last_nt;
        boost::signal<bool (FPT t, Integer nt),
                      suzerain::functional::all> signal;
    };

    const FPT min_dt;
    const FPT max_dt;
    FPT current_t;
    Integer current_nt;
    boost::ptr_list<Entry> entries;
};

template< typename FPT, typename Integer >
AbstractTimeController<FPT,Integer>::AbstractTimeController(FPT initial_t,
                                                            FPT min_dt,
                                                            FPT max_dt)
    : min_dt(min_dt), max_dt(max_dt),
      current_t(initial_t), current_nt(0),
      entries()
{
    // NOP
}

template< typename FPT, typename Integer >
AbstractTimeController<FPT,Integer>::~AbstractTimeController()
{
    // boost::ptr_list<Entry> semantics automatically destroys entries
}

template< typename FPT, typename Integer >
template< typename Callback >
void AbstractTimeController<FPT,Integer>::addCallback(FPT every_dt,
                                                      Integer every_nt,
                                                      Callback callback)
{
    // Linear search for existing Entry with same every_dt, every_nt
    typename boost::ptr_list<Entry>::iterator iter       = entries.begin();
    typename boost::ptr_list<Entry>::const_iterator end  = entries.end();
    for (;iter != end; ++iter) {
        if ((*iter).every_dt == every_dt && (*iter).every_nt == every_nt) {
            break;
        }
    }

    if (iter != end) {
        // Found a match.  Connect callback to the entry.
        (*iter).signal.connect(callback);
    } else {
        //  No exact match found.  Create a new Entry and track it.
        Entry *e    = new Entry;      // Allocate Entry on heap
        e->every_dt = every_dt;
        e->every_nt = every_nt;
        e->last_t   = current_t;
        e->last_nt  = current_nt;
        e->signal.connect(callback);
        entries.push_back(e);         // Transfer Entry ownership
    }
}

template< typename FPT, typename Integer >
FPT AbstractTimeController<FPT,Integer>::advanceTime(FPT final_t,
                                                     Integer final_nt)
{
    while (current_t < final_t) {
        FPT max_dt = final_t - current_t;
        for (typename boost::ptr_list<Entry>::iterator iter = entries.begin();
             iter != entries.end();
             ++iter) {
            // FIXME Continue implementing
        }
    }

    return current_t;
}

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMECONTROLLER_HPP
