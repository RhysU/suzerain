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
    AbstractTimeController(FPT initial_t = 0,
                           FPT min_dt = 1e-8,
                           FPT max_dt = std::numeric_limits<FPT>::max());

    virtual ~AbstractTimeController();

    FPT current_t() const { return current_t_; }

    Integer current_nt() const { return current_nt_; }

    FPT min_dt() const { return min_dt_; }

    void min_dt(FPT new_min_dt) { min_dt_ = new_min_dt; }

    FPT max_dt() const { return max_dt_; }

    void max_dt(FPT new_max_dt) { max_dt_ = new_max_dt; }

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
        FPT     next_t;
        Integer next_nt;
        boost::signal<bool (FPT t, Integer nt),
                      suzerain::functional::all> signal; // noncopyable
    };

    typedef boost::ptr_vector<Entry> EntryList;

    FPT min_dt_;
    FPT max_dt_;
    FPT current_t_;
    Integer current_nt_;
    EntryList entries_;
};

template< typename FPT, typename Integer >
AbstractTimeController<FPT,Integer>::AbstractTimeController(FPT initial_t,
                                                            FPT min_dt,
                                                            FPT max_dt)
    : min_dt_(min_dt), max_dt_(max_dt),
      current_t_(initial_t), current_nt_(0),
      entries_(5)
{
    // NOP
}

template< typename FPT, typename Integer >
AbstractTimeController<FPT,Integer>::~AbstractTimeController()
{
    // boost::ptr_vector<Entry> semantics automatically destroys entries_
}

template< typename FPT, typename Integer >
template< typename Callback >
void AbstractTimeController<FPT,Integer>::addCallback(FPT every_dt,
                                                      Integer every_nt,
                                                      Callback callback)
{
    if (every_dt <= 0) throw std::invalid_argument("every_dt <= 0");
    if (every_nt <= 0) throw std::invalid_argument("every_nt <= 0");

    // Linear search for existing Entry with same every_dt, every_nt
    typename EntryList::iterator iter      = entries_.begin();
    typename EntryList::const_iterator end = entries_.end();
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
        e->next_t   = current_t_ + every_dt;
        e->next_nt  = current_nt_ + every_nt;
        e->signal.connect(callback);
        entries_.push_back(e);         // Transfer Entry ownership
    }
}

template< typename FPT, typename Integer >
FPT AbstractTimeController<FPT,Integer>::advanceTime(const FPT final_t,
                                                     const Integer final_nt)
{
    while (current_t_ < final_t && current_nt_ < final_nt) {

        // Determine maximum possible step size allowed by all criteria
        FPT possible_dt = std::min(max_dt_, final_t - current_t_);
        for (typename EntryList::iterator iter = entries_.begin();
             iter != entries_.end();
             ++iter) {
            possible_dt = std::min(possible_dt, (*iter).next_t - current_t_);
        }

        // Take time step and advance simulation time
        const FPT actual_dt = this->stepTime(possible_dt);
        current_t_  += actual_dt;
        current_nt_ += 1;

        // Check callbacks
        for (typename EntryList::iterator iter = entries_.begin();
             iter != entries_.end();
             ++iter) {

            // Callback required?
            if (SUZERAIN_UNLIKELY(    current_t_ == (*iter).next_t
                                   || current_nt_ == (*iter).next_nt)) {

                // Update Entry with time of next required callback
                (*iter).next_t  = current_t_  + (*iter).every_dt;
                (*iter).next_nt = current_nt_ + (*iter).every_nt;

                // Perform callback with abort if callback returns false
                if (!((*iter).signal(current_t_, current_nt_))) {
                    goto abort;
                }
            }
        }

        // Abort if the step size was too small, but only if driven by physics
        if (SUZERAIN_UNLIKELY(possible_dt >= min_dt_ && actual_dt < min_dt_)) {
            goto abort;
        }
    }

abort:

    return current_t_;
}

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMECONTROLLER_HPP
