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
 * Provides higher-level control mechanisms based atop abstract
 * time integration schemes.
 */

namespace suzerain
{

namespace timestepper
{

/**
 * Provides an abstract base class for high-level time advancement control
 * including periodic callbacks.  Such callbacks could be used, for example, to
 * periodically write restart files, write statistics, and output status.
 * Subclasses must only implement the stepTime() method with the relevant
 * time integration logic to take a single step.
 *
 * @tparam FPT Floating point type used to track the current time.
 * @tparam Integer Integer type used to track the current time step.
 */
template<
    typename FPT     = double,
    typename Integer = unsigned long
>
class AbstractTimeController
{

public:

    /**
     * Construct an instance with the given parameters.
     *
     * @param initial_t Initial simulation time.
     * @param min_dt    Initial minimum acceptable time step.
     *                  See min_dt() for the associated semantics.
     * @param max_dt    Initial maximum acceptable time step.
     *                  See max_dt() for the associated semantics.
     */
    AbstractTimeController(
            FPT initial_t = 0,
            FPT min_dt = std::numeric_limits<FPT>::epsilon(),
            FPT max_dt = std::numeric_limits<FPT>::max());

    /**
     * Virtual destructor as appropriate for an abstract base class.
     **/
    virtual ~AbstractTimeController();

    //@{

    /**
     * Retrieve the minimum acceptable time step.  Physically
     * driven time steps below this value will cause advanceTime()
     * to abort.
     *
     * @return The current minimum acceptable time step.
     */
    FPT min_dt() const { return min_dt_; }

    /**
     * Set the minimum acceptable time step.
     *
     * @param new_min_dt Minimum acceptable time step to set.
     * @see min_dt() for more details.
     */
    void min_dt(FPT new_min_dt) { min_dt_ = new_min_dt; }

    //@}

    //@{

    /**
     * Retrieve the maximum acceptable time step.  The controller
     * will never take a time step in excess of this size.
     *
     * @return The current maximum acceptable time step.
     */
    FPT max_dt() const { return max_dt_; }

    /**
     * Set the maximum acceptable time step.
     *
     * @param new_max_dt Minimum acceptable time step to set.
     * @see max_dt() for more details.
     */
    void max_dt(FPT new_max_dt) { max_dt_ = new_max_dt; }

    //@}

    /**
     * Register a periodic callback to be used during advanceTime().  As the
     * controller marches the simulation, \c callback will be invoked \c
     * every_dt time steps or after \c every_t simulation time passes,
     * whichever comes first.
     *
     * Callbacks must be functions or functors taking two arguments and
     * returning a boolean value.  The first argument is of type \c FPT and
     * will contain the current simulation time.  The second argument is of
     * type \c Integer wand will contain the current simulation time step.  One
     * example signature would be
     * @code
     *   bool my_callback(FPT t, Integer nt);
     * @endcode
     * The callback should return \c true if the controller should continue
     * advancing.  If the callback returns \c false, the controller will
     * immediately stop advancing time.
     *
     * @param every_dt The maximum simulation time duration
     *                 between callbacks.
     * @param every_nt The simulation time step count
     *                 between callbacks.
     * @param callback The callback to invoke.
     *
     * @see <tt>std::numeric_limits<T>::max()</tt> if you want to
     *      specify either no criteria for \c every_dt or \c every_nt.
     * @see <a href="http://www.boost.org/doc/html/ref.html">Boost.Ref</a>
     *      if you need to provide a stateful or noncopyable functor
     *      to \c callback.
     */
    template<typename CallbackType>
    void addCallback(FPT every_dt,
                     Integer every_nt,
                     CallbackType callback);

    /**
     * Advance the simulation in time using stepTime() and perform required
     * callbacks along the way.  The simulation will advance until either
     * current_t() reaches \c final_t or current_nt() reaches \c final_nt,
     * whichever comes first.  If any physically-determined time step
     * is smaller than min_dt() or any registered callback returns \c false,
     * the controller will immediately stop.
     *
     * @param final_t  Maximum simulation time.
     *                 If you want to advance some relative amount \c dt,
     *                 try supplying <tt>current_t() + dt</tt>.
     * @param final_nt Maximum simulation time step.
     *                 If you want to advance some relative step count \c nt,
     *                 try supplying <tt>current_nt() + nt</tt>.
     *
     * @return True if the controller completed the request successfully.
     *         False if the controller aborted for some reason.
     */
    bool advanceTime(FPT final_t,
                     Integer final_nt = std::numeric_limits<Integer>::max());

    //@{

    /**
     * Retrieve the current simulation time.
     *
     * @return The current simulation time.
     */
    FPT current_t() const { return current_t_; }

    /**
     * Retrieve the current simulation time step.
     *
     * @return The current simulation discrete time step.
     */
    Integer current_nt() const { return current_nt_; }

    //@}

protected:

    /**
     * Advance simulation time by a single step of at most \c max_dt.
     * Subclasses implement this method to hook in relevant time integration
     * logic.  Physically-driven accuracy or stability requirements may limit
     * the actual step taken to be smaller than \c max_dt.
     *
     * @param max_dt Maximum step size to take.
     *
     * @return The actual step size taken by the time integration scheme.
     */
    virtual FPT stepTime(FPT max_dt) const = 0;

private:
    // Mark Entry as noncopyable to avoid accidental performance hits
    struct Entry : public boost::noncopyable {
        FPT     every_dt, next_t;
        Integer every_nt, next_nt;
        boost::function<bool (FPT t, Integer nt)> callback;
    };

    // ptr_container to explicitly manage noncopyable instances
    // ptr_vector to keep data as contiguous as possible in memory
    typedef boost::ptr_vector<Entry> EntryList;

    FPT min_dt_, max_dt_, current_t_;
    Integer current_nt_;
    EntryList entries_;

    template< typename T >
    static T add_and_coerce_overflow_to_max(T a, T b)
    {
        // Overflow occurs when max < a + b or equivalently max - b < a.
        // This overflow detection works if-and-only-if b is nonnegative!
        if (SUZERAIN_UNLIKELY(std::numeric_limits<T>::max() - b < a)) {
            return std::numeric_limits<T>::max();
        } else {
            return a + b;
        }
    }
};

template< typename FPT, typename Integer >
AbstractTimeController<FPT,Integer>::AbstractTimeController(FPT initial_t,
                                                            FPT min_dt,
                                                            FPT max_dt)
    : min_dt_(min_dt), max_dt_(max_dt), current_t_(initial_t),
      current_nt_(0),
      entries_(7)
{
    // NOP
}

template< typename FPT, typename Integer >
AbstractTimeController<FPT,Integer>::~AbstractTimeController()
{
    // boost::ptr_vector<Entry> semantics automatically destroys entries_
}

template< typename FPT, typename Integer >
template< typename CallbackType >
void AbstractTimeController<FPT,Integer>::addCallback(FPT every_dt,
                                                      Integer every_nt,
                                                      CallbackType callback)
{
    if (every_dt <= 0) throw std::invalid_argument("every_dt <= 0");
    if (every_nt <= 0) throw std::invalid_argument("every_nt <= 0");
    assert(min_dt_ <= max_dt_);

    Entry *e    = new Entry;      // Allocate Entry on heap
    e->every_dt = every_dt;
    e->every_nt = every_nt;
    e->next_t   = add_and_coerce_overflow_to_max(current_t_,  every_dt);
    e->next_nt  = add_and_coerce_overflow_to_max(current_nt_, every_nt);
    e->callback = callback;
    entries_.push_back(e);         // Transfer Entry memory ownership
}

template< typename FPT, typename Integer >
bool AbstractTimeController<FPT,Integer>::advanceTime(const FPT final_t,
                                                      const Integer final_nt)
{
    using std::min;

    // Maintain the next simulation time something interesting must happen
    FPT next_event_t = std::numeric_limits<FPT>::max();

    // Find the simulation time of the first callback
    for (typename EntryList::iterator iter = entries_.begin();
         iter != entries_.end();
         ++iter) {
        next_event_t = min(next_event_t, (*iter).next_t);
    }

    // Advance time until done or we abort for some reason
    while (current_t_ < final_t && current_nt_ < final_nt) {

        // Determine maximum possible step size allowed by all criteria
        const FPT possible_dt
            = min(max_dt_, min(final_t, next_event_t) - current_t_);
        assert(possible_dt > 0);

        // Take time step and then advance simulation time
        const FPT actual_dt = this->stepTime(possible_dt);
        assert(actual_dt <= possible_dt);
        current_t_  += actual_dt;
        current_nt_ += 1;

        // Check callbacks and determine next callback simulation time
        next_event_t = std::numeric_limits<FPT>::max();
        for (typename EntryList::iterator iter = entries_.begin();
             iter != entries_.end();
             ++iter) {

            // Callback required?  If so, do it.
            if (SUZERAIN_UNLIKELY(    current_t_  == (*iter).next_t
                                   || current_nt_ == (*iter).next_nt)) {

                // Update Entry with time of next required callback
                (*iter).next_t = add_and_coerce_overflow_to_max(
                        current_t_,  (*iter).every_dt);
                (*iter).next_nt = add_and_coerce_overflow_to_max(
                        current_nt_, (*iter).every_nt);

                // Perform callback
                if (!((*iter).callback(current_t_, current_nt_))) {
                    // Stop advancing if callback says so
                    return false;
                }
            }

            // Update next_event_t based on this callback's needs
            next_event_t = min(next_event_t, (*iter).next_t);
        }

        // Abort if the step size was too small, but only if driven by physics
        if (SUZERAIN_UNLIKELY(possible_dt >= min_dt_ && actual_dt < min_dt_)) {
            return false;
        }
    }

    // Successfully advanced up to provided criteria
    return true;
}

} // namespace timestepper

} // namespace suzerain

#endif // __SUZERAIN_TIMECONTROLLER_HPP
