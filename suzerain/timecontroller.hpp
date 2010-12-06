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

/** @file
 * Provides higher-level control mechanisms based atop abstract
 * time integration schemes.
 */

namespace suzerain
{

namespace timestepper
{

/**
 * Provides high-level time advancement control, including periodic callbacks,
 * built atop a simple time stepping interface.  Such callbacks could be
 * used, for example, to periodically write restart files, write statistics,
 * and output status.
 *
 * @tparam FPT Floating point type used to track the current time.
 * @tparam Integer Integer type used to track the current time step.
 */
template<
    typename FPT     = double,
    typename Integer = unsigned long
>
class TimeController
{

public:

    /**
     * Construct a TimeController with the given parameters.
     *
     * The \c stepper argument is how one provides the time stepper that will
     * be driven by the constructed TimeController instance.  It must be a
     * function or functor compatible with <tt>boost::function<FPT (FPT)></tt>.
     * When the stepper is invoked, the argument will be the maximum possible
     * time step that the stepper can take.  The <tt>stepper</tt>'s return
     * value must be the actual time step size taken.
     *
     * @param stepper   Time step logic wrapped by this controller.
     * @param initial_t Initial simulation time.
     * @param min_dt    Initial minimum acceptable time step.
     *                  See min_dt() for the associated semantics.
     * @param max_dt    Initial maximum acceptable time step.
     *                  See max_dt() for the associated semantics.
     *
     * @see <a href="http://www.boost.org/doc/html/function.html">
     *      Boost.Function</a> for more details on
     *      <tt>boost::function</tt>.
     * @see <a href="http://www.boost.org/doc/html/ref.html">Boost.Ref</a>
     *      if you need to provide a stateful or noncopyable functor
     *      as the \c stepper argument.
     */
    template<typename StepperType>
    TimeController(StepperType stepper,
                   FPT initial_t = 0,
                   FPT min_dt = std::numeric_limits<FPT>::epsilon(),
                   FPT max_dt = std::numeric_limits<FPT>::max());

    //@{

    /**
     * Retrieve the minimum acceptable time step.  Physically
     * driven time steps below this value will cause advance()
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
     * Register a one-time callback to be invoked during advance()
     * when the simulation reaches \c what_t or \c what_nt, whichever
     * comes first.
     *
     * The argument \c callback must be a function or functor compatible with
     * <tt>boost::function<bool (FPT, Integer)></tt>.  When invoked, the first
     * argument will contain the current_t() and the second argument will
     * contain current_nt().  The callback must return \c true if the
     * controller should continue advancing.  If the callback returns \c false,
     * the controller will immediately stop advancing time.
     *
     * @param what_t  The simulation time to perform the callback.
     * @param what_nt The simulation time step to perform the callback.
     * @param callback The callback to invoke.
     *
     * @see <tt>std::numeric_limits<T>::max()</tt> if you want to
     *      specify either no criteria for \c what_t or \c what_nt.
     * @see <a href="http://www.boost.org/doc/html/ref.html">Boost.Ref</a>
     *      if you need to provide a stateful or noncopyable functor
     *      as the \c callback argument.
     */
    template<typename CallbackType>
    void add_callback(FPT what_t,
                      Integer what_nt,
                      CallbackType callback);

    /**
     * Register a periodic callback to be invoked during advance().  As the
     * controller marches the simulation, \c callback will be invoked \c
     * every_dt time steps or after \c every_t simulation time passes,
     * whichever comes first.
     *
     * @param every_dt The maximum simulation time duration
     *                 between callbacks.
     * @param every_nt The simulation time step count
     *                 between callbacks.
     * @param callback The callback to invoke.  See add_callback() for
     *                 a discussion of this argument's semantics.
     *
     * @see <tt>std::numeric_limits<T>::max()</tt> if you want to
     *      specify either no criteria for \c every_dt or \c every_nt.
     */
    template<typename CallbackType>
    void add_periodic_callback(FPT every_dt,
                               Integer every_nt,
                               CallbackType callback);

    /**
     * Advance the simulation in time using the time stepper set and
     * construction time and perform required callbacks along the way.  The
     * simulation will advance until either current_t() reaches \c final_t or
     * current_nt() reaches \c final_nt, whichever comes first.  If any
     * physically-determined time step is smaller than min_dt() or any
     * registered callback returns \c false, the controller will immediately
     * stop.
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
    bool advance(FPT final_t,
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

private:

    // Mark Entry as noncopyable to avoid accidental performance hits
    struct Entry : public boost::noncopyable {
        bool periodic;
        FPT     every_dt, next_t;
        Integer every_nt, next_nt;
        boost::function<bool (FPT t, Integer nt)> callback;
    };

    // ptr_container to explicitly manage noncopyable instances with minimal
    // overhead and to automatically free contained resources on destruction.
    // ptr_vector to keep data as contiguous as possible in memory.
    typedef boost::ptr_vector<Entry> EntryList;

    typename boost::function<FPT (FPT)> stepper_;
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
template< typename StepperType >
TimeController<FPT,Integer>::TimeController(StepperType stepper,
                                            FPT initial_t,
                                            FPT min_dt,
                                            FPT max_dt)
    : stepper_(stepper),
      min_dt_(min_dt), max_dt_(max_dt), current_t_(initial_t),
      current_nt_(0),
      entries_(7)
{
    // NOP
}

template< typename FPT, typename Integer >
template< typename CallbackType >
void TimeController<FPT,Integer>::add_callback(FPT what_t,
                                               Integer what_nt,
                                               CallbackType callback)
{
    Entry *e    = new Entry;      // Allocate Entry on heap
    e->periodic = false;
    e->every_dt = 0;
    e->every_nt = 0;
    e->next_t   = what_t;
    e->next_nt  = what_nt;
    e->callback = callback;
    entries_.push_back(e);        // Transfer Entry memory ownership
}

template< typename FPT, typename Integer >
template< typename CallbackType >
void TimeController<FPT,Integer>::add_periodic_callback(FPT every_dt,
                                                        Integer every_nt,
                                                        CallbackType callback)
{
    if (every_dt <= 0) throw std::invalid_argument("every_dt <= 0");
    if (every_nt <= 0) throw std::invalid_argument("every_nt <= 0");

    Entry *e    = new Entry;      // Allocate Entry on heap
    e->periodic = true;
    e->every_dt = every_dt;
    e->every_nt = every_nt;
    e->next_t   = add_and_coerce_overflow_to_max(current_t_,  every_dt);
    e->next_nt  = add_and_coerce_overflow_to_max(current_nt_, every_nt);
    e->callback = callback;
    entries_.push_back(e);        // Transfer Entry memory ownership
}

template< typename FPT, typename Integer >
bool TimeController<FPT,Integer>::advance(const FPT final_t,
                                          const Integer final_nt)
{
    assert(min_dt_ <= max_dt_);
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

        // Take time step and record new simulation time
        const FPT actual_dt = stepper_(possible_dt);
        assert(actual_dt <= possible_dt);
        current_t_  += actual_dt;
        current_nt_ += 1;

        // Check callbacks and determine next callback simulation time
        next_event_t = std::numeric_limits<FPT>::max();
        typename EntryList::iterator iter = entries_.begin();
        while (iter != entries_.end()) {

            // Callback required?
            if (SUZERAIN_UNLIKELY(    current_t_  == (*iter).next_t
                                   || current_nt_ == (*iter).next_nt)) {

                // Perform required callback
                // Must perform state updates prior to any possible abort
                const bool keep_advancing
                    = (*iter).callback(current_t_, current_nt_);

                // Remove single-shot Entry from further consideration
                if (SUZERAIN_UNLIKELY(!((*iter).periodic))) {
                    iter = entries_.erase(iter);
                    if (SUZERAIN_UNLIKELY(!keep_advancing)) {
                        return false;
                    }
                    continue;
                }

                // Update periodic Entry with time of next required callback
                (*iter).next_t = add_and_coerce_overflow_to_max(
                        current_t_,  (*iter).every_dt);
                (*iter).next_nt = add_and_coerce_overflow_to_max(
                        current_nt_, (*iter).every_nt);

                if (SUZERAIN_UNLIKELY(!keep_advancing)) {
                    return false;
                }
            }

            // Update next_event_t based on this Entry's needs
            next_event_t = min(next_event_t, (*iter).next_t);

            // Next!
            ++iter;
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
