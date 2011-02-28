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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/variance.hpp>

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
 * @tparam TimeType Floating point type used to track the current time.
 * @tparam StepType Integer type used to track the current time step.
 */
template<
    typename TimeType = double,
    typename StepType = std::size_t
>
class TimeController
{

public:

    /** Type used to express simulation time quantities. */
    typedef TimeType time_type;

    /** Type used to express discrete simulation step quantities. */
    typedef StepType step_type;

    /**
     * Construct a TimeController with the given parameters.
     *
     * The \c stepper argument is how one provides the time stepper that will
     * be driven by the constructed TimeController instance.  It must be a
     * function or functor compatible with <tt>boost::function<time_type
     * (time_type)></tt>.  When the stepper is invoked, the argument will be
     * the maximum possible time step that the stepper can take.  The
     * <tt>stepper</tt>'s return value must be the actual time step size taken.
     *
     * @param stepper   Time step logic wrapped by this controller.
     * @param initial_t Initial simulation time.
     * @param min_dt    Initial minimum acceptable time step.  Specifying
     *                  zero, the default, is equivalent to providing
     *                  <tt>std::numeric_limits<time_type>::epsilon()</tt>.
     *                  See min_dt() for the associated semantics.
     * @param max_dt    Initial maximum acceptable time step.  Specifying
     *                  zero, the default, is equivalent to providing
     *                  <tt>std::numeric_limits<time_type>::max()</tt>.
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
                   time_type initial_t = 0,
                   time_type min_dt = 0,
                   time_type max_dt = 0);

    /**
     * Virtual destructor since others may subclass this logic.
     */
    virtual ~TimeController() {};

    //@{

    /**
     * Retrieve the minimum acceptable time step.  Physically
     * driven time steps below this value will cause advance()
     * to abort.
     *
     * @return The current minimum acceptable time step.
     */
    time_type min_dt() const { return min_dt_; }

    /**
     * Set the minimum acceptable time step.
     *
     * @param new_min_dt Minimum acceptable time step to set.
     * @see min_dt() for more details.
     */
    void min_dt(time_type new_min_dt) { min_dt_ = new_min_dt; }

    //@}

    //@{

    /**
     * Retrieve the maximum acceptable time step.  The controller
     * will never take a time step in excess of this size.
     *
     * @return The current maximum acceptable time step.
     */
    time_type max_dt() const { return max_dt_; }

    /**
     * Set the maximum acceptable time step.
     *
     * @param new_max_dt Minimum acceptable time step to set.
     * @see max_dt() for more details.
     */
    void max_dt(time_type new_max_dt) { max_dt_ = new_max_dt; }

    //@}

    /**
     * Register a one-time callback to be invoked during advance()
     * when the simulation reaches \c what_t or \c what_nt, whichever
     * comes first.
     *
     * The argument \c callback must be a function or functor compatible with
     * <tt>boost::function<bool (time_type, step_type)></tt>.  When invoked,
     * the first argument will contain the current_t() and the second argument
     * will contain current_nt().  The callback must return \c true if the
     * controller should continue advancing.  If the callback returns \c false,
     * the controller will immediately stop advancing time.
     *
     * @param what_t  The simulation time to perform the callback.
     *                The provided value must be larger than current_t().
     * @param what_nt The simulation time step to perform the callback.
     *                The provided value must be larger than current_nt().
     * @param callback The callback to invoke.
     *
     * @see <tt>std::numeric_limits<T>::max()</tt> if you want to
     *      specify either no criteria for \c what_t or \c what_nt.
     * @see <a href="http://www.boost.org/doc/html/ref.html">Boost.Ref</a>
     *      if you need to provide a stateful or noncopyable functor
     *      as the \c callback argument.
     */
    template<typename CallbackType>
    void add_callback(time_type what_t,
                      step_type what_nt,
                      CallbackType callback);

    /**
     * Register a periodic callback to be invoked during advance().  As the
     * controller marches the simulation, \c callback will be invoked \c
     * every_dt time steps or after \c every_t simulation time passes,
     * whichever comes first.
     *
     * @param every_dt The maximum simulation time duration between callbacks.
     * @param every_nt The simulation time step count between callbacks.
     * @param callback The callback to invoke.  See add_callback() for
     *                 a discussion of this argument's semantics.
     *
     * @see <tt>std::numeric_limits<T>::max()</tt> if you want to
     *      specify either no criteria for \c every_dt or \c every_nt.
     * @see <a href="http://www.boost.org/doc/html/ref.html">Boost.Ref</a>
     *      if you need to provide a stateful or noncopyable functor
     *      as the \c callback argument.
     */
    template<typename CallbackType>
    void add_periodic_callback(time_type every_dt,
                               step_type every_nt,
                               CallbackType callback);

    /**
     * Advance the simulation in time using the time stepper set at
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
     * @see The method step() for an easy way to specify a fixed number of
     *      time steps.
     */
    bool advance(time_type final_t,
                 step_type final_nt = std::numeric_limits<step_type>::max());

    /**
     * Advance the simulation in time using the time stepper set at
     * construction time and perform required callbacks along the way.  The
     * simulation will advance until \c count_nt time steps have been
     * completed.  If any physically-determined time step is smaller than
     * min_dt() or any registered callback returns \c false, the controller
     * will immediately stop.
     *
     * @param count_nt Number of time steps to take.
     *
     * @return True if the controller completed the request successfully.
     *         False if the controller aborted for some reason.
     * @see The method advance() for a richer interface for controlling
     *      time advancement.
     */
    bool step(step_type count_nt = 1)
    {
        return advance(std::numeric_limits<time_type>::max(),
                       add_and_coerce_overflow_to_max(count_nt, current_nt()));
    }

    //@{

    /**
     * Retrieve the current simulation time.
     *
     * @return The current simulation time.
     */
    time_type current_t() const { return current_t_; }

    /**
     * Retrieve the current simulation time step.
     *
     * @return The current simulation discrete time step.
     */
    step_type current_nt() const {
        return boost::accumulators::extract::count(dt_stats);
    }

    //@}

    //@{

    /**
     * Retrieve the minimum time step taken during time advancement
     *
     * @return The minimum time step taken during time advancement
     */
    time_type taken_min() const {
        return boost::accumulators::extract::min(dt_stats);
    }

    /**
     * Retrieve the mean time step taken during time advancement
     *
     * @return The mean time step taken during time advancement
     */
    time_type taken_mean() const {
        return boost::accumulators::extract::mean(dt_stats);
    }

    /**
     * Retrieve the maximum time step taken during time advancement
     *
     * @return The maximum time step taken during time advancement
     */
    time_type taken_max() const {
        return boost::accumulators::extract::max(dt_stats);
    }

    /**
     * Retrieve the standard deviation of the time steps taken during
     * time advancement.
     *
     * @return the standard deviation of the time steps taken during
     * time advancement.
     */
    time_type taken_stddev() const {
        return std::sqrt(boost::accumulators::extract::variance(dt_stats));
    }

    //@}

private:

    // Mark Entry as noncopyable to avoid accidental performance hits
    struct Entry : public boost::noncopyable {
        bool periodic;
        time_type every_dt, next_t;
        step_type every_nt, next_nt;
        boost::function<bool (time_type t, step_type nt)> callback;
    };

    // ptr_container to explicitly manage noncopyable instances with minimal
    // overhead and to automatically free contained resources on destruction.
    // ptr_vector to keep data as contiguous as possible in memory.
    typedef boost::ptr_vector<Entry> EntryList;

    typename boost::function<time_type (time_type)> stepper_;
    time_type min_dt_, max_dt_, current_t_;
    EntryList entries_;

    // Maintain running statistics on the actual time step sizes taken
    // Also stores the current time step counter as a convenient side effect
    boost::accumulators::accumulator_set<
            time_type,
            boost::accumulators::features<
                boost::accumulators::tag::count,
                boost::accumulators::tag::min,
                boost::accumulators::tag::mean,
                boost::accumulators::tag::max,
                boost::accumulators::tag::lazy_variance
            >
        > dt_stats;

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

template< typename TimeType, typename StepType >
template< typename StepperType >
TimeController<TimeType,StepType>::TimeController(StepperType stepper,
                                                  time_type initial_t,
                                                  time_type min_dt,
                                                  time_type max_dt)
    : stepper_(stepper),
      min_dt_(min_dt != 0 ? min_dt : std::numeric_limits<time_type>::epsilon()),
      max_dt_(max_dt != 0 ? max_dt : std::numeric_limits<time_type>::max()),
      current_t_(initial_t),
      entries_(7)
{
    // NOP
}

template< typename TimeType, typename StepType >
template< typename CallbackType >
void TimeController<TimeType,StepType>::add_callback(time_type what_t,
                                                     step_type what_nt,
                                                     CallbackType callback)
{
    // Callbacks established in the past are likely usage errors (<).  The
    // advance() method does not perform callbacks prior to advancing
    // simulation time, so we cannot handle "immediate" callback requests (==).
    // Combined, these provide the (<=) constraints just below.
    if (what_t <= current_t_) {
        throw std::invalid_argument("what_t <= current_t()");
    }
    if (what_nt <= current_nt()) {
        throw std::invalid_argument("what_nt <= current_nt()");
    }

    Entry *e    = new Entry;      // Allocate Entry on heap
    e->periodic = false;
    e->every_dt = 0;
    e->every_nt = 0;
    e->next_t   = what_t;
    e->next_nt  = what_nt;
    e->callback = callback;
    entries_.push_back(e);        // Transfer Entry memory ownership
}

template< typename TimeType, typename StepType >
template< typename CallbackType >
void TimeController<TimeType,StepType>::add_periodic_callback(
        time_type every_dt,
        step_type every_nt,
        CallbackType callback)
{
    if (every_dt <= 0) throw std::invalid_argument("every_dt <= 0");
    if (every_nt <= 0) throw std::invalid_argument("every_nt <= 0");

    Entry *e    = new Entry;      // Allocate Entry on heap
    e->periodic = true;
    e->every_dt = every_dt;
    e->every_nt = every_nt;
    e->next_t   = add_and_coerce_overflow_to_max(current_t_,   every_dt);
    e->next_nt  = add_and_coerce_overflow_to_max(current_nt(), every_nt);
    e->callback = callback;
    entries_.push_back(e);        // Transfer Entry memory ownership
}

template< typename TimeType, typename StepType >
bool TimeController<TimeType,StepType>::advance(const time_type final_t,
                                                const step_type final_nt)
{
    assert(min_dt_ <= max_dt_);
    using std::min;

    // Maintain the next simulation time something interesting must happen
    time_type next_event_t = std::numeric_limits<time_type>::max();

    // Find the simulation time of the first callback
    for (typename EntryList::iterator iter = entries_.begin();
         iter != entries_.end();
         ++iter) {
        next_event_t = min(next_event_t, (*iter).next_t);
    }

    // Advance time until done or we abort for some reason
    while (current_t_ < final_t && current_nt() < final_nt) {

        // Determine maximum possible step size allowed by all criteria
        const time_type possible_dt
            = min(max_dt_, min(final_t, next_event_t) - current_t_);
        assert(possible_dt > 0);

        // Take step, accumulate statistics, and record new simulation time
        const time_type actual_dt = stepper_(possible_dt);
        assert(actual_dt <= possible_dt);
        dt_stats(actual_dt);
        current_t_ += actual_dt;

        // Check callbacks and determine next callback simulation time
        next_event_t = std::numeric_limits<time_type>::max();
        typename EntryList::iterator iter = entries_.begin();
        while (iter != entries_.end()) {

            // Callback required?
            if (SUZERAIN_UNLIKELY(    current_t_   == (*iter).next_t
                                   || current_nt() == (*iter).next_nt)) {

                // Perform required callback
                // Must perform state updates prior to any possible abort
                const bool keep_advancing
                    = (*iter).callback(current_t_, current_nt());

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
                        current_t_,   (*iter).every_dt);
                (*iter).next_nt = add_and_coerce_overflow_to_max(
                        current_nt(), (*iter).every_nt);

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
