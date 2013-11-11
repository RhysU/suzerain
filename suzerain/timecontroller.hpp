//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_TIMECONTROLLER_HPP
#define SUZERAIN_TIMECONTROLLER_HPP

/** @file
 * Provides high-level control mechanisms based atop step-based
 * time integration schemes.
 */

#include <suzerain/common.hpp>
#include <suzerain/timers.h>

namespace suzerain
{

/**
 * Provides high-level time advancement control, including periodic callbacks,
 * built atop a simple time stepping interface.  Such callbacks could be
 * used, for example, to periodically write restart files, write statistics,
 * and output status.
 *
 * @tparam TimeType Floating point type used to track the current time.
 * @tparam StepType Integer type used to track the current time step.
 *                  Using types with greater ranges than <tt>std::size_t</tt>
 *                  is not recommended.
 * @tparam StopType Type used to represent success or failure within callbacks
 *                  and time advancement routines.  Provided so that this
 *                  class may use C++'s obvious <tt>bool</tt> choice but may
 *                  be easily wrapped in a C-based API using <tt>int</tt>s.
 */
template<
    typename TimeType = real_t,
    typename StepType = std::size_t,
    typename StopType = bool
>
class timecontroller
{

public:

    //@{

    /** Type used to express simulation time quantities. */
    typedef TimeType time_type;

    /** Type used to express discrete simulation step quantities. */
    typedef StepType step_type;

    /** Type used by callbacks to tell the controller to stop advancing. */
    typedef StopType stop_type;

    /**
     * Obtain a simulation time equivalent to running "forever".  This is
     * useful for specifying callbacks with indefinite time periods (that is,
     * callbacks never performed based on the current simulation time).  It is
     * equivalent to using <tt>std::numeric_limits<TimeType>::max()</tt>.
     *
     * @return A simulation time equivalent to running forever.
     */
    static TimeType forever_t() {
        return std::numeric_limits<TimeType>::max();
    }

    /**
     * Obtain a simulation time step count equivalent to running forever.
     * This is useful for specifying callbacks with indefinite time step
     * periods (that is, callbacks never performed based on the current
     * simulation timestep).  It is equivalent to using
     * <tt>std::numeric_limits<StepType>::max()</tt>.
     *
     * @return A simulation time step count equivalent to running forever.
     */
    static StepType forever_nt() {
        return std::numeric_limits<StepType>::max();
    }

    //@}

    /**
     * Construct a timecontroller with the given parameters.
     *
     * The \c stepper argument is how one provides the time stepper that will
     * be driven by the constructed timecontroller instance.  It must be a
     * function or functor compatible with <tt>boost::function<TimeType
     * (TimeType)></tt>.  When the stepper is invoked, the argument will be
     * the maximum possible time step that the stepper can take.  The
     * <tt>stepper</tt>'s return value must be the actual time step size taken.
     *
     * @param stepper   Time step logic wrapped by this controller.
     * @param initial_t Initial simulation time.
     * @param min_dt    Initial minimum acceptable time step.  Specifying
     *                  zero, the default, is equivalent to providing
     *                  <tt>std::numeric_limits<TimeType>::epsilon()</tt>.
     *                  See min_dt() for the associated semantics.
     * @param max_dt    Initial maximum acceptable time step.  Specifying
     *                  zero, the default, is equivalent to providing
     *                  <tt>std::numeric_limits<TimeType>::max()</tt>.
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
    timecontroller(StepperType stepper,
                   TimeType initial_t = 0,
                   TimeType min_dt = 0,
                   TimeType max_dt = 0);

    /**
     * Virtual destructor since others may subclass this logic.
     */
    virtual ~timecontroller() {};

    //@{

    /**
     * Retrieve the minimum acceptable time step.  Physically
     * driven time steps below this value will cause advance()
     * to abort.
     *
     * @return The current minimum acceptable time step.
     */
    TimeType min_dt() const { return min_dt_; }

    /**
     * Set the minimum acceptable time step.
     *
     * @param new_min_dt Minimum acceptable time step to set.
     * @see min_dt() for more details.
     */
    void min_dt(TimeType new_min_dt) { min_dt_ = new_min_dt; }

    //@}

    //@{

    /**
     * Retrieve the maximum acceptable time step.  The controller
     * will never take a time step in excess of this size.
     *
     * @return The current maximum acceptable time step.
     */
    TimeType max_dt() const { return max_dt_; }

    /**
     * Set the maximum acceptable time step.
     *
     * @param new_max_dt Minimum acceptable time step to set.
     * @see max_dt() for more details.
     */
    void max_dt(TimeType new_max_dt) { max_dt_ = new_max_dt; }

    //@}

    /**
     * Register a one-time callback to be invoked during advance()
     * when the simulation reaches \c what_t or \c what_nt, whichever
     * comes first.
     *
     * The argument \c callback must be a function or functor compatible with
     * <tt>boost::function<StopType (TimeType, StepType)></tt>.  When
     * invoked, the first argument will contain the current_t() and the second
     * argument will contain current_nt().  The callback must return \c true if
     * the controller should continue advancing.  If the callback returns \c
     * false, the controller will immediately stop advancing time.
     *
     * @param what_t  The simulation time to perform the callback.
     *                The provided value must be larger than current_t().
     *                Use forever_t() if you want no callbacks based on
     *                the current simulation time.
     * @param what_nt The simulation time step to perform the callback.
     *                The provided value must be larger than current_nt().
     *                Use forever_nt() if you want no callbacks based on
     *                the current time step.
     * @param callback The callback to invoke.
     *
     * Specifying both <tt>what_t == forever_t()</tt> and <tt>what_nt ==
     * forever_nt()</tt> causes the callback to be discarded.
     *
     * @see <a href="http://www.boost.org/doc/html/ref.html">Boost.Ref</a>
     *      if you need to provide a stateful or noncopyable functor
     *      as the \c callback argument.
     */
    template<typename CallbackType>
    void add_callback(TimeType what_t,
                      StepType what_nt,
                      CallbackType callback);

    /**
     * Register a periodic callback to be invoked during advance().  As the
     * controller marches the simulation, \c callback will be invoked \c
     * every_dt time steps or after \c every_t simulation time passes,
     * whichever comes first.
     *
     * @param every_dt The maximum simulation time duration between callbacks.
     *                 Use forever_t() if you want no callbacks based on
     *                 the current simulation time.
     * @param every_nt The simulation time step count between callbacks.
     *                 Use forever_nt() if you want no callbacks based on
     *                 the current time step.
     * @param callback The callback to invoke.  See add_callback() for
     *                 a discussion of this argument's semantics.
     *
     * Specifying both <tt>every_dt == forever_t()</tt> and <tt>every_nt ==
     * forever_nt()</tt> causes the callback to be discarded.
     *
     * @see <a href="http://www.boost.org/doc/html/ref.html">Boost.Ref</a>
     *      if you need to provide a stateful or noncopyable functor
     *      as the \c callback argument.
     */
    template<typename CallbackType>
    void add_periodic_callback(TimeType every_dt,
                               StepType every_nt,
                               CallbackType callback);

    /**
     * Register a periodic callback to be invoked during advance().  As the
     * controller marches the simulation, \c callback will be invoked \c
     * every_dt time steps or after \c every_t simulation time passes,
     * whichever comes first.  The first such invocation will occur no
     * more than \c first_dt time units after this registration call.
     *
     * @param first_dt The maximum simulation time duration before the
     *                 first callback.  Use forever_t() if you want no
     *                 callbacks based on the current simulation time.
     * @param every_dt The maximum simulation time duration between callbacks.
     *                 Use forever_t() if you want no callbacks based on
     *                 the current simulation time.
     * @param every_nt The simulation time step count between callbacks.
     *                 Use forever_nt() if you want no callbacks based on
     *                 the current time step.
     * @param callback The callback to invoke.  See add_callback() for
     *                 a discussion of this argument's semantics.
     *
     * Specifying <tt>every_dt == forever_t()</tt> and <tt>every_nt ==
     * forever_nt()</tt> causes the callback to be discarded.  If a
     * once-only callback is desired, instead use \ref add_callback.
     *
     * @see <a href="http://www.boost.org/doc/html/ref.html">Boost.Ref</a>
     *      if you need to provide a stateful or noncopyable functor
     *      as the \c callback argument.
     */
    template<typename CallbackType>
    void add_periodic_callback(TimeType first_dt,
                               TimeType every_dt,
                               StepType every_nt,
                               CallbackType callback);

    /**
     * Advance the simulation in time using the time stepper set at
     * construction time and perform required callbacks along the way.  The
     * simulation will advance until either current_t() reaches \c final_t or
     * current_nt() reaches \c final_nt, whichever comes first.  If any
     * physically-determined time step is smaller than min_dt() or any
     * registered callback returns \c false, the controller will immediately
     * stop.  The former condition can be diagnosed by comparing current_dt()
     * to min_dt().  The controller will also stop if it encounters a step size
     * for which <tt>isfinite</tt> is \c false.
     *
     * @param final_t  Maximum simulation time.
     *                 If you want to advance some relative amount \c dt,
     *                 try supplying <tt>current_t() + dt</tt>.  For unlimited
     *                 simulation time, use forever_t().
     * @param final_nt Maximum simulation time step.
     *                 If you want to advance some relative step count \c nt,
     *                 try supplying <tt>current_nt() + nt</tt>.  For unlimited
     *                 time step counts, use forever_nt().
     *
     * @return True if the controller completed the request successfully.
     *         False if the controller aborted for some reason.
     * @see The method step() for an easy way to specify a fixed number of
     *      time steps.
     */
    StopType advance(
            TimeType final_t  = std::numeric_limits<TimeType>::max(),
            StepType final_nt = std::numeric_limits<StepType>::max());

    /**
     * Advance the simulation in time using the time stepper set at
     * construction time and perform required callbacks along the way.  The
     * simulation will advance until \c count_nt time steps have been
     * completed.  If any physically-determined time step is smaller than
     * min_dt() or any registered callback returns \c false, the controller
     * will immediately stop.  The former condition can be diagnosed by
     * comparing current_dt() to min_dt(). The controller will also stop if it
     * encounters a step size for which <tt>isfinite</tt> is \c false.
     *
     * @param count_nt Number of time steps to take.  For unlimited
     *                 time step counts, use forever_nt().
     *
     * @return True if the controller completed the request successfully.
     *         False if the controller aborted for some reason.
     * @see The method advance() for a richer interface for controlling
     *      time advancement.
     */
    StopType step(StepType count_nt = 1)
    {
        return advance(std::numeric_limits<TimeType>::max(),
                       add_and_coerce_overflow_to_max(count_nt, current_nt()));
    }

    //@{

    /**
     * Retrieve the current simulation time.
     *
     * @return The current simulation time.
     */
    TimeType current_t() const { return current_t_; }

    /**
     * Retrieve the current simulation time step.  This is the size of the most
     * recent time step taken and is undefined prior to any time advancement
     * occurring via advance() or step().
     *
     * @return The current simulation time step.
     */
    TimeType current_dt() const { return current_dt_; }

    /**
     * Retrieve the current simulation time step.
     *
     * @return The current simulation discrete time step.
     */
    StepType current_nt() const {
        return boost::accumulators::extract::count(dt_stats);
    }

    //@}

    //@{

    /**
     * Retrieve the minimum time step taken during time advancement
     *
     * @return The minimum time step taken during time advancement
     */
    TimeType taken_min() const {
        return boost::accumulators::extract::min(dt_stats);
    }

    /**
     * Retrieve the mean time step taken during time advancement
     *
     * @return The mean time step taken during time advancement
     */
    TimeType taken_mean() const {
        return boost::accumulators::extract::mean(dt_stats);
    }

    /**
     * Retrieve the maximum time step taken during time advancement
     *
     * @return The maximum time step taken during time advancement
     */
    TimeType taken_max() const {
        return boost::accumulators::extract::max(dt_stats);
    }

    /**
     * Retrieve the standard deviation of the time steps taken during
     * time advancement.
     *
     * @return the standard deviation of the time steps taken during
     * time advancement.
     */
    TimeType taken_stddev() const {
        return std::sqrt(boost::accumulators::extract::variance(dt_stats));
    }

    //@}

private:

    // Mark Entry as noncopyable to avoid accidental performance hits
    struct Entry : public boost::noncopyable {
        bool periodic;
        TimeType every_dt, next_t;
        StepType every_nt, next_nt;
        boost::function<StopType (TimeType t, StepType nt)> callback;
    };

    // ptr_container to explicitly manage noncopyable instances with minimal
    // overhead and to automatically free contained resources on destruction.
    // ptr_vector to keep data as contiguous as possible in memory.
    typedef boost::ptr_vector<Entry> EntryList;

    typename boost::function<TimeType (TimeType)> stepper_;
    TimeType min_dt_, max_dt_, current_t_, current_dt_;
    EntryList entries_;

    // Maintain running statistics on the actual time step sizes taken
    // Also stores the current time step counter as a convenient side effect
    boost::accumulators::accumulator_set<
            TimeType,
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

template< typename TimeType, typename StepType, typename StopType >
template< typename StepperType >
timecontroller<TimeType,StepType,StopType>::timecontroller(
        StepperType stepper,
        TimeType initial_t,
        TimeType min_dt,
        TimeType max_dt)
    : stepper_(stepper),
      min_dt_(min_dt != 0 ? min_dt : std::numeric_limits<TimeType>::epsilon()),
      max_dt_(max_dt != 0 ? max_dt : std::numeric_limits<TimeType>::max()),
      current_t_(initial_t),
      entries_(7)
{
}

template< typename TimeType, typename StepType, typename StopType >
template< typename CallbackType >
void timecontroller<TimeType,StepType,StopType>::add_callback(
        TimeType what_t,
        StepType what_nt,
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

    // Avoid runtime costs for callbacks that should never occur.
    if (    what_t  == std::numeric_limits<TimeType>::max()
         && what_nt == std::numeric_limits<StepType>::max()) {
        return;
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

template< typename TimeType, typename StepType, typename StopType >
template< typename CallbackType >
void timecontroller<TimeType,StepType,StopType>::add_periodic_callback(
        TimeType every_dt,
        StepType every_nt,
        CallbackType callback)
{
    // Register the callback using first_dt = every_dt
    return add_periodic_callback(every_dt, every_dt, every_nt, callback);
}

template< typename TimeType, typename StepType, typename StopType >
template< typename CallbackType >
void timecontroller<TimeType,StepType,StopType>::add_periodic_callback(
        TimeType first_dt,
        TimeType every_dt,
        StepType every_nt,
        CallbackType callback)
{
    // first_dt after every_dt to improve messages when delegation target
    if (every_dt <= 0) throw std::invalid_argument("every_dt <= 0");
    if (first_dt <= 0) throw std::invalid_argument("first_dt <= 0");
    if (every_nt <= 0) throw std::invalid_argument("every_nt <= 0");

    // Avoid runtime costs for callbacks that should never occur.
    // Notice check occurs before add_and_coerce_overflow_to_max.
    if (    every_dt == std::numeric_limits<TimeType>::max()
         && every_nt == std::numeric_limits<StepType>::max()) {
        return;
    }

    Entry *e    = new Entry;      // Allocate Entry on heap
    e->periodic = true;
    e->every_dt = every_dt;
    e->every_nt = every_nt;
    e->next_t   = add_and_coerce_overflow_to_max(current_t_,   first_dt);
    e->next_nt  = add_and_coerce_overflow_to_max(current_nt(), every_nt);
    e->callback = callback;
    entries_.push_back(e);        // Transfer Entry memory ownership
}

template< typename TimeType, typename StepType, typename StopType >
StopType timecontroller<TimeType,StepType,StopType>::advance(
        const TimeType final_t,
        const StepType final_nt)
{
    SUZERAIN_TIMER_SCOPED("timecontroller::advance");

    SUZERAIN_ENSURE(min_dt_ <= max_dt_);
    using std::min;

    // Maintain the next simulation time something interesting must happen
    TimeType next_event_t = std::numeric_limits<TimeType>::max();

    // Find the simulation time of the first callback
    for (typename EntryList::iterator iter = entries_.begin();
         iter != entries_.end();
         ++iter) {
        next_event_t = min(next_event_t, (*iter).next_t);
    }

    // Advance time until done or we abort for some reason
    while (current_t_ < final_t && current_nt() < final_nt) {

        // Determine maximum possible step size allowed by all criteria
        const TimeType possible_dt
            = min(max_dt_, min(final_t, next_event_t) - current_t_);
        SUZERAIN_ENSURE(possible_dt > 0);

        // Take a single step, record new step and time, and update statistics
        current_dt_ = stepper_(possible_dt);
        current_t_ += current_dt_;
        dt_stats(current_dt_);

        // Sanity check the time step we just took
        if (SUZERAIN_UNLIKELY(!(boost::math::isfinite)(current_dt_))) {
            return false; // By design !isfinite(current_t_)
        }
        SUZERAIN_ENSURE(current_dt_ <= possible_dt);

        // Check callbacks and determine next callback simulation time
        next_event_t = std::numeric_limits<TimeType>::max();
        typename EntryList::iterator iter = entries_.begin();
        while (iter != entries_.end()) {

            // Callback required?
            if (SUZERAIN_UNLIKELY(    current_t_   == (*iter).next_t
                                   || current_nt() == (*iter).next_nt)) {

                // Perform required callback
                // Must perform state updates prior to any possible abort
                const StopType keep_advancing
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
        if (SUZERAIN_UNLIKELY(   current_dt_ < min_dt_
#pragma warning(push,disable:1572)
                              && current_dt_ != possible_dt)) {
#pragma warning(pop)
            return false; // On return current_dt() < min_dt()
        }
    }

    // Successfully advanced up to provided criteria
    return StopType(1);
}

} // namespace suzerain

#endif // SUZERAIN_TIMECONTROLLER_HPP
