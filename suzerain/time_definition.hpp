/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * time_definition.hpp: classes handling advance definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_TIME_DEFINITION_HPP
#define __SUZERAIN_TIME_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/validation.hpp>
#include <suzerain/exprparse.hpp>

/** @file
 * Provides classes handling time advancement settings.
 */

namespace suzerain {

namespace problem {

/**
 * Encapsulates flags related to time advancement.  Includes details on how
 * far the simulation should be advanced as long as how frequently status
 * updates should occur.
 */
template< typename FPT = double >
class TimeDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default values.  All of these can
     * be overridden by command line options.  Values not appearing here
     * <i>will</i> appear on the command line with integers and floats
     * defaulting to zero and NaN, respectively.
     *
     * @param advance_dt  Maximum amount of physical time to advance
     *                    the simulation.
     * @param advance_nt  Maximum number of discrete time steps to
     *                    advance the simulation.
     * @param advance_wt  Maximum amount of wall time to advance
     *                    the simulation.
     *                    advance the simulation.
     * @param status_dt   Maximum physical time between status updates.
     * @param status_nt   Maximum number of discrete time steps between
     *                    status updates.
     * @param min_dt      Minimum allowable physically-driven time step.
     * @param max_dt      Maximum allowable physically-driven time step.
     */
    TimeDefinition(FPT advance_dt,
                   int advance_nt,
                   FPT advance_wt,
                   FPT status_dt,
                   int status_nt,
                   FPT min_dt,
                   FPT max_dt);

    /**
     * Construct an instance with the given default values.  All of these can
     * be overridden by command line options.  Values not appearing here
     * <i>will not</i> appear on the command line with integers and floats
     * defaulting to zero and NaN, respectively.
     *
     * @param evmagfactor Safety factor used to pre-multiply a scheme's
     *                    maximum pure real and imaginary eigenvalue
     *                    magnitudes.  Usually in <tt>(0,1]</tt>, this
     *                    increases the conservativeness of the time stepping.
     */
    TimeDefinition(const char * evmagfactor);

    /** Maximum amount of physical time to advance the simulation. */
    FPT advance_dt;

    /** Maximum number of discrete time steps to advance the simulation. */
    int advance_nt;

    /** Maximum amount of wall time to advance the simulation. */
    FPT advance_wt;

    /** Maximum physical time between status updates. */
    FPT status_dt;

    /** Maximum number of discrete time steps between status updates. */
    int status_nt;

    /** Minimum allowable physically-driven time step */
    FPT min_dt;

    /** Maximum allowable physically-driven time step */
    FPT max_dt;

    /**
     * Factor used to pre-multiply a scheme's pure real and imaginary
     * eigenvalue magnitudes.
     *
     * @see suzerain::timestepper::lowstorage::SMR91Method for more details.
     */
    FPT evmagfactor;

private:

    /** Prepare one-time options for time advancement */
    void initialize_advancement(FPT default_advance_dt,
                                int default_advance_nt,
                                FPT default_advance_wt,
                                FPT default_status_dt,
                                int default_status_nt,
                                FPT default_min_dt,
                                FPT default_max_dt);

    /** Prepare repeatedly-used options similar to scenario parameters */
    void initialize_scenario(const char * default_evmagfactor);

    /** Helper used to parse string-based options */
    template<typename T>
    static void parse_option(const std::string &s,
                             T *value, void (*validator)(T, const char *),
                             const char *name)
    {
#pragma warning(push,disable:2259)
        const T t = suzerain::exprparse<FPT>(s, name);
#pragma warning(pop)
        validator(t, name);
        *value = t;
    }

    /** Helper used to parse <tt>dd:hh:mm:ss.ss</tt>-like options */
    static void parse_walltime(const std::string &s,
                               FPT *value,
                               const char *name)
    {
        // Parses dd:hh:mm:ss.ss, ..., mm:ss.ss, or ss.ss into seconds.
        // Additionally allows any of dd, hh, mm, or ss.ss to be expressions.

        using boost::algorithm::is_any_of;
        using boost::algorithm::split;
        using boost::algorithm::trim;
        using boost::bind;
        using std::back_inserter;
        using std::for_each;
        using std::invalid_argument;
        using std::locale;
        using std::remove;
        using std::string;
        using std::transform;
        using std::vector;

        // Split colon-separated spec into whitespace-trimmed, non-empty tokens
        vector<string> tokens;
        split(tokens, s, is_any_of(":"));
        for_each(tokens.begin(), tokens.end(), bind(trim<string>, _1, locale()));
        tokens.erase(remove(tokens.begin(), tokens.end(), ""), tokens.end());

        // Check incoming string lengths
        if (tokens.size() < 1) throw invalid_argument(
                string(name) + " did not specify a wall time");
        if (tokens.size() > 4) throw invalid_argument(
                string(name) + " not of format [dd:[hh:[mm:]]]ss.ss");

        // Parse dd, hh, mm, and ss into components[0], [1], [2], and [3]
        vector<FPT> components(4 - tokens.size(), 0);  // Zero missing values
        components.reserve(4);                         // Preallocate storage
        FPT (*f)(const std::string&, const char *) = &suzerain::exprparse<FPT>;
        transform(tokens.begin(), tokens.end(), back_inserter(components),
                  bind(f, _1, name));

        // Convert components to floating seconds
        FPT t =      components[0];  // days
        t = 24 * t + components[1];  // hours
        t = 60 * t + components[2];  // minutes
        t = 60 * t + components[3];  // seconds

        // Validate and store result
        suzerain::validation::ensure_nonnegative(t, name);
        *value = t;
    }

};

template< typename FPT >
void TimeDefinition<FPT>::initialize_advancement(
        FPT default_advance_dt,
        int default_advance_nt,
        FPT default_advance_wt,
        FPT default_status_dt,
        int default_status_nt,
        FPT default_min_dt,
        FPT default_max_dt)
{
    advance_dt = default_advance_dt;
    advance_nt = default_advance_nt;
    advance_wt = default_advance_wt;
    status_dt  = default_status_dt;
    status_nt  = default_status_nt;
    min_dt     = default_min_dt;
    max_dt     = default_max_dt;

    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::value;
    using std::string;
    using suzerain::validation::ensure_nonnegative;
    using suzerain::validation::ensure_positive;

    this->add_options()
        ("advance_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<FPT>, _1, &advance_dt,
                            &ensure_nonnegative<FPT>, "advance_dt"))
            ->default_value(lexical_cast<string>(advance_dt)),
         "Maximum amount of physical time to advance the simulation")
        ("advance_nt", value<string>(NULL)
            ->notifier(bind(&parse_option<int>, _1, &advance_nt,
                            &ensure_nonnegative<int>, "advance_nt"))
            ->default_value(lexical_cast<string>(advance_nt)),
         "Maximum number of discrete time steps to advance the simulation")
        ("advance_wt", value<string>(NULL)
            ->notifier(bind(&parse_walltime, _1, &advance_wt, "advance_wt"))
            ->default_value(lexical_cast<string>(advance_wt)),
            "Maximum amount of wall time to advance the simulation"
            " as [dd:[hh:[mm:]]]ss.ss")
        ("status_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<FPT>, _1, &status_dt,
                            &ensure_nonnegative<FPT>, "status_dt"))
            ->default_value(lexical_cast<string>(status_dt)),
         "Maximum physical time between status updates")
        ("status_nt", value<string>(NULL)
            ->notifier(bind(&parse_option<int>, _1, &status_nt,
                            &ensure_nonnegative<int>, "status_nt"))
            ->default_value(lexical_cast<string>(status_nt)),
         "Maximum number of discrete time steps between status updates")
        ("min_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<FPT>, _1, &min_dt,
                            &ensure_nonnegative<FPT>, "min_dt"))
            ->default_value(lexical_cast<string>(min_dt)),
         "Minimum allowable physically-driven time step")
        ("max_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<FPT>, _1, &max_dt,
                            &ensure_nonnegative<FPT>, "max_dt"))
            ->default_value(lexical_cast<string>(max_dt)),
         "Maximum allowable physically-driven time step")
    ;
}

template< typename FPT >
void TimeDefinition<FPT>::initialize_scenario(
        const char * default_evmagfactor)
{
    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    std::auto_ptr<boost::program_options::typed_value<std::string> > p;
    using suzerain::validation::ensure_positive;

    // evmagfactor
    evmagfactor = std::numeric_limits<FPT>::quiet_NaN();
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_option<FPT>, _1, &evmagfactor,
                            &ensure_positive<FPT>, "evmagfactor"));
    if (default_evmagfactor) p->default_value(default_evmagfactor);
    this->add_options()("evmagfactor", p.release(),
         "Safety factor in (0,1] used to adjust time step aggressiveness");
}

template< typename FPT >
TimeDefinition<FPT>::TimeDefinition(FPT advance_dt,
                                    int advance_nt,
                                    FPT advance_wt,
                                    FPT status_dt,
                                    int status_nt,
                                    FPT min_dt,
                                    FPT max_dt)
    : IDefinition("Time advancement parameters")
{
    initialize_advancement(advance_dt,
                           advance_nt,
                           advance_wt,
                           status_dt,
                           status_nt,
                           min_dt,
                           max_dt);

    initialize_scenario(NULL);
}

template< typename FPT >
TimeDefinition<FPT>::TimeDefinition(const char * evmagfactor)
    : IDefinition("Time advancement parameters")
{
    advance_dt = std::numeric_limits<FPT>::quiet_NaN();
    advance_nt = 0;
    advance_wt = std::numeric_limits<FPT>::quiet_NaN();
    status_dt  = std::numeric_limits<FPT>::quiet_NaN();
    status_nt  = 0;
    min_dt     = std::numeric_limits<FPT>::quiet_NaN();
    max_dt     = std::numeric_limits<FPT>::quiet_NaN();

    initialize_scenario(evmagfactor);
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_TIME_DEFINITION_HPP
