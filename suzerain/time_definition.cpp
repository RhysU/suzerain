//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// time_definition.cpp: classes handling advance definitions
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/validation.hpp>
#include <suzerain/exprparse.hpp>

/** @file
 * Provides classes handling time advancement settings.
 */

namespace suzerain {

/** Helper used to parse string-based options */
template<typename T>
static void parse_option(const std::string &s,
                         T *value, void (*validator)(T, const char *),
                         const char *name)
{
#pragma warning(push,disable:2259)
    const T t = exprparse<real_t>(s, name);
#pragma warning(pop)
    validator(t, name);
    *value = t;
}

/** Helper used to parse <tt>dd:hh:mm:ss.ss</tt>-like options */
static void parse_walltime(const std::string &s,
                           real_t *value,
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
    vector<real_t> components(4 - tokens.size(), 0);  // Zero missing values
    components.reserve(4);                            // Preallocate storage
    real_t (*f)(const std::string&, const char *) = &exprparse<real_t>;
    transform(tokens.begin(), tokens.end(), back_inserter(components),
              bind(f, _1, name));

    // Convert components to floating seconds
    real_t t =   components[0];  // days
    t = 24 * t + components[1];  // hours
    t = 60 * t + components[2];  // minutes
    t = 60 * t + components[3];  // seconds

    // Validate and store result
    validation::ensure_nonnegative(t, name);
    *value = t;
}

namespace problem {

void time_definition::initialize_advancement(
        real_t default_advance_dt,
        int    default_advance_nt,
        real_t default_advance_wt,
        real_t default_status_dt,
        int    default_status_nt,
        real_t default_min_dt,
        real_t default_max_dt)
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
    using validation::ensure_nonnegative;
    using validation::ensure_positive;

    this->add_options()
        ("advance_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<real_t>, _1, &advance_dt,
                            &ensure_nonnegative<real_t>, "advance_dt"))
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
            ->notifier(bind(&parse_option<real_t>, _1, &status_dt,
                            &ensure_nonnegative<real_t>, "status_dt"))
            ->default_value(lexical_cast<string>(status_dt)),
         "Maximum physical time between status updates")
        ("status_nt", value<string>(NULL)
            ->notifier(bind(&parse_option<int>, _1, &status_nt,
                            &ensure_nonnegative<int>, "status_nt"))
            ->default_value(lexical_cast<string>(status_nt)),
         "Maximum number of discrete time steps between status updates")
        ("min_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<real_t>, _1, &min_dt,
                            &ensure_nonnegative<real_t>, "min_dt"))
            ->default_value(lexical_cast<string>(min_dt)),
         "Minimum allowable physically-driven time step")
        ("max_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<real_t>, _1, &max_dt,
                            &ensure_nonnegative<real_t>, "max_dt"))
            ->default_value(lexical_cast<string>(max_dt)),
         "Maximum allowable physically-driven time step")
    ;
}

void time_definition::initialize_scenario(
        const char * default_evmagfactor)
{
    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    std::auto_ptr<boost::program_options::typed_value<std::string> > p;
    using validation::ensure_positive;

    // evmagfactor
    evmagfactor = std::numeric_limits<real_t>::quiet_NaN();
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_option<real_t>, _1, &evmagfactor,
                            &ensure_positive<real_t>, "evmagfactor"));
    if (default_evmagfactor) p->default_value(default_evmagfactor);
    this->add_options()("evmagfactor", p.release(),
         "Safety factor in (0,1] used to adjust time step aggressiveness");
}

time_definition::time_definition(real_t advance_dt,
                               int    advance_nt,
                               real_t advance_wt,
                               real_t status_dt,
                               int    status_nt,
                               real_t min_dt,
                               real_t max_dt)
    : definition_base("Time advancement parameters")
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

time_definition::time_definition(const char * evmagfactor)
    : definition_base("Time advancement parameters")
{
    advance_dt = std::numeric_limits<real_t>::quiet_NaN();
    advance_nt = 0;
    advance_wt = std::numeric_limits<real_t>::quiet_NaN();
    status_dt  = std::numeric_limits<real_t>::quiet_NaN();
    status_nt  = 0;
    min_dt     = std::numeric_limits<real_t>::quiet_NaN();
    max_dt     = std::numeric_limits<real_t>::quiet_NaN();

    initialize_scenario(evmagfactor);
}

} // namespace problem

} // namespace suzerain
