//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

/** @file
 * @copydoc definition_time.hpp
 */

#include <suzerain/support/definition_time.hpp>

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

/** Helper used to parse size_t-based options */
static void parse_size_t(const std::string& s,
                         std::size_t* value,
                         const char* name)
{
    using std::floor;
#pragma warning(push,disable:2259)
    const real_t t = floor(exprparse<real_t>(s, name) + real_t(1)/2);
#pragma warning(pop)
    validation::ensure_nonnegative(t, name);
    *value = t;
}

/** Helper used to parse string-based options */
template<typename T>
static void parse_option(const std::string& s,
                         T* value,
                         void (*validator)(T, const char*),
                         const char* name)
{
#pragma warning(push,disable:2259)
    const T t = exprparse<real_t>(s, name);
#pragma warning(pop)
    validator(t, name);
    *value = t;
}

/** Helper used to parse <tt>dd:hh:mm:ss.ss</tt>-like options */
static void parse_walltime(const std::string& s,
                           real_t* value,
                           const char* name)
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
    real_t (*f)(const std::string&, const char*) = &exprparse<real_t>;
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

definition_time::definition_time(const real_t      advance_dt,
                                 const std::size_t advance_nt,
                                 const real_t      advance_wt,
                                 const real_t      status_dt,
                                 const std::size_t status_nt,
                                 const bool        status_final,
                                 const real_t      min_dt,
                                 const real_t      max_dt)
    : advance_dt  (advance_dt)
    , advance_nt  (advance_nt)
    , advance_wt  (advance_wt)
    , status_dt   (status_dt)
    , status_nt   (status_nt)
    , status_final(status_final)
    , min_dt      (min_dt)
    , max_dt      (max_dt)
    , evmagfactor (std::numeric_limits<real_t>::quiet_NaN())
    , constructor (constructor1)
{
}

definition_time::definition_time(const real_t evmagfactor)
    : advance_dt  (std::numeric_limits<real_t>::quiet_NaN())
    , advance_nt  (0)
    , advance_wt  (std::numeric_limits<real_t>::quiet_NaN())
    , status_dt   (std::numeric_limits<real_t>::quiet_NaN())
    , status_nt   (0)
    , status_final(true)
    , min_dt      (std::numeric_limits<real_t>::quiet_NaN())
    , max_dt      (std::numeric_limits<real_t>::quiet_NaN())
    , evmagfactor (evmagfactor)
    , constructor (constructor2)
{
}

// Descriptions used in options_description and possibly save/load.
static const char description_advance_dt[]
        = "Maximum amount of physical time to advance the simulation";

static const char description_advance_nt[]
        = "Maximum number of discrete time steps to advance the simulation";

static const char description_advance_wt[]
        = "Maximum amount of wall time to advance the simulation as [dd:[hh:[mm:]]]ss.ss";

static const char description_status_dt[]
        = "Maximum physical time between status updates";

static const char description_status_nt[]
        = "Maximum number of discrete time steps between status updates";

static const char description_status_final[]
        = "Should a final status update occur after advance completes?";

static const char description_min_dt[]
        = "Minimum allowable physically-driven time step";

static const char description_max_dt[]
        = "Maximum allowable physically-driven time step";

static const char description_evmagfactor[]
        = "Safety factor in (0,1] used to adjust time step aggressiveness";

boost::program_options::options_description
definition_time::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::numeric_limits;
    using std::string;
    using validation::ensure_nonnegative;
    using validation::ensure_positive;

    boost::program_options::options_description retval(
            "Time advancement parameters");

    // See constructor definition for rational behind conditional
    if (this->constructor != constructor2) {
        retval.add_options()
        ("advance_dt", value<string>(NULL)
         ->notifier(bind(&parse_option<real_t>, _1, &advance_dt,
                         &ensure_nonnegative<real_t>, "advance_dt"))
         ->default_value(lexical_cast<string>(advance_dt)),
         description_advance_dt)
        ("advance_nt", value<string>(NULL)
         ->notifier(bind(&parse_size_t, _1, &advance_nt, "advance_nt"))
         ->default_value(lexical_cast<string>(advance_nt)),
         description_advance_nt)
        ("advance_wt", value<string>(NULL)
         ->notifier(bind(&parse_walltime, _1, &advance_wt, "advance_wt"))
         ->default_value(lexical_cast<string>(advance_wt)),
         description_advance_wt)
        ("status_dt", value<string>(NULL)
         ->notifier(bind(&parse_option<real_t>, _1, &status_dt,
                         &ensure_nonnegative<real_t>, "status_dt"))
         ->default_value(lexical_cast<string>(status_dt)),
         description_status_dt)
        ("status_nt", value<string>(NULL)
         ->notifier(bind(&parse_size_t, _1, &status_nt, "status_nt"))
         ->default_value(lexical_cast<string>(status_nt)),
         description_status_nt)
        ("status_final", value<bool>(&status_final)
         ->default_value(status_final),
         description_status_final)
        ("min_dt", value<string>(NULL)
         ->notifier(bind(&parse_option<real_t>, _1, &min_dt,
                         &ensure_nonnegative<real_t>, "min_dt"))
         ->default_value(lexical_cast<string>(min_dt)),
         description_min_dt)
        ("max_dt", value<string>(NULL)
         ->notifier(bind(&parse_option<real_t>, _1, &max_dt,
                         &ensure_nonnegative<real_t>, "max_dt"))
         ->default_value(lexical_cast<string>(max_dt)),
         description_max_dt)
        ;
    }

    // Complicated add_options() calls allow changing the default value of
    // evmagfactor displayed depending on the constructor.  NaN used as a NOP
    // value by clients.  Validation routines used below silently allow NaNs.
    auto_ptr<typed_value<string> > p(value<string>());
    p->notifier(bind(&parse_option<real_t>, _1, &evmagfactor,
                     &ensure_positive<real_t>, "evmagfactor"));
    if (constructor == constructor2)
    {
        p->default_value(lexical_cast<string>(evmagfactor));
    }
    retval.add_options()
        ("evmagfactor", p.release(), description_evmagfactor)
    ;

    return retval;
}

void definition_time::save(
        const esio_handle h) const
{
    DEBUG0("Saving some, but not all, definition_time parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "evmagfactor", &evmagfactor, 0, description_evmagfactor);
}

void definition_time::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading some, but not all, definition_time parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    // evmagfactor not present in pre-revision 24043 restart files
    if (!(boost::math::isnan)(evmagfactor)) {
        if (verbose)
            INFO0("Overriding timedef using evmagfactor = " << evmagfactor);
    } else if (ESIO_NOTFOUND == esio_line_size(h, "evmagfactor", NULL)) {
        evmagfactor = 0.72;
        if (verbose)
            INFO0("Employing default evmagfactor = " << evmagfactor);
    } else {
        esio_line_read(h, "evmagfactor", &evmagfactor, 0);
    }
}

} // namespace support

} // namespace suzerain
