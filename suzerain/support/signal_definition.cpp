//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
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

/** @file
 * @copydoc signal_definition.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/signal_definition.hpp>

#include <suzerain/error.h>
#include <suzerain/os.h>

namespace suzerain {

namespace support {

// Helper method to turn, e.g., "USR1, sigusr2, 3 , Term "
// into signal numbers, e.g. {SIGUSR1, SIGUSR2, 3, SIGTERM}.
static std::vector<int> parse_spec(const std::string& name,
                                   const std::string& spec)
{
    using boost::algorithm::is_any_of;
    using boost::algorithm::split;
    using boost::algorithm::trim;
    using boost::bind;
    using std::for_each;
    using std::locale;
    using std::remove;
    using std::string;
    using std::vector;

    // Split comma-separated spec into whitespace-trimmed, non-empty tokens
    vector<string> tokens;
    split(tokens, spec, is_any_of(","));
    for_each(tokens.begin(), tokens.end(), bind(trim<string>, _1, locale()));
    tokens.erase(remove(tokens.begin(), tokens.end(), ""), tokens.end());

    vector<int> retval;

    // Convert symbolic names or integer strings into signal numbers
    for (vector<string>::iterator iter = tokens.begin();
         iter != tokens.end(); ++iter)
    {

        const string& token = *iter;

        // First, try to lookup a named signal number
        int signum = suzerain_signal_number(token.c_str());
        if (signum == SUZERAIN_FAILURE)
        {
            // That failing, next, try to parse token as a positive integer
            errno = 0;
            char* end;
            const long l = strtol(token.c_str(), &end, 0);
            if (errno != 0 || *end != '\0' || l <= 0 || l > INT_MAX)
            {
                std::ostringstream oss;
                oss << "Signal specification '"
                    << token
                    << "' provided for "
                    << name
                    << " is neither a known signal"
                       " nor a strictly positive integer";
                throw std::invalid_argument(oss.str());
            }
            signum = l;
        }
        retval.push_back(signum);
    }

    return retval;
}

signal_definition::signal_definition(const std::string& specstatus,
                                     const std::string& specrestart,
                                     const std::string& specstatistics,
                                     const std::string& specteardown,
                                     const std::string& spechalt)
    : status    (parse_spec("signal_definition(specstatus,...)",
                            specstatus))
    , restart   (parse_spec("signal_definition(...,specrestart,...)",
                            specrestart))
    , statistics(parse_spec("signal_definition(...,specstatistics,...)",
                            specstatistics))
    , teardown  (parse_spec("signal_definition(...,specteardown)",
                            specteardown))
    , halt      (parse_spec("signal_definition(...,spechalt)",
                            spechalt))
{
}

boost::program_options::options_description
signal_definition::options_description()
{
    using boost::algorithm::join;
    using boost::program_options::options_description;
    using boost::program_options::value;
    using std::back_insert_iterator;
    using std::string;
    using std::transform;
    using std::vector;

    options_description retval("Actions to take on receipt of various signals");

    // Build human-readable default strings based on the current values
    vector<string> names_status;
    transform(status.begin(), status.end(),
              back_insert_iterator<vector<string> >(names_status),
              &suzerain_signal_name);

    vector<string> names_restart;
    transform(restart.begin(), restart.end(),
              back_insert_iterator<vector<string> >(names_restart),
              &suzerain_signal_name);

    vector<string> names_statistics;
    transform(statistics.begin(), statistics.end(),
              back_insert_iterator<vector<string> >(names_statistics),
              &suzerain_signal_name);

    vector<string> names_teardown;
    transform(teardown.begin(), teardown.end(),
              back_insert_iterator<vector<string> >(names_teardown),
              &suzerain_signal_name);

    vector<string> names_halt;
    transform(halt.begin(), halt.end(),
              back_insert_iterator<vector<string> >(names_halt),
              &suzerain_signal_name);

    // Construct the collection of options to be returned
    retval.add_options()
    ("signal_status", value<std::string>()
     ->default_value(join(names_status, ","))
     ->notifier(boost::bind(&signal_definition::parse_status, this, _1)),
     "Show status information on any signal in this comma-separated list")

    ("signal_restart", value<std::string>()
     ->default_value(join(names_restart, ","))
     ->notifier(boost::bind(&signal_definition::parse_restart, this, _1)),
     "Write restart file on any signal in this comma-separated list")

    ("signal_statistics", value<std::string>()
     ->default_value(join(names_statistics, ","))
     ->notifier(boost::bind(&signal_definition::parse_statistics, this, _1)),
     "Write statistics file on any signal in this comma-separated list")

    ("signal_teardown", value<std::string>()
     ->default_value(join(names_teardown, ","))
     ->notifier(boost::bind(&signal_definition::parse_teardown, this, _1)),
     "Tear down simulation on any signal in this comma-separated list")

    ("signal_halt", value<std::string>()
     ->default_value(join(names_halt, ","))
     ->notifier(boost::bind(&signal_definition::parse_halt, this, _1)),
     "Halt simulation on any signal in this comma-separated list")
    ;

    return retval;
}

void signal_definition::parse_status(const std::string& spec)
{
    std::vector<int> tmp = parse_spec("--signal_status", spec);
    this->status.swap(tmp);
}

void signal_definition::parse_restart(const std::string& spec)
{
    std::vector<int> tmp = parse_spec("--signal_restart", spec);
    this->restart.swap(tmp);
}

void signal_definition::parse_statistics(const std::string& spec)
{
    std::vector<int> tmp = parse_spec("--signal_statistics", spec);
    this->statistics.swap(tmp);
}

void signal_definition::parse_teardown(const std::string& spec)
{
    std::vector<int> tmp = parse_spec("--signal_teardown", spec);
    this->teardown.swap(tmp);
}

void signal_definition::parse_halt(const std::string& spec)
{
    std::vector<int> tmp = parse_spec("--signal_halt", spec);
    this->halt.swap(tmp);
}

} // end namespace support

} // end namespace suzerain
