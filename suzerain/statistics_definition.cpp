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
 * statistics_definition.hpp: classes handling statistics output definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/statistics_definition.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace problem {

template<typename T>
static void parse_option(const std::string &s,
                         T *value, void (*validator)(T, const char *),
                         const char *name)
{
#pragma warning(push,disable:2259)
    const T t = suzerain::exprparse<double>(s, name);
#pragma warning(pop)
    validator(t, name);
    *value = t;
}

StatisticsDefinition::StatisticsDefinition(
        const std::string& destination,
        int retain,
        double dt,
        int nt)
    : IDefinition("Statistics sampling parameters"),
      destination(destination),
      retain(retain),
      dt(dt),
      nt(nt)
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::value;
    using boost::program_options::bool_switch;
    using std::string;
    using suzerain::validation::ensure_nonnegative;
    using suzerain::validation::ensure_positive;

    this->add_options()
        ("statistics_destination", value(&this->destination)
            ->default_value(this->destination),
         "Archiving destination to use when committing statistics files.  "
         "One or more #'s must be present and will be replaced by a sequence number.  "
         "Any trailing \"XXXXXX\" will be used to generate a unique template.")
        ("statistics_retain", value<string>(NULL)
            ->notifier(bind(&parse_option<int>, _1, &this->retain,
                            &ensure_nonnegative<int>, "statistics_retain"))
            ->default_value(lexical_cast<string>(this->retain)),
         "Maximum number of committed statistics files to retain")
        ("statistics_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<double>, _1, &this->dt,
                            &ensure_nonnegative<double>, "statistics_dt"))
            ->default_value(lexical_cast<string>(this->dt)),
         "Maximum amount of simulation time between sampling statistics")
        ("statistics_nt", value<string>(NULL)
            ->notifier(bind(&parse_option<int>, _1, &this->nt,
                            &ensure_nonnegative<int>, "statistics_nt"))
            ->default_value(lexical_cast<string>(this->nt)),
         "Maximum number of time steps between sampling statistics")
    ;
}

} // namespace problem

} // namespace suzerain
