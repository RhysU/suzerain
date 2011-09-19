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
 * restart_definition.cpp: implementations handling restart definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/restart_definition.hpp>
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

RestartDefinition::RestartDefinition(
        const std::string& metadata,
        const std::string& uncommitted,
        const std::string& desttemplate,
        int retain,
        double dt,
        int nt)
    : IDefinition("Restart-related parameters"),
      metadata(metadata),
      uncommitted(uncommitted),
      desttemplate(desttemplate),
      retain(retain),
      dt(dt),
      nt(nt)
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::value;
    using std::string;
    using suzerain::validation::ensure_nonnegative;
    using suzerain::validation::ensure_positive;

    this->add_options()
        ("metadata", value(&this->metadata)
            ->default_value(this->metadata),
         "Path to use when saving common restart metadata.  "
         "Any trailing \"XXXXXX\" will be used to generate a unique name."
         )
        ("uncommitted", value(&this->uncommitted)
            ->default_value(this->uncommitted),
         "Path to use when saving uncommitted restart data.  "
         "Any trailing \"XXXXXX\" will be used to generate a unique name.")
        ("desttemplate", value(&this->desttemplate)
            ->default_value(this->desttemplate),
         "Restart archiving pattern to use when committing restart files.  "
         "One or more #'s must be present and will be replaced by a sequence number.  "
         "Any trailing \"XXXXXX\" will be used to generate a unique template.")
        ("retain", value<string>(NULL)
            ->notifier(bind(&parse_option<int>, _1, &this->retain,
                            &ensure_nonnegative<int>, "retain"))
            ->default_value(lexical_cast<string>(this->retain)),
         "Maximum number of committed restart files to retain")
        ("restart_dt", value<string>(NULL)
            ->notifier(bind(&parse_option<double>, _1, &this->dt,
                            &ensure_nonnegative<double>, "restart_dt"))
            ->default_value(lexical_cast<string>(this->dt)),
         "Maximum amount of simulation time between restart files")
        ("restart_nt", value<string>(NULL)
            ->notifier(bind(&parse_option<int>, _1, &this->nt,
                            &ensure_nonnegative<int>, "restart_nt"))
            ->default_value(lexical_cast<string>(this->nt)),
         "Maximum number of time steps between restart files")
    ;
}

} // namespace problem

} // namespace suzerain
