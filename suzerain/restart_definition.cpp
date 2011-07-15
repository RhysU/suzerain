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

#include <suzerain/common.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace problem {

RestartDefinition::RestartDefinition(
        const std::string& metadata,
        const std::string& uncommitted,
        const std::string& desttemplate,
        int retain,
        double restart_dt,
        int restart_nt)
    : IDefinition("Restart-related parameters"),
      metadata(metadata),
      uncommitted(uncommitted),
      desttemplate(desttemplate),
      retain(retain),
      restart_dt(restart_dt),
      restart_nt(restart_nt)
{
    namespace po = ::boost::program_options;

    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_positive;
    using ::suzerain::validation::ensure_nonnegative;

    ::std::pointer_to_binary_function<int,const char*,void>
        ptr_fun_ensure_positive_int(ensure_positive<int>);
    ::std::pointer_to_binary_function<int,const char*,void>
        ptr_fun_ensure_nonnegative_int(ensure_nonnegative<int>);
    ::std::pointer_to_binary_function<double,const char*,void>
        ptr_fun_ensure_nonnegative_double(ensure_nonnegative<double>);

    this->add_options()
        ("metadata", po::value<std::string>(&this->metadata)
            ->default_value(this->metadata),
         "Path to use when saving common restart metadata.  "
         "Any trailing \"XXXXXX\" will be used to generate a unique name."
         )
        ("uncommitted", po::value(&this->uncommitted)
            ->default_value(this->uncommitted),
         "Path to use when saving uncommitted restart data.  "
         "Any trailing \"XXXXXX\" will be used to generate a unique name.")
        ("desttemplate", po::value(&this->desttemplate)
            ->default_value(this->desttemplate),
         "Restart archiving pattern to use when committing restart files.  "
         "One or more #'s must be present and will be replaced by a sequence number.  "
         "Any trailing \"XXXXXX\" will be used to generate a unique template.")
        ("retain", po::value(&this->retain)
            ->notifier(bind2nd(ptr_fun_ensure_positive_int, "retain"))
            ->default_value(this->retain),
         "Maximum number of committed restart files to retain")
        ("restart_dt", po::value(&this->restart_dt)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_double, "restart_dt"))
            ->default_value(this->restart_dt),
         "Maximum amount of simulation time between restart files")
        ("restart_nt", po::value(&this->restart_nt)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_int, "restart_nt"))
            ->default_value(this->restart_nt),
         "Maximum number of time steps between restart files")
    ;
}

} // namespace problem

} // namespace suzerain
