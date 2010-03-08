/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * profile.cpp: miscellaneous utilities for profiling parallel transposes
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/validation.hpp>
#include <suzerain/profile.hpp>

namespace suzerain {

ProfileDefinition::ProfileDefinition()
    : options_(std::string("Profiling definition")),
      howmany_(2),                       // Default is one complex field
      backward_(0),
      forward_(0),
      roundtrip_(1),                     // One round trip
      nrep_(1)                           // One repetition
{
    namespace po = ::boost::program_options;

    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;

    options_.add_options()
        ("howmany", po::value<int>(&howmany_)
            ->notifier(bind2nd(ptr_fun(ensure_positive<int>),"howmany"))
            ->default_value(howmany_),
        "Number of real-valued components per state value")
        ("backward", po::value<int>(&backward_)
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"backward"))
            ->default_value(backward_),
         "Number of state fields to only convert from wave to physical space")
        ("forward", po::value<int>(&forward_)
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"forward"))
            ->default_value(forward_),
         "Number of state fields to only convert from physical to wave space")
        ("roundtrip", po::value<int>(&roundtrip_)
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<int>),"roundtrip"))
            ->default_value(roundtrip_),
         "Number of state fields to only convert from physical to wave space")
        ("nrep", po::value<int>(&nrep_)
            ->notifier(bind2nd(ptr_fun(ensure_positive<int>),"nrep"))
            ->default_value(nrep_),
         "Number of repetitions to perform for timing purposes");
    ;
}

} // namespace suzerain
