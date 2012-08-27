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
// profile.cpp: miscellaneous utilities for profiling parallel transposes
// $Id$

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
