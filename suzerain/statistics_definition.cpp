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
// statistics_definition.hpp: classes handling statistics output definitions
// $Id$

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
