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
// restart_definition.cpp: implementations handling restart definitions
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/validation.hpp>

namespace suzerain
{

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

restart_definition::restart_definition(
    const std::string& metadata,
    const std::string& uncommitted,
    const std::string& destination,
    std::size_t retain,
    real_t dt,
    std::size_t nt)
    : definition_base("Restart-related parameters"),
      metadata(metadata),
      uncommitted(uncommitted),
      destination(destination),
      retain(retain),
      dt(dt),
      nt(nt),
      physical(false)
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::value;
    using boost::program_options::bool_switch;
    using std::string;
    using validation::ensure_nonnegative;

    this->add_options()
    ("metadata", value(&this->metadata)
     ->default_value(this->metadata),
     "Path to use when saving common metadata for output files.  "
     "Any trailing \"XXXXXX\" will be used to generate a unique name.")
    ("uncommitted", value(&this->uncommitted)
     ->default_value(this->uncommitted),
     "Path to use when saving uncommitted output files.  "
     "Any trailing \"XXXXXX\" will be used to generate a unique name.")
    ("restart_destination", value(&this->destination)
     ->default_value(this->destination),
     "Archiving destination to use when committing restart files.  "
     "One or more #'s must be present and will be replaced by a sequence number.  "
     "Any trailing \"XXXXXX\" will be used to generate a unique template.")
    ("restart_retain", value<string>(NULL)
     ->notifier(bind(&parse_size_t, _1, &this->retain, "restart_retain"))
     ->default_value(lexical_cast<string>(this->retain)),
     "Maximum number of committed restart files to retain")
    ("restart_dt", value<string>(NULL)
     ->notifier(bind(&parse_option<real_t>, _1, &this->dt,
                     &ensure_nonnegative<real_t>, "restart_dt"))
     ->default_value(lexical_cast<string>(this->dt)),
     "Maximum amount of simulation time between restart files")
    ("restart_nt", value<string>(NULL)
     ->notifier(bind(&parse_size_t, _1, &this->nt, "restart_nt"))
     ->default_value(lexical_cast<string>(this->nt)),
     "Maximum number of time steps between restart files")
    ("restart_physical", bool_switch(&this->physical),
     "Specify flag to save restart fields as primitive variables "
     "stored at collocation points in physical space")
    ;
}

} // namespace suzerain
