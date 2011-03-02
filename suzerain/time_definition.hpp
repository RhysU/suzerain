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
 * advance_definition.hpp: classes handling advance definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_TIME_DEFINITION_HPP
#define __SUZERAIN_TIME_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling time advancement settings.
 */

namespace suzerain {

namespace problem {

/**
 * Encapsulates flags related to time advancement.  Includes details on how
 * far the simulation should be advanced as long as how frequently status
 * updates should occur.
 */
template< typename FPT = double >
class TimeDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default values.
     * All of these can be overridden by command line options.
     *
     * @param default_advance_dt Maximum amount of physical time to advance
     *                           the simulation.
     * @param default_advance_nt Maximum number of discrete time steps to
     *                           advance the simulation.
     * @param default_status_dt  Maximum physical time between status updates.
     * @param default_status_dt  Maximum number of discrete time steps between
     *                   status updates.
     */
    explicit TimeDefinition(FPT default_advance_dt = 0,
                            int default_advance_nt = 0,
                            FPT default_status_dt  = 0,
                            int default_status_nt  = 0);

    /** Maximum amount of physical time to advance the simulation. */
    FPT advance_dt;

    /** Maximum number of discrete time steps to advance the simulation. */
    int advance_nt;

    /** Maximum physical time between status updates. */
    FPT status_dt;

    /** Maximum number of discrete time steps between status updates. */
    int status_nt;
};

template< typename FPT >
TimeDefinition<FPT>::TimeDefinition(FPT default_advance_dt,
                                    int default_advance_nt,
                                    FPT default_status_dt,
                                    int default_status_nt)
    : IDefinition("Time advancement parameters"),
      advance_dt(default_advance_dt),
      advance_nt(default_advance_nt),
      status_dt(default_status_dt),
      status_nt(default_status_nt)
{
    namespace po = ::boost::program_options;

    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;

    ::std::pointer_to_binary_function<int,const char*,void>
        ptr_fun_ensure_nonnegative_int(ensure_nonnegative<int>);
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ptr_fun_ensure_nonnegative_FPT(ensure_nonnegative<FPT>);

    this->add_options()
        ("advance_dt", po::value(&advance_dt)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "advance_dt"))
            ->default_value(advance_dt),
         "Maximum amount of physical time to advance the simulation")
        ("advance_nt", po::value(&advance_nt)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_int, "advance_nt"))
            ->default_value(advance_nt),
         "Maximum number of discrete time steps to advance the simulation")
        ("status_dt", po::value(&status_dt)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT, "status_dt"))
            ->default_value(status_dt),
         "Maximum physical time between status updates")
        ("status_nt", po::value(&status_nt)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_int, "status_nt"))
            ->default_value(status_nt),
         "Maximum number of discrete time steps between status updates")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_TIME_DEFINITION_HPP
