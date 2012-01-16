/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 The PECOS Development Team
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
#ifndef __SUZERAIN_STATISTICS_DEFINITION_HPP
#define __SUZERAIN_STATISTICS_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides classes handling arguments used to control statistics file output.
 */

namespace suzerain {

namespace problem {

/**
 * Encapsulates flags related to statistics output behavior for the simulation.
 * Includes the sampling rate and particulars around writing and archiving
 * statistics files.
 */
class StatisticsDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default values.
     * All of these can be overridden by command line options.
     *
     * @param destination  Archiving pattern to use when
     *                     committing statistics output files.
     * @param retain       Maximum number of committed statistics files
     *                     to retain.
     * @param nt           Number of simulation steps between sampling.
     * @param dt           Amount of simulation time between sampling.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      \c desttemplate and \c retain.
     */
    explicit StatisticsDefinition(
            const std::string& destination,
            int retain = (1<<15),
            double dt  = 1,
            int nt     = 0);

    /**
     * The archiving pattern to use when committing statistics output files.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      destination.
     */
    std::string destination;

    /**
     * The maximum number of committed statistics files to retain.
     */
    int retain;

    /**
     * The maximum amount of simulation time between sampling statistics.
     */
    double dt;

    /**
     * The maximum number of simulation steps to take between sampling
     * statistics.
     */
    int nt;
};

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_STATISTICS_DEFINITION_HPP
