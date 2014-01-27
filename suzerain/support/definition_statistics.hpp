//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_SUPPORT_STATISTICS_DEFINITION_HPP
#define SUZERAIN_SUPPORT_STATISTICS_DEFINITION_HPP

/** @file
 * Classes handling control of statistics file output.
 */

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>

namespace suzerain {

namespace support {

/**
 * Encapsulates flags related to statistics output behavior for the simulation.
 * Includes the sampling rate and particulars around writing and archiving
 * statistics files.
 */
class definition_statistics
    : public virtual definition_base
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
     * @param nt           Number of simulation steps between samples.
     * @param dt           Amount of simulation time between samples.
     * @param final        Should a final sample be taken after
     *                     time advance successfully completes?
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      \c desttemplate and \c retain.
     */
    explicit definition_statistics(
            const std::string& destination            ,
            const std::size_t  retain      = (1 << 15),
            const real_t       dt          = 0        ,
            const std::size_t  nt          = 0        ,
            const bool         final       = false    );

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

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
    std::size_t retain;

    /**
     * The maximum amount of simulation time between statistical samples.
     */
    real_t dt;

    /**
     * The maximum number of simulation steps to take between statistical
     * samples.
     */
    std::size_t nt;

    /**
     * Should a final sample be taken after time advance successfully completes?
     */
    bool final;
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_STATISTICS_DEFINITION_HPP
