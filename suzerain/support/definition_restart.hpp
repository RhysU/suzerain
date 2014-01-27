//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

#ifndef SUZERAIN_SUPPORT_RESTART_DEFINITION_HPP
#define SUZERAIN_SUPPORT_RESTART_DEFINITION_HPP

/** @file
 * Provides classes handling restart definitions, which are runtime
 * arguments used to control simulation restart behavior.
 */

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>

namespace suzerain {

namespace support {

/**
 * Encapsulates flags related to restart file behavior for the simulation.
 * Includes the particulars around writing and archiving restart files.
 */
class definition_restart
    : public virtual definition_base
{
public:
    /**
     * Construct an instance with the given default values.
     * All of these can be overridden by command line options.
     *
     * @param metadata     Restart file path to use when saving
     *                     common restart file metadata.
     * @param uncommitted  Restart file path to use when saving
     *                     uncommitted restart data.
     * @param destination  Restart archiving pattern to use when
     *                     committing restart files.
     * @param retain       Maximum number of committed restart files
     *                     to retain.
     * @param nt           Number of simulation steps between restart writes.
     * @param dt           Amount of simulation time between restart writes.
     * @param final        Should a final restart write occur after
     *                     time advance successfully completes?
     * @param physical     Should files save state in physical space?
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      \c destination and \c retain.
     */
    definition_restart(const std::string& metadata                ,
                       const std::string& uncommitted             ,
                       const std::string& destination             ,
                       const std::size_t  retain       = (1 << 15),
                       const real_t       dt           = 0        ,
                       const std::size_t  nt           = 0        ,
                       const bool         final        = true     ,
                       const bool         physical     = false    );

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

    /**
     * The restart file path to use when saving common restart file
     * metadata.
     */
    std::string metadata;

    /**
     * The file path to use when saving uncommitted restart data.
     */
    std::string uncommitted;

    /**
     * The archiving pattern to use when committing restart files.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      destination.
     */
    std::string destination;

    /**
     * The maximum number of committed restart files to retain.
     */
    std::size_t retain;

    /**
     * The maximum amount of simulation time between writing restart files.
     */
    real_t dt;

    /**
     * The maximum number of simulation steps to take between writing restart
     * files.
     */
    std::size_t nt;

    /**
     * Should a final restart write occur after time advance successfully
     * completes?
     */
    bool final;

    /**
     * Save restart fields as collocation point values in physical space?
     */
    bool physical;
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_RESTART_DEFINITION_HPP
