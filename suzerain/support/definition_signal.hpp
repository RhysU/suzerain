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

#ifndef SUZERAIN_SUPPORT_SIGNAL_DEFINITION_HPP
#define SUZERAIN_SUPPORT_SIGNAL_DEFINITION_HPP

/** @file
 * Classes handling signal processing definitions.
 */

#include <suzerain/common.hpp>
#include <suzerain/support/definition_base.hpp>

namespace suzerain {

namespace support {

/**
 * Parses and contains signal processing details to allow a program to take
 * runtime-configurable responses to POSIX signals per \c signal.h.
 */
class signal_definition
    : public virtual definition_base
{
public:

    /**
     * Construct an instance with the given default behaviors.
     *
     * @param specstatus     Comma-separated signal names indicating a status
     *                       message should be displayed.
     * @param specrestart    Comma-separated signal names indicating a restart
     *                       file should be written.
     * @param specstatistics Comma-separated signal names indicating statistics
     *                       should be sampled and output to file.
     * @param specteardown   Comma-separated signal names indicating the
     *                       simulation should be torn down.
     * @param spechalt       Comma-separated signal names indicating the
     *                       simulation should be halted.
     */
    explicit signal_definition(
            const std::string& specstatus      = "HUP",
            const std::string& specrestart     = "HUP",
            const std::string& specstatistics  = "",
            const std::string& specteardown    = "USR1,USR2,TERM",
            const std::string& spechalt        = "INT");

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

    /** Signal numbers indicating a status message should be displayed. */
    std::vector<int> status;

    /** Signal numbers indicating a restart file should be written. */
    std::vector<int> restart;

    /** Signal numbers indicating a statistics file should be written. */
    std::vector<int> statistics;

    /** Signal numbers indicating the simulation should be torn down. */
    std::vector<int> teardown;

    /** Signal numbers indicating the simulation should be halted. */
    std::vector<int> halt;

private:

    /** Parse a status specification into this->status */
    void parse_status(const std::string& spec);

    /** Parse a status specification into this->restart */
    void parse_restart(const std::string& spec);

    /** Parse a status specification into this->statistics */
    void parse_statistics(const std::string& spec);

    /** Parse a status specification into this->teardown */
    void parse_teardown(const std::string& spec);

    /** Parse a status specification into this->halt */
    void parse_halt(const std::string& spec);
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_SIGNAL_DEFINITION_HPP
