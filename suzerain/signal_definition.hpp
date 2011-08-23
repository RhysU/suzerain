/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011 The PECOS Development Team
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
 * signal_definition.hpp: classes handling signal processing definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_SIGNAL_DEFINITION_HPP
#define __SUZERAIN_SIGNAL_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides classes handling signal processing definitions.
 */

namespace suzerain {

namespace problem {

/**
 * Parses and contains signal processing details to allow a program to take
 * runtime-configurable responses to POSIX signals per \c signal.h.
 */
class SignalDefinition : public IDefinition
{
public:

    /**
     * Construct an instance with the given default behaviors.
     *
     * @param specstatus   Comma-separated signal names indicating a status
     *                     message should be displayed.
     * @param specrestart  Comma-separated signal names indicating a restart
     *                     file should be written.
     * @param specteardown Comma-separated signal names indicating the
     *                     simulation should be exited.
     */
    explicit SignalDefinition(
            const std::string& specstatus   = "HUP",
            const std::string& specrestart  = "HUP",
            const std::string& specteardown = "INT,USR1,USR2,TERM");

    /** Signal numbers indicating a status message should be displayed. */
    std::vector<int> status;

    /** Signal numbers indicating a restart file should be written. */
    std::vector<int> restart;

    /** Signal numbers indicating the simulation should be exited. */
    std::vector<int> teardown;

private:

    /** Parse a status specification into this->status */
    void parse_status(const std::string &spec);

    /** Parse a status specification into this->restart */
    void parse_restart(const std::string &spec);

    /** Parse a status specification into this->teardown */
    void parse_teardown(const std::string &spec);
};

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_SIGNAL_DEFINITION_HPP
