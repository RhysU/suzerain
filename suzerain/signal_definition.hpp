//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// signal_definition.hpp: classes handling signal processing definitions
// $Id$

#ifndef SUZERAIN_SIGNAL_DEFINITION_HPP
#define SUZERAIN_SIGNAL_DEFINITION_HPP

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
     * @param specstatus     Comma-separated signal names indicating a status
     *                       message should be displayed.
     * @param specrestart    Comma-separated signal names indicating a restart
     *                       file should be written.
     * @param specstatistics Comma-separated signal names indicating statistics
     *                       should be sampled and output to file.
     * @param specteardown   Comma-separated signal names indicating the
     *                       simulation should be exited.
     */
    explicit SignalDefinition(
            const std::string& specstatus      = "HUP",
            const std::string& specrestart     = "HUP",
            const std::string& specstatistics  = "",
            const std::string& specteardown    = "INT,USR1,USR2,TERM");

    /** Signal numbers indicating a status message should be displayed. */
    std::vector<int> status;

    /** Signal numbers indicating a restart file should be written. */
    std::vector<int> restart;

    /** Signal numbers indicating a statistics file should be written. */
    std::vector<int> statistics;

    /** Signal numbers indicating the simulation should be exited. */
    std::vector<int> teardown;

private:

    /** Parse a status specification into this->status */
    void parse_status(const std::string &spec);

    /** Parse a status specification into this->restart */
    void parse_restart(const std::string &spec);

    /** Parse a status specification into this->statistics */
    void parse_statistics(const std::string &spec);

    /** Parse a status specification into this->teardown */
    void parse_teardown(const std::string &spec);
};

} // namespace problem

} // namespace suzerain

#endif // SUZERAIN_SIGNAL_DEFINITION_HPP
