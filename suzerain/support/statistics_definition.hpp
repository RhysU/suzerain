//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
class statistics_definition : public definition_base
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
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      \c desttemplate and \c retain.
     */
    explicit statistics_definition(
            const std::string& destination            ,
            const std::size_t  retain      = (1 << 15),
            const real_t       dt          = 0        ,
            const std::size_t  nt          = 0        );

    /** @copydoc support::definition_base::options_description() */
    boost::program_options::options_description options_description();

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
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_STATISTICS_DEFINITION_HPP
