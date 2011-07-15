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
 * restart_definition.hpp: classes handling restart definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_RESTART_DEFINITION_HPP
#define __SUZERAIN_RESTART_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling restart definitions, which are runtime
 * arguments used to control simulation restart behavior.
 */

namespace suzerain {

namespace problem {

/**
 * Encapsulates flags related to restart file behavior for the simulation.
 * Includes the particulars around writing and archiving restart files.
 */
template< class String = std::string >
class RestartDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default values.
     * All of these can be overridden by command line options.
     *
     * @param default_metadata     Restart file path to use when saving
     *                             common restart file metadata.
     * @param default_uncommitted  Restart file path to use when saving
     *                             uncommitted restart data.
     * @param default_desttemplate Restart archiving pattern to use when
     *                             committing restart files.
     * @param default_retain       Maximum number of committed restart files
     *                             to retain.
     * @param default_restart_nt   Number of simulation steps between restart
     *                             writes.
     * @param default_restart_dt   Amount of simulation time between restart
     *                             writes.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      default_desttemplate and default_retain.
     */
    explicit RestartDefinition(
            const String& default_metadata     = "",
            const String& default_uncommitted  = "",
            const String& default_desttemplate = "",
            int default_retain                 = 1,
            double default_restart_dt          = 0,
            int default_restart_nt             = 0);

    /**
     * The restart file path to use when saving common restart file
     * metadata.
     *
     * @return the metadata file path to use.
     */
    String& metadata() { return metadata_; }

    /** @copydoc RestartDefinition::metadata() */
    const String& metadata() const { return metadata_; }

    /**
     * The file path to use when saving uncommitted restart data.
     *
     * @return the file path to use when saving uncommitted restart data.
     */
    String& uncommitted() { return uncommitted_; }

    /** @copydoc RestartDefinition::uncommitted() */
    const String& uncommitted() const { return uncommitted_; }

    /**
     * The archiving pattern to use when committing restart files.
     *
     * @return the archiving pattern to use when committing restart files.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      desttemplate.
     */
    String& desttemplate() { return desttemplate_; }

    /** @copydoc RestartDefinition::desttemplate() */
    const String& desttemplate() const { return desttemplate_; }

    /**
     * The maximum number of committed restart files to retain.
     *
     * @return the maximum number of committed restart files to retain.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      retain.
     */
    int& retain() { return retain_; }

    /** @copydoc RestartDefinition::retain() */
    const int& retain() const { return retain_; }

    /**
     * The maximum amount of simulation time between writing restart
     * files.
     *
     * @return the maximum amount of simulation time between writing restart
     * files.
     */
    double& restart_dt() { return restart_dt_; }

    /** @copydoc RestartDefinition::restart_dt() */
    const double& restart_dt() const { return restart_dt_; }

    /**
     * The maximum number of simulation steps to take between
     * writing restart files.
     *
     * @return the maximum number of simulation steps to take between
     * writing restart files.
     */
    int& restart_nt() { return restart_nt_; }

    /** @copydoc RestartDefinition::restart_nt() */
    const int& restart_nt() const { return restart_nt_; }

private:

    /** Stores the metadata write path */
    String metadata_;

    /** Stores the uncommitted restart path */
    String uncommitted_;

    /** Stores the committed restart path */
    String desttemplate_;

    /** Stores the maximum committed retain count */
    int retain_;

    /** Stores restart writing frequency in simulation time */
    double restart_dt_;

    /** Stores restart writing frequency in time steps */
    int restart_nt_;
};

template< class String >
RestartDefinition<String>::RestartDefinition(
        const String& default_metadata,
        const String& default_uncommitted,
        const String& default_desttemplate,
        int default_retain,
        double default_restart_dt,
        int default_restart_nt)
    : IDefinition("Restart-related parameters"),
      metadata_(default_metadata),
      uncommitted_(default_uncommitted),
      desttemplate_(default_desttemplate),
      retain_(default_retain),
      restart_dt_(default_restart_dt),
      restart_nt_(default_restart_nt)
{
    namespace po = ::boost::program_options;

    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_positive;
    using ::suzerain::validation::ensure_nonnegative;

    ::std::pointer_to_binary_function<int,const char*,void>
        ptr_fun_ensure_positive_int(ensure_positive<int>);
    ::std::pointer_to_binary_function<int,const char*,void>
        ptr_fun_ensure_nonnegative_int(ensure_nonnegative<int>);
    ::std::pointer_to_binary_function<double,const char*,void>
        ptr_fun_ensure_nonnegative_double(ensure_nonnegative<double>);

    this->add_options()
        ("metadata", po::value<String>(&metadata_)
            ->default_value(metadata_),
         "Path to use when saving common restart metadata")
        ("uncommitted", po::value<String>(&uncommitted_)
            ->default_value(uncommitted_),
         "Path to use when saving uncommitted restart data")
        ("desttemplate", po::value<String>(&desttemplate_)
            ->default_value(desttemplate_),
         "Restart archiving pattern to use when committing restart files")
        ("retain", po::value(&retain_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_int, "retain"))
            ->default_value(retain_),
         "Maximum number of committed restart files to retain")
        ("restart_dt", po::value(&restart_dt_)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_double, "restart_dt"))
            ->default_value(restart_dt_),
         "Maximum amount of simulation time between restart files")
        ("restart_nt", po::value(&restart_nt_)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_int, "restart_nt"))
            ->default_value(restart_nt_),
         "Maximum number of time steps between restart files")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_RESTART_DEFINITION_HPP
