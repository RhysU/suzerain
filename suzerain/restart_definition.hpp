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
 * Includes whether or not to load a restart file and the particulars around
 * writing and archiving restart files.
 */
template< class String = std::string >
class RestartDefinition : public IDefinition
{
public:
    /**
     * Construct an instance with the given default values.
     * All of these can be overridden by command line options.
     *
     * @param default_load         Restart file path to load at startup.
     *                             Leave zero-length to specify new simulation.
     * @param default_metadata     Restart file path to use when saving
     *                             common restart file metadata.
     * @param default_uncommitted  Restart file path to use when saving
     *                             uncommitted restart data.
     * @param default_desttemplate Restart archiving pattern to use when
     *                             committing restart files.
     * @param default_retain       Maximum number of committed restart files
     *                             to retain.
     * @param default_every_nt     Number of simulation steps between restart
     *                             writes.
     * @param default_every_dt     Amount of simulation time between restart
     *                             writes.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      default_desttemplate and default_retain.
     */
    explicit RestartDefinition(
            const String& default_load         = "",
            const String& default_metadata     = "",
            const String& default_uncommitted  = "",
            const String& default_desttemplate = "",
            int default_retain                 = 1,
            double default_every_dt            = 0,
            int default_every_nt               = 0);

    /**
     * Retrieve the restart file path to load on simulation startup.  A
     * zero-length value indicates that no restart file should be loaded.
     *
     * @return the restart file path to load.
     */
    const String& load() const { return load_; }

    /**
     * Retrieve the restart file path to use when saving common restart file
     * metadata.
     *
     * @return the metadata file path to use.
     */
    const String& metadata() const { return metadata_; }

    /**
     * Retrieve the file path to use when saving uncommitted restart data.
     *
     * @return the file path to use when saving uncommitted restart data.
     */
    const String& uncommitted() const { return uncommitted_; }

    /**
     * Retrieve the archiving pattern to use when committing restart files.
     *
     * @return the archiving pattern to use when committing restart files.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      desttemplate.
     */
    const String& desttemplate() const { return desttemplate_; }

    /**
     * Retrieve the maximum number of committed restart files to retain.
     *
     * @return the maximum number of committed restart files to retain.
     *
     * @see ESIO's esio_file_close_restart() for the semantics of
     *      retain.
     */
    int retain() const { return retain_; }

    /**
     * Retrieve the maximum amount of simulation time between writing restart
     * files.
     *
     * @return the maximum amount of simulation time between writing restart
     * files.
     */
    double every_dt() const { return every_dt_; }

    /**
     * Retrieve the maximum number of simulation steps to take between
     * writing restart files.
     *
     * @return the maximum number of simulation steps to take between
     * writing restart files.
     */
    int every_nt() const { return every_nt_; }

private:

    String load_;         /**< Stores the restart load path */
    String metadata_;     /**< Stores the metadata write path */
    String uncommitted_;  /**< Stores the uncommitted restart path */
    String desttemplate_; /**< Stores the committed restart path */
    int retain_;          /**< Stores the maximum committed retain */
    double every_dt_;      /**< Stores restart writing frequency in simulation time */
    int every_nt_;        /**< Stores restart writing frequency in time steps */
};

template< class String >
RestartDefinition<String>::RestartDefinition(
        const String& default_load,
        const String& default_metadata,
        const String& default_uncommitted,
        const String& default_desttemplate,
        int default_retain,
        double default_every_dt,
        int default_every_nt)
    : IDefinition("Restart-related parameters"),
      load_(default_load),
      metadata_(default_metadata),
      uncommitted_(default_uncommitted),
      desttemplate_(default_desttemplate),
      retain_(default_retain),
      every_dt_(default_every_dt),
      every_nt_(default_every_nt)
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
        ("load", po::value<String>(&load_)
            ->default_value(load_),
         "Restart file to load on startup")
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
        ("every_dt", po::value(&every_dt_)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_double, "every_dt"))
            ->default_value(every_dt_),
         "Maximum amount of simulation time between restart files")
        ("every_nt", po::value(&every_nt_)
            ->notifier(bind2nd(ptr_fun_ensure_nonnegative_int, "every_nt"))
            ->default_value(every_nt_),
         "Maximum number of time steps between restart files")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_RESTART_DEFINITION_HPP
