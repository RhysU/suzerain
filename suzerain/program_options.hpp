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
 * program_options.hpp: handles parsing program options from CLI, files
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_PROGRAM_OPTIONS_HPP
#define __SUZERAIN_PROGRAM_OPTIONS_HPP

#include <suzerain/common.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/problem.hpp>

/** @file
 * Provides classes handling program options parsing given a problem
 * definition.
 */

namespace suzerain {

/**
 * Performs options parsing for an application, including command line and
 * options file processing.  Includes the following functionality:
 *
 *   - Allows automatically adding options given an IDefinition instance
 *   - Allows adding options using
 *     <tt>boost::program_options::options_description_easy_init</tt>
 *   - Provides uniform processing of both command line arguments and files
 *   - Command line options take precedence over file-based options
 *   - Earlier file options take precedence over later file options
 *   - Warns whenever an unrecognized option is provided
 *   - Provides <tt>--help</tt> and <tt>--version</tt> functionality
 *   - Optionally provides a program description for <tt>--help</tt>
 *   - Reads any input files in a scalable way by reading them on
 *     the root process and then broadcasting the contents to all ranks.
 */
class ProgramOptions
{
public:

    /**
     * Default constructor which does not supply a program description.
     */
    ProgramOptions() : variables_(), options_() {}

    /**
     * Constructor providing a program synopsis, an argument synopsis,
     * and a more complete application description.
     *
     * @param application_synopsis Application synopsis for <tt>--help</tt>.
     * @param argument_synopsis    Argument synopsis to display for
     *                             <tt>--help</tt> option.  For example,
     *                             "[FILE]..." or "SOURCE... DIRECTORY".
     * @param description          Extent description of the application to
     *                             be displayed at the bottom of
     *                             <tt>--help</tt>.
     */
    ProgramOptions(const std::string &application_synopsis,
                   const std::string &argument_synopsis = "",
                   const std::string &description = "")
        : variables_(),
          options_(),
          application_synopsis_(application_synopsis),
          argument_synopsis_(argument_synopsis),
          application_description_(description) {}

    /**
     * Adds all the options stored within the given IDefinition to the
     * program's known command line options.  Any notification callbacks
     * defined therein will be called after options processing completes.
     *
     * @param definition Definition instance that will receive any
     *        resulting information from program options processing.
     *
     * @return <tt>*this</tt> to support chaining.
     *
     * @see suzerain::problem::IDefinition for the necessary contract.
     */
    ProgramOptions& add_definition(suzerain::problem::IDefinition &definition)
    {
        options_.add(definition.options());
        return *this;
    }

    /**
     * Allows callers to add new options using the Boost.Program_options
     * <tt>options_description_easy_init</tt> mechanism.
     *
     * @return A "live" <tt>options_description_easy_init</tt> instance
     *         tied to this instance.
     *
     * @see The <a href="http://www.boost.org/doc/html/program_options.html">
     *      Boost.Program_options documention</a> for more details.
     */
    boost::program_options::options_description_easy_init add_options()
    {
        return options_.add_options();
    }

    /**
     * Obtain the underlying <tt>options_description</tt> maintaining
     * the arguments.  Useful to add options specifications more complex
     * than those available through #add_options.
     *
     * @return the underlying <tt>options_description</tt> instance.
     */
    boost::program_options::options_description& options() { return options_; }

    /**
     * Process <tt>main</tt>'s <tt>argc</tt> and <tt>argv</tt> according to the
     * previously provided options and IDefinition instances.  This method
     * performs all necessary parsing, file handling, and notification
     * callbacks.  It will throw exceptions on errors.  If either
     * <tt>--help</tt> or <tt>--version</tt> is encountered the appropriate
     * message will be written on <tt>std::cout</tt> and <tt>exit(0)</tt> will
     * be called.  Any input files are read by only rank zero of <tt>comm</tt>
     * and are then broadcast to all other ranks.  Error messages are only
     * output by rank zero.
     *
     * @param argc  Should contain <tt>main</tt>'s \c argc.
     * @param argv  Should contain <tt>main</tt>'s \c argv.
     * @param comm  MPI Communicator over which the processing is collective.
     * @param debug Stream on which debugging messages will be output.
     * @param info  Stream on which informational messages will be output.
     * @param warn  Stream on which warning messages will be output.
     * @param error Stream on which fatal error messages will be output.
     *
     * @return A vector containing all non-option arguments.
     */
    std::vector<std::string> process(int argc,
                                     char **argv,
                                     MPI_Comm comm,
                                     std::ostream &debug,
                                     std::ostream &info,
                                     std::ostream &warn,
                                     std::ostream &error);

    /**
     * A convenience method suppressing debugging messages and supplying
     * <tt>std::cout</tt> and <tt>std::cerr</tt> to <tt>process</tt>.
     *
     * @param argc Should contain <tt>main</tt>'s \c argc.
     * @param argv Should contain <tt>main</tt>'s \c argv.
     * @param comm MPI Communicator over which the processing is collective.
     *
     * @return A vector containing all non-option arguments.
     */
    std::vector<std::string> process(int argc,
                                     char **argv,
                                     MPI_Comm comm = MPI_COMM_WORLD)
    {
        boost::onullstream nullstream;
        return process(argc,
                       argv,
                       comm,
                       nullstream, // Debug messages
                       std::cout,  // Info messages
                       std::cerr,  // Warn messages
                       std::cerr); // Error messages
    }

    /**
     * Provides access to the variable map used to store options.  Useful if
     * complicated processing or recording a retroactive response file is
     * necessary.  It derives from
     * <tt>std::map<std::string,boost::program_options::variable_value></tt>.
     *
     * @return The Boost.Program_options <tt>variables_map</tt> used by
     *         this instance.
     */
    boost::program_options::variables_map& variables() { return variables_; }

    /** @return the "--verbose" flag count from the command line.  */
    int verbose() { return verbose_; }

    /** @return the "--verbose-all" flag count from the command line.  */
    int verbose_all() { return verbose_all_; }

protected:

    /**
     * Protected implementation of processing logic.  Implementation kept
     * separate from public API to simplify using of
     * <tt>boost::onullstream</tt> on non-zero ranks.
     *
     * @see process for the public API.
     */
    std::vector<std::string> process_internal(int argc,
                                              char **argv,
                                              MPI_Comm comm,
                                              std::ostream &debug,
                                              std::ostream &info,
                                              std::ostream &warn,
                                              std::ostream &error);

    /**
     * The Boost.Program_options variables_map used internally.
     * It is populated during #process invocation.
     */
    boost::program_options::variables_map variables_;

    /**
     * The Boost.Program_options options_description used internally.
     * It should be populated prior to #process invocation.
     */
    boost::program_options::options_description options_;

    /** The application synopsis used for user-oriented messages. */
    std::string application_synopsis_;

    /**
     * The synopsis used for arguments in the usage message.
     * For examples, look at the <tt>cp(1)</tt> SYNOPSIS section
     * after the string "[OPTION]...".
     */
    std::string argument_synopsis_;

    /** The application description used for user-oriented messages. */
    std::string application_description_;

    /** Maintain the "--verbose" flag count from the command line */
    int verbose_;

    /** Maintain the "--verbose-all" flag count from the command line */
    int verbose_all_;
};

} // namespace suzerain

#endif // __SUZERAIN_PROGRAM_OPTIONS_HPP
