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

#ifndef SUZERAIN_SUPPORT_PROGRAM_OPTIONS_HPP
#define SUZERAIN_SUPPORT_PROGRAM_OPTIONS_HPP

/** @file
 * Classes handling program options parsing.
 */

#include <suzerain/common.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/support/definition_base.hpp>

// https://svn.boost.org/trac/boost/ticket/7568
SUZERAIN_GCC_DIAG_OFF(unused-parameter);
#include <boost/program_options.hpp>
SUZERAIN_GCC_DIAG_ON(unused-parameter);

namespace suzerain {

namespace support {

/**
 * Performs options parsing for an application, including command line and
 * options file processing.  Includes the following functionality:
 *
 *   - Allows automatically adding options given an definition_base instance
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
class program_options
{
public:

    /**
     * Default constructor which does not supply a program description.
     */
    program_options();

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
     * @param version              Version information to be displayed
     *                             when <tt>--version</tt> is used.
     */
    program_options(const std::string &application_synopsis,
                    const std::string &argument_synopsis = "",
                    const std::string &description = "",
                    const std::string &version = "");

    /**
     * Adds all the options stored within the given definition_base to the
     * program's known command line options.  Any notification callbacks
     * defined therein will be called after options processing completes.
     *
     * @param definition Definition instance that will receive any
     *        resulting information from program options processing.
     *
     * @return <tt>*this</tt> to support chaining.
     *
     * @see problem::definition_base for the necessary contract.
     */
    program_options& add_definition(definition_base &definition);

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
    boost::program_options::options_description& options()
    {
        return options_;
    }

    /**
     * Process <tt>main</tt>'s <tt>argc</tt> and <tt>argv</tt> according to the
     * previously provided options and definition_base instances.  This method
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
                                     MPI_Comm comm = MPI_COMM_WORLD);

    /**
     * Provides access to the variable map used to store options.  Useful if
     * complicated processing or recording a retroactive response file is
     * necessary.  It derives from
     * <tt>std::map<std::string,boost::program_options::variable_value></tt>.
     *
     * @return The Boost.Program_options <tt>variables_map</tt> used by
     *         this instance.
     */
    boost::program_options::variables_map& variables()
    {
        return variables_;
    }

    /**
     * Provides constant access to the variable map used to store options.
     *
     * @return The Boost.Program_options <tt>variables_map</tt> used by
     *         this instance.
     */
    const boost::program_options::variables_map& variables() const
    {
        return variables_;
    }

    /** @return the "--verbose" flag count from the command line.  */
    int verbose()
    {
        return verbose_;
    }

    /** @return the "--verbose-all" flag count from the command line.  */
    int verbose_all() {
        return verbose_all_;
    }

    /**
     * Ensures that \c opt1 and \c opt2 were not both specified.
     * Should only be invoked after process().
     *
     * @param opt1 First option name to check
     * @param opt2 Second option name to check
     *
     * @throw std::invalid_argument if both opt1 and opt2 were supplied.
     */
    template<typename String>
    void conflicting_options(const String opt1, const String opt2) const
    {
        // Routine is Copyright Vladimir Prus 2002-2004.
        // Distributed under the Boost Software License, Version 1.0.
        // (See http://www.boost.org/LICENSE_1_0.txt)
        if (   variables_.count(opt1) && !variables_[opt1].defaulted()
            && variables_.count(opt2) && !variables_[opt2].defaulted())
            throw std::invalid_argument(std::string("Conflicting options '")
                                        + opt1 + "' and '" + opt2 + "'.");
    }

    /**
     * Ensures that if \c for_what was specified then \c required_option was
     * specified too.  Should only be invoked after process().
     *
     * @param for_what First option name to check
     * @param required_option Second option name to check
     *
     * @throw std::invalid_argument if \c for_what was specified without \c
     * required_option.
     */
    template<typename String>
    void option_dependency(const String for_what,
                           const String required_option) const
    {
        // Routine is Copyright Vladimir Prus 2002-2004.
        // Distributed under the Boost Software License, Version 1.0.
        // (See http://www.boost.org/LICENSE_1_0.txt)
        if (variables_.count(for_what) && !variables_[for_what].defaulted())
            if (   variables_.count(required_option) == 0
                || variables_[required_option].defaulted())
                throw std::invalid_argument(std::string("Option '") + for_what
                            + "' requires option '" + required_option + "'.");
    }

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

    /** The application version used for user-oriented messages. */
    std::string application_version_;

    /** Maintain the "--verbose" flag count from the command line */
    int verbose_;

    /** Maintain the "--verbose-all" flag count from the command line */
    int verbose_all_;
};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_PROGRAM_OPTIONS_HPP
