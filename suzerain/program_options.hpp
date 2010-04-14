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
 */
class ProgramOptions
{
public:

    /**
     * Default constructor which does not supply a program description.
     */
    ProgramOptions();

    /**
     * Constructor providing a program description.
     *
     * @param description Description to display for <tt>--help</tt> option
     */
    ProgramOptions(std::string &description);

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
    ProgramOptions& add_definition(suzerain::problem::IDefinition &definition);

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
    boost::program_options::options_description_easy_init add_options();

    /**
     * Obtain the underlying <tt>options_description</tt> maintaining
     * the arguments.  Useful to add options specifications more complex
     * than those available through #add_options.
     *
     * @return the underlying <tt>options_description</tt> instance.
     */
    boost::program_options::options_description& options();


    /**
     * Process <tt>main</tt>'s <tt>argc</tt> and <tt>argv</tt> according
     * to the previously provided options and IDefinition instances.
     * This method performs all necessary parsing, file handling,
     * and notification callbacks.  It will throw exceptions on errors.
     * If either <tt>--help</tt> or <tt>--version</tt> is encountered the
     * appropriate message will be written on <tt>std::cout</tt> and
     * <tt>exit(0)</tt> will be called.
     *
     * @param argc Should contain <tt>main</tt>'s \c argc.
     * @param argv Should contain <tt>main</tt>'s \c argv.
     * @param debug Stream on which debugging messages will be output.
     * @param info Stream on which informational messages will be output.
     * @param warn Stream on which warning messages will be output.
     * @param error Stream on which fatal error messages will be output.
     */
    template< typename charT,
              typename DebugOStream,
              typename InfoOStream,
              typename WarnOStream,
              typename ErrorOStream >
    void ProgramOptions::process(int argc,
                                 charT **argv,
                                 DebugOStream &debug,
                                 InfoOStream &info,
                                 WarnOStream &warn,
                                 ErrorOStream &error);

    /**
     * A convenience method suppressing debugging messages and
     * supplying <tt>std::cout</tt> and <tt>std::cerr</tt> to
     * <tt>process</tt>.
     *
     * @param argc Should contain <tt>main</tt>'s \c argc.
     * @param argv Should contain <tt>main</tt>'s \c argv.
     */
    template< typename charT >
    void process(int argc, charT **argv);

    /**
     * Provides access to the variable map used to store options.  Useful if
     * complicated processing or recording a retroactive response file is
     * necessary.  It derives from
     * <tt>std::map<std::string,boost::program_options::variable_value></tt>.
     *
     * @return The Boost.Program_options <tt>variables_map</tt> used by
     *         this instance.
     */
    boost::program_options::variables_map& variables();

protected:

    /**
     * Print version information on the given stream.  Includes
     * the application name, the build date, time and compiler information.
     *
     * @param out Stream on which to write the information.
     * @param application_name Application name to be written.
     */
    template< typename OStream >
    void print_version(OStream&out, const std::string &application_name);

    /**
     * Print help information on the given stream.  Includes
     * the application name, an optional description from construction time,
     * and the full options understood by the application.
     *
     * @param out Stream on which to write the information.
     * @param application_name Application name to be written.
     */
    template< typename OStream >
    void print_help(OStream &out, const std::string &application_name);

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

    /** The application description used for user-oriented messages. */
    std::string application_description_;
};

ProgramOptions::ProgramOptions()
    : variables_(),
      options_("Configuration options"),
      application_description_()
{
    // NOP
}

ProgramOptions::ProgramOptions(std::string &description)
    : variables_(),
      options_("Configuration options"),
      application_description_(description)
{
    // NOP
}

ProgramOptions& ProgramOptions::add_definition(
        suzerain::problem::IDefinition &definition)
{
    options_.add(definition.options());
    return *this;
}

inline
boost::program_options::options_description_easy_init
ProgramOptions::add_options() {
    return options_.add_options();
};

inline
boost::program_options::options_description&
ProgramOptions::options() {
    return options_;
};

template< typename charT,
          typename DebugOStream,
          typename InfoOStream,
          typename WarnOStream,
          typename ErrorOStream >
void ProgramOptions::process(int argc,
                             charT **argv,
                             DebugOStream &debug,
                             InfoOStream &info,
                             WarnOStream &warn,
                             ErrorOStream &error)
{
    typedef std::basic_string<charT> string_charT;
    namespace po = boost::program_options;

    // Prepare options allowed only on command line
    po::options_description desc_clionly("Program information");
    desc_clionly.add_options()
        ("help,h",    "show usage information")
        ("version,v", "print version string")
    ;

    // Prepare options allowed on command line and in configuration file
    // These are never shown to the user
    po::options_description desc_hidden("Hidden options");
    desc_hidden.add_options()
        ("input-file", po::value< std::vector<std::string> >(), "input file")
    ;

    // Build the options acceptable on the CLI, in a file, and in help message
    po::options_description opts_cli;
    opts_cli.add(options_).add(desc_hidden).add(desc_clionly);

    po::options_description opts_file;
    opts_file.add(options_).add(desc_hidden);

    po::options_description opts_visible;
    opts_visible.add(options_).add(desc_clionly);

    // Have positional parameters act like input-file
    po::positional_options_description opts_positional;

    opts_positional.add("input-file", -1);

    // Parse all the command line options
    po::basic_parsed_options<charT> parsed_cli
        = po::basic_command_line_parser<charT>(argc, argv)
                .options(opts_cli)
                .positional(opts_positional)
                .allow_unregistered()
                .run();
    po::store(parsed_cli, variables_);

    // Warn whenever an unrecognized option is encountered on the CLI
    // Warn instead of balk because MPI stacks may utilize argc/argv
    {
        std::vector< string_charT > unrecognized
            = po::collect_unrecognized(parsed_cli.options,
                                       po::exclude_positional);
        BOOST_FOREACH( string_charT option, unrecognized ) {
            warn << "Unrecognized option '" << option << "'"
                 << " on command line"
                 << std::endl;
        }
    }

    // Process command-line only parameters
    if (variables_.count("help")) {
        print_help(std::cout, argv[0]);
        exit(0);
    }
    if (variables_.count("version")) {
        print_version(std::cout, argv[0]);
        exit(0);
    }

    // Parse any input files provided on the command line
    // Earlier files shadow/override settings found in later files
    if (variables_.count("input-file")) {
        using std::vector;
        BOOST_FOREACH(const string_charT &filename,
                      variables_["input-file"].as< vector<string_charT> >()) {

            debug << "Reading additional options from file '"
                  << filename << "'" << std::endl;

            std::basic_ifstream<charT> ifs(
                    filename.c_str(), std::ifstream::in);
            po::basic_parsed_options<charT> parsed_file
                = po::parse_config_file(ifs, opts_file, true);
            po::store(parsed_file, variables_);

            // TODO Display appropriate warning on unrecognized options
            // collect_unrecognized/parse_config_file broken in Boost 1.40
            // Refer to https://svn.boost.org/trac/boost/ticket/3775
            vector< string_charT > unrecognized
                = po::collect_unrecognized(parsed_file.options,
                                           po::exclude_positional);
            BOOST_FOREACH( string_charT option, unrecognized ) {
                warn << "Unrecognized option '" << option << "'"
                     << " in file '" << filename << "'"
                     << std::endl;
            }

            // TODO Display appropriate debug message on shadowed options
        }
    }

    // Perform all notification callbacks because all processing is done
    po::notify(variables_);
}

template< typename charT >
void ProgramOptions::process(int argc, charT **argv)
{
    boost::basic_onullstream<charT> nullstream;
    return process(argc,
                   argv,
                   nullstream, // Debug messages
                   std::cout,  // Info messages
                   std::cerr,  // Warn messages
                   std::cerr); // Error messages
}

inline
boost::program_options::variables_map& ProgramOptions::variables() {
    return variables_;
};

template< typename OStream >
void ProgramOptions::print_version(OStream& out,
                                   const std::string &application_name)
{
#ifdef PACKAGE_STRING
    out << PACKAGE_STRING << " ";
#endif
    out << application_name
        << " (built " __DATE__ " " __TIME__
#if defined(__INTEL_COMPILER)
        << " using Intel "
        << __INTEL_COMPILER << " " << __INTEL_COMPILER_BUILD_DATE
#elif defined(__GNUC__)
        << " using GNU "
        << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__
#else
        << " using an unknown compiler"
#endif
        << ")"
        << std::endl;
}

template< typename OStream >
void ProgramOptions::print_help(OStream &out,
                                const std::string &application_name)
{
    out << "\nUsage: " << application_name << " [OPTION] [FILE]...\n";

    if (!application_description_.empty()) {
        out << application_description_ << '\n';
    }

    out << '\n' << options_ << std::endl;
}

} // namespace suzerain

#endif // __SUZERAIN_PROGRAM_OPTIONS_HPP
