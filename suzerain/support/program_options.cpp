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

/** @file
 * @copydoc program_options.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/support/program_options.hpp>

#include <cerrno>
#include <cstdlib>
#include <boost/test/utils/nullstream.hpp>
#include <libgen.h>
#include <limits.h>

#include <suzerain/version.hpp>

// Parse Unix-like verbosity flags (http://stackoverflow.com/questions/5486753)
static std::pair<std::string, std::string> verbosity(const std::string& s)
{
    std::pair<std::string, std::string> retval;

    if (s.find("-v") == 0) {

        size_t value = 1;
        try {
            value = boost::lexical_cast<size_t>(s.substr(2));
        } catch(...) {
            while (s[1+value] == 'v') ++value;
        }
        retval.first = "verbose";
        retval.second = boost::lexical_cast<std::string>(value);

    } else if (s.find("-V") == 0) {

        size_t value = 1;
        try {
            value = boost::lexical_cast<size_t>(s.substr(2));
        } catch(...) {
            while (s[1+value] == 'V') ++value;
        }
        retval.first = "verbose-all";
        retval.second = boost::lexical_cast<std::string>(value);

    }

    return retval;
}

namespace suzerain {

namespace support {

program_options::program_options()
    : variables_()
    , options_()
    , application_synopsis_()
    , argument_synopsis_()
    , application_description_()
    , application_version_()
    , verbose_()
    , verbose_all_()
{
}

program_options::program_options(
        const std::string &application_synopsis,
        const std::string &argument_synopsis,
        const std::string &description,
        const std::string &version)
    : variables_()
    , options_()
    , application_synopsis_(application_synopsis)
    , argument_synopsis_(argument_synopsis)
    , application_description_(description)
    , application_version_(version)
    , verbose_()
    , verbose_all_()
{
}

program_options& program_options::add_definition(
        definition_base &definition)
{
    options_.add(definition.options_description());
    return *this;
}

std::vector<std::string> program_options::process(
        int argc,
        char **argv,
        MPI_Comm comm,
        std::ostream &debug,
        std::ostream &info,
        std::ostream &warn,
        std::ostream &error)
{
    const int rank = mpi::comm_rank(comm);
    if (rank == 0) {
        return process_internal(argc, argv, comm,
                                debug, info, warn, error);
    } else {
        boost::onullstream nullstream;
        return process_internal(argc, argv, comm,
                nullstream, nullstream, nullstream, nullstream);
    }
}

std::vector<std::string> program_options::process(
        int argc,
        char **argv,
        MPI_Comm comm)
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

std::vector<std::string> program_options::process_internal(
        int argc, char **argv, MPI_Comm comm,
        std::ostream &debug, std::ostream &info,
        std::ostream &warn, std::ostream &error)
{
    SUZERAIN_UNUSED(debug); // Possibly unused
    SUZERAIN_UNUSED(info);  // Possibly unused
    SUZERAIN_UNUSED(warn);  // Possibly unused
    SUZERAIN_UNUSED(error); // Possibly unused
    using std::endl;
    using std::istream;
    using std::streambuf;
    using std::string;
    using std::vector;
    namespace po = boost::program_options;
    const int rank = mpi::comm_rank(comm);

    // Obtain the program name using basename(argv[0])
    std::string program_name;
    {
        shared_ptr<char> name_copy(strdup(argv[0]), &free);
        program_name = basename(name_copy.get());
    }

    // Prepare options allowed only on command line
    po::options_description desc_clionly("Program information");
    desc_clionly.add_options()
        ("help,h",
         "Show usage information")
        ("version",
         "Print version string")
        ("verbose,v", po::value<std::vector<std::string> >(),
         "Increase verbosity of messages from rank zero")
        ("verbose-all,V", po::value<std::vector<std::string> >(),
         "Increase verbosity of all rank-specific messages")
    ;

    // Prepare a response-file option iff non-trivial options available
    if (options_.options().size() > 0) {
        desc_clionly.add_options()
            ("response-file", po::value< vector<string> >()->composing(),
             "File to additionally read for options")
        ;
    }

    // Prepare options allowed on command line and in configuration file
    // These are never shown to the user.
    // TODO Do we want __positional__ available in a response file?  Eh...
    po::options_description desc_hidden("Hidden options");
    desc_hidden.add_options()
        ("__positional__", po::value< vector<string> >(), "positional args")
    ;

    // Guess at the option width to use (krufty and bash-ish, currently).
    // FIXME Make terminal-width detection work as expected.
    // Problem is not in the code below, but in how linewidth is used
    // by Boost.  See http://lists.boost.org/boost-users/2014/01/81007.php.
    unsigned linewidth = po::options_description::m_default_line_length;
    {
        using namespace std;
        if (getenv("COLUMNS")) {
            errno = 0;
            const long t = strtol(getenv("COLUMNS"), (char **) NULL, 10);
            if (    !errno
                && t >= numeric_limits<unsigned>::min()
                && t <= numeric_limits<unsigned>::max()) {
                linewidth = t;
            }
        }
    }

    // Build the options acceptable on the CLI, in a file, and in help message
    po::options_description opts_cli    (linewidth);
    po::options_description opts_file   (linewidth);
    po::options_description opts_visible(linewidth);
    opts_cli    .add(options_).add(desc_hidden).add(desc_clionly);
    opts_file   .add(options_).add(desc_hidden)                  ;
    opts_visible.add(options_)                 .add(desc_clionly);

    // Collect positional arguments into "__positional__"
    po::positional_options_description opts_positional;
    opts_positional.add("__positional__", -1);

    // Parse all the command line options
    po::parsed_options parsed_cli = po::command_line_parser(argc, argv)
                                        .options(opts_cli)
                                        .positional(opts_positional)
                                        .allow_unregistered()
                                        .extra_parser(verbosity)
                                        .run();
    po::store(parsed_cli, variables_);

    // Warn whenever an unrecognized option is encountered on the CLI
    // Warn instead of balk because MPI stacks may utilize argc/argv
    {
        vector<string> unrecognized = po::collect_unrecognized(
                parsed_cli.options, po::exclude_positional);
        BOOST_FOREACH( string option, unrecognized ) {
            warn << "Unrecognized option '" << option << "'"
                 << " on command line"
                 << endl;
        }
    }

    // Process --help
    if (variables_.count("help")) {
        info << "Usage: " << program_name
             << " [OPTION]... " << argument_synopsis_
             << '\n';

        if (!application_synopsis_.empty()) {
            info << application_synopsis_ << '\n';
        }

        if (opts_visible.options().size() > 0) {
            if (options_.options().size() > 0) {
                info << '\n' << "Options:" << '\n';
            }
            info << opts_visible << '\n';
        }

        if (!application_description_.empty()) {
            info << '\n' << application_description_ << '\n';
        }

        exit(0);
    }

    // Process --version using suzerain::version(...)
    if (variables_.count("version")) {
        info << version(program_name, application_version_) << std::endl;
        exit(0);
    }

    // Parse any input files provided on the command line
    // Earlier files shadow/override settings found in later files
    if (variables_.count("response-file")) {
        BOOST_FOREACH(const string &filename,
                      variables_["response-file"].as< vector<string> >()) {

            debug << "Reading additional options from file '"
                  << filename << "'" << endl;

            // Specialized istream::failure subclass for this routine
            struct failure : public istream::failure {
                failure(const string &op, const string& file, int errnum)
                    : istream::failure("Failure during " + op
                            + " for response-file '" + file
                            + "': " + strerror(errnum) ) {}

                failure(const string& file) : istream::failure(
                        "Failure processing response-file '" + file + "'") {}
            };

            long len;
            scoped_array<char> buf;

            // Rank zero slurps the file into appropriately-sized buf
            if (rank == 0) {
                errno = 0;
                FILE * fp = NULL;
                try {
                    // Open file and obtain length
                    fp = fopen(filename.c_str(), "rb");
                    if (!fp || errno)
                        throw failure("fopen", filename, errno);
                    if (fseek(fp, 0, SEEK_END))
                        throw failure("fseek", filename, errno);
                    len = ftell(fp);
                    if (len == -1)
                        throw failure("ftell", filename, errno);
                    if (fseek(fp, 0, SEEK_SET))
                        throw failure("fseek", filename, errno);

                    // Allocate buffer and slurp contents
                    buf.reset(new char[len]);
                    if (fread(buf.get(), 1, len, fp) != (size_t) len) {
                        throw failure("fread", filename, errno);
                    }

                    // Close appears after catch
                } catch (istream::failure &f) {
                    error << f.what() << endl;
                    len = -1;
                }
                if (fp) fclose(fp);
            }

            // Buffer size and then buf contents are broadcast to other ranks
            MPI_Bcast(&len, 1, MPI_LONG, 0, comm);
            if (len == -1)
                throw failure(filename);
            if (rank > 0)
                buf.reset(new char[len]);
            MPI_Bcast(buf.get(), len, MPI_CHAR, 0, MPI_COMM_WORLD);

            // Create an istream from the buffer contents a la
            // bytes.com/topic/c/answers/582365-making-istream-char-array
            struct membuf : public streambuf
            {
                membuf(streambuf::char_type *b, size_t n)
                {
                    this->setg(b, b, b + n);
                }
            };
            membuf mb(buf.get(), len);
            istream is(&mb);

            // Parse and store from the buffer-based istream
            po::parsed_options parsed_file
                = po::parse_config_file(is, opts_file, true);
            po::store(parsed_file, variables_);

            // TODO Display appropriate warning on unrecognized options
            // collect_unrecognized/parse_config_file broken in Boost 1.40
            // Refer to https://svn.boost.org/trac/boost/ticket/3775
            {
                vector< string > unrecognized = po::collect_unrecognized(
                        parsed_file.options, po::exclude_positional);
                BOOST_FOREACH( string option, unrecognized ) {
                    warn << "Unrecognized option '" << option << "'"
                        << " in response-file '" << filename << "'"
                        << endl;
                }
            }

            // TODO Display appropriate debug message on shadowed options
        }
    }

    // Compute --verbose and --verbose-all levels by summation
    // following http://stackoverflow.com/questions/5486753
    verbose_ = 0;
    if (variables_.count("verbose")) {
        const std::vector<std::string> &values
            = variables_["verbose"].as<std::vector<std::string> >();
        for (std::size_t i = 0; i < values.size(); ++i) {
            verbose_ += boost::lexical_cast<int>(values[i]);
        }
    }
    verbose_all_ = 0;
    if (variables_.count("verbose-all")) {
        const std::vector<std::string> &values
            = variables_["verbose-all"].as<std::vector<std::string> >();
        for (std::size_t i = 0; i < values.size(); ++i) {
            verbose_all_ += boost::lexical_cast<int>(values[i]);
        }
    }

    // Perform all notification callbacks because all processing is done
    po::notify(variables_);

    // Return all position arguments to the caller for his/her use.  Run into
    // any_cast issues if returning non-existent values, so be a bit careful
    // about what exactly we return.
    if (variables_.count("__positional__")) {
        return variables_["__positional__"].as< vector<string> >();
    } else {
        return vector<string>();
    }
}

} // end namespace support

} // end namespace suzerain
