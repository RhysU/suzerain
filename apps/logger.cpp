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
 * logger.cpp: logging tools built atop log4cxx
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <sys/mman.h>
#include <suzerain/error.h>
#include <suzerain/mpi.hpp>
#include <suzerain/os.h>
#include "logger.hpp"
#include <log4cxx/file.h>
#include <log4cxx/helpers/pool.h>
#include <log4cxx/logmanager.h>

#if !defined(LOG4CXX)
#define LOG4CXX 1
#endif
#include <log4cxx/helpers/aprinitializer.h>

namespace logger {

namespace detail {

// Workaround http://old.nabble.com/Static-destruction-fiasco--td31026705.html
// and the related https://issues.apache.org/jira/browse/LOGCXX-338
static struct Log4cxxSingletonBugWorkaroundsType {
    Log4cxxSingletonBugWorkaroundsType() {
        ::log4cxx::helpers::APRInitializer::initialize();
    }
// FIXME Does not correct crashes on exit
//  ~Log4cxxSingletonBugWorkaroundsType() {
//      if (!::log4cxx::helpers::APRInitializer::isDestructed) {
//          ::log4cxx::LogManager::shutdown();
//      }
//  }
} workarounds;

// Global logging instance
::log4cxx::LoggerPtr loggerptr;

// World rank determined within initialize method
int worldrank = 0;

} // end namespace detail

// Provides nice, rank-specific logger names
static void initialize_logger_using_world_rank()
{
    using detail::loggerptr;
    using detail::worldrank;
    using log4cxx::Level;
    using log4cxx::LevelPtr;
    using log4cxx::Logger;
    using suzerain::mpi::comm_rank;
    using suzerain::mpi::comm_rank_identifier;

    loggerptr = Logger::getLogger(comm_rank_identifier(MPI_COMM_WORLD));

    // If debugging is disabled, ensure level is at least INFO on ranks 1+
    LevelPtr level = loggerptr->getLevel();
    if (comm_rank(MPI_COMM_WORLD) > 0 && !loggerptr->isDebugEnabled()) {
        if (!level || !level->isGreaterOrEqual(Level::getInfo())) {
            level = Level::getInfo();
        }
    }
    loggerptr->setLevel(level);
}

// Generate and return a temporary filename template for mkstemp
static char * new_mkstemp_config_template()
{
    static const char affix[] = "/log4cxx.properties.XXXXXX";
    const char *tmpdir = suzerain_temporary_directory();

    char * retval = new char[strlen(tmpdir) + sizeof(affix) + 1];
    strcpy(retval, tmpdir);
    strcat(retval, affix);

    return retval;
}

// Default log4cxx configuration to use when none can be found
static const char default_log4cxx_configuration[] =
    "# See Configuration section at http://logging.apache.org/log4cxx/index.html\n"
    "log4j.rootLogger=DEBUG, CONSOLE\n"
    "log4j.appender.CONSOLE=org.apache.log4j.ConsoleAppender\n"
    "log4j.appender.CONSOLE.layout=org.apache.log4j.PatternLayout\n"
    "log4j.appender.CONSOLE.layout.ConversionPattern=%-5p %8r %c  %m%n\n";

void initialize(MPI_Comm)
{
    // Programmatically enforce execution before log4cxx auto-configuration
    if (log4cxx::LogManager::getLoggerRepository()->isConfigured()) {
        throw std::logic_error(
            "Apache log4cxx subsystem initialized before logger::initialize()."
            "\nThis is a CODING ERROR and will drastically hurt scalability.");
    }

    static const char log4j_config_envvar[]   = "log4j.configuration";
    static const char log4cxx_config_envvar[] = "LOG4CXX_CONFIGURATION";
    using detail::worldrank;
    using detail::loggerptr;

    // Ensure MPI is ready to go
    suzerain::mpi::ensure_mpi_initialized();

    // Store our rank within MPI_COMM_WORLD for INFO0-like macro usage
    worldrank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);

    // Buffer length and data storing the to-be-broadcast configuration
    char *buf    = NULL;  // mmap on rank 0, char[] otherwise
    int   buflen = 0;

    // Initialize rank zero prior to all other ranks so that
    // only rank zero hits the filesystem searching for configuration files.
    if (worldrank == 0) {

        std::string configpath;
        ::log4cxx::helpers::Pool pool;
        ::log4cxx::File file;
        int fd;

        // Lookup logic for configpath based on log4cxx source:
        //   1) Check environment variable LOG4CXX_CONFIGURATION
        //   2) Check environment variable log4j.configuration
        //   3) Check existence of log4cxx.{xml,properties}
        //   4) Check existence of log4j.{xml,properties}
        // See log4cxx/trunk/src/main/cpp/defaultconfigurator.cpp
        if (getenv(log4cxx_config_envvar)) {
            configpath = getenv(log4cxx_config_envvar);
        }
        if (configpath.empty() && getenv(log4j_config_envvar)) {
            configpath = getenv(log4j_config_envvar);
        }
        if (configpath.empty()) {
            if (file.setPath("log4cxx.xml").exists(pool)) {
                configpath = "log4cxx.xml";
            } else if (file.setPath("log4cxx.properties").exists(pool)) {
                configpath = "log4cxx.properties";
            } else if (file.setPath("log4j.xml").exists(pool)) {
                configpath = "log4j.xml";
            } else if (file.setPath("log4j.properties").exists(pool)) {
                configpath = "log4j.properties";
            }
        }

        // Initialize the logger on rank zero and open a file descriptor
        // containing the configuration data to be broadcast
        if (!configpath.empty() && file.setPath(configpath).exists(pool)) {

            // Configuration file was found, use it to initialize logging
            setenv(log4cxx_config_envvar, configpath.c_str(), 1);
            initialize_logger_using_world_rank();

            // Open the configuration file in read-only mode for broadcast
            buflen = boost::numeric_cast<int>(file.length(pool));
            fd = open(configpath.c_str(), O_RDONLY);

            DEBUG("Logging initialized on rank zero using file "
                  << configpath);
        } else {

            // Configuration file was not found, create one.

            // Create a unique temporary configuration file (best effort)
            boost::scoped_array<char> tmpl(new_mkstemp_config_template());
            fd = mkstemp(tmpl.get());
            if (fd != -1) {
                FILE *f = fdopen(fd, "w+");
                if (f) {
                    fwrite(default_log4cxx_configuration,
                           1, sizeof(default_log4cxx_configuration), f);
                    buflen = boost::numeric_cast<int>(ftell(f));
                    fclose(f);
                }
                close(fd);
            }

            // Initialize logging using the temporary file contents
            setenv(log4cxx_config_envvar, tmpl.get(), 1);
            initialize_logger_using_world_rank();
            WARN("Logging system found no configuration file.  Using default:\n"
                 << default_log4cxx_configuration);

            // Re-open the configuration file in read-only mode for broadcast
            fd = open(tmpl.get(), O_RDONLY);
            unlink(tmpl.get());  // preemptive unlink for automatic cleanup

        }

        // mmap configuration file so we can use it as a broadcast buffer
        if (buflen > 0) {
            errno = 0;
            buf = (char *) mmap(NULL, buflen, PROT_READ, MAP_PRIVATE, fd, 0);
            if (buf == MAP_FAILED) {
                WARN("Logging configuration file not broadcast (mmap failed: "
                     << strerror(errno) << ")");
                buflen = 0;
            }
            close(fd);
        }
    }

    // Broadcast buffer length to higher ranks
    SUZERAIN_MPICHKR(MPI_Bcast(&buflen, 1, MPI_INT, 0, MPI_COMM_WORLD));

    // Abort remaining fancy initialization steps if buflen is zero
    if (buflen == 0) {
        if (worldrank == 0) {
            WARN("Scalable logging initialization procedure aborted");
        } else {
            initialize_logger_using_world_rank();
        }
        return;  // Run away!
    }

    // Broadcast configuration file to non-zero ranks
    if (worldrank > 0) buf = new char[buflen];
    SUZERAIN_MPICHKR(MPI_Bcast(buf, buflen, MPI_CHAR, 0, MPI_COMM_WORLD));
    if (worldrank == 0) munmap(buf, buflen);

    // Non-zero ranks write file, initialize logging from it, and remove it
    // High paranoia as either bugs and/or race conditions have arisen herein
    if (worldrank > 0) {

        // Always delete[] previously allocated buf on scope exit
        boost::scoped_array<char> guard_buf(buf);

        // RAII helper to create, close, and unlink a temporary file
        struct guard_tmpfd {
            const char * t;
            int const i;
#ifdef _GNU_SOURCE
            guard_tmpfd(char *t) : t(t), i(mkostemp(t, O_SYNC)) {}
#else
            guard_tmpfd(char *t) : t(t), i(mkstemp(t)) {}
#endif
            operator int()          const { return i; }
            operator const char *() const { return t; }
            operator bool()         const { return i != -1; }
            ~guard_tmpfd() { if (*this) { close(i); unlink(t); } delete[] t; }
        };

        // RAII helper to close a FILE*
        struct guard_pFILE {
            FILE * const p;
            guard_pFILE(FILE *p) : p(p) {}
            operator FILE*() const { return p; }
            operator bool()  const { return p != NULL; }
            ~guard_pFILE() { if (*this) { fclose(p); } }
        };

        // Write the received configuration to a temporary file with known name
        const guard_tmpfd fd(new_mkstemp_config_template());
        if (fd) {
            const guard_pFILE f(fdopen(fd, "w+"));
            if (f) {
                fwrite(buf, 1, buflen, f);
                fflush(f);
            }
        }

        // Initialize log4cxx using the named temporary file
        setenv(log4cxx_config_envvar, fd, 1);
        initialize_logger_using_world_rank();

        // Resources cleaned up by RAII helper class destructors
    }

    // Assert that all ranks indeed have initialized logging subsystems
    assert(log4cxx::LogManager::getLoggerRepository()->isConfigured());
}

} // end namespace logger
