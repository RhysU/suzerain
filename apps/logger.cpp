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
#include <suzerain/error.h>
#include <suzerain/mpi.hpp>

#include "logger.hpp"

#include <log4cxx/basicconfigurator.h>
#include <log4cxx/file.h>
#include <log4cxx/helpers/bytearrayinputstream.h>
#include <log4cxx/helpers/bytebuffer.h>
#include <log4cxx/helpers/exception.h>
#include <log4cxx/helpers/fileinputstream.h>
#include <log4cxx/helpers/pool.h>
#include <log4cxx/helpers/properties.h>
#include <log4cxx/logmanager.h>
#include <log4cxx/appender.h>
#include <log4cxx/spi/filter.h>
#include <log4cxx/spi/loggerrepository.h>
#include <log4cxx/propertyconfigurator.h>

#if !defined(LOG4CXX)
#define LOG4CXX 1
#endif
#include <log4cxx/helpers/aprinitializer.h>

using namespace log4cxx;
using namespace log4cxx::helpers;

// A custom HierarchyEventListener implementation for debugging
class HEL : public virtual spi::HierarchyEventListener,
            public virtual ObjectImpl
{
public:
    DECLARE_LOG4CXX_OBJECT(HEL)
    BEGIN_LOG4CXX_CAST_MAP()
            LOG4CXX_CAST_ENTRY(HEL)
            LOG4CXX_CAST_ENTRY(HierarchyEventListener)
    END_LOG4CXX_CAST_MAP()

    virtual void addAppenderEvent(
        const LoggerPtr &logger, const AppenderPtr &appender)
    {
//      std::cerr << logger->getName() << ' '
//                << appender->getName() << std::endl;
    }

    virtual void removeAppenderEvent(
        const LoggerPtr &logger, const AppenderPtr &appender) { /* NOP */ }
};

IMPLEMENT_LOG4CXX_OBJECT(HEL)

namespace logger {

// Workaround http://old.nabble.com/Static-destruction-fiasco--td31026705.html
static struct Log4cxxWorkaroundsType {
    Log4cxxWorkaroundsType() {
        APRInitializer::initialize();
    }
} workarounds;

LoggerPtr rankzero;

LoggerPtr all;

static spi::HierarchyEventListenerPtr help(new HEL());

static void initialize_loggers_within_comm_world(int worldrank)
{
    // "rankzero" is a disguise for the log4cxx root logger.
    rankzero = Logger::getRootLogger();
    LevelPtr defaultlevel = rankzero->getLevel();
    if (worldrank > 0) rankzero->setLevel(Level::getOff());

    // "all" is a descendent of "rankzero" active on all ranks
    // If not configured, "all" takes rankzero's level from rank 0.
    std::ostringstream oss;
    oss << 'r' << worldrank;
    all = Logger::getLogger(oss.str());
    if (!all->getLevel()) all->setLevel(defaultlevel);
}

// Default log4cxx configuration to use when none can be found
static const char default_log4cxx_config[] =
    "# See \"Configuration\" at http://logging.apache.org/log4cxx/index.html\n"
    "log4j.rootLogger=DEBUG, CONSOLE\n"
    "log4j.appender.CONSOLE=org.apache.log4j.ConsoleAppender\n"
    "log4j.appender.CONSOLE.layout=org.apache.log4j.PatternLayout\n"
    "log4j.appender.CONSOLE.layout.ConversionPattern=%-5p %8r %-10c %m%n\n";

void initialize(MPI_Comm)
{
    static const char log4j_config_envvar[]   = "log4j.configuration";
    static const char log4cxx_config_envvar[] = "LOG4CXX_CONFIGURATION";

    // Programmatically enforce execution before log4cxx auto-configuration.
    if (LogManager::getLoggerRepository()->isConfigured()) {
        throw std::logic_error(
            "Apache log4cxx subsystem initialized before logger::initialize()."
            "\nThis is a CODING ERROR and will drastically hurt scalability.");
    }

    // Ensure MPI is ready to go
    suzerain::mpi::ensure_mpi_initialized();

    // Obtain our rank within MPI_COMM_WORLD
    const int worldrank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);

    // Add our rank-aware HierarchyEventListener to the repository
    LogManager::getLoggerRepository()->addHierarchyEventListener(help);

    // Check that we have not accidentally triggered configuration.
    // Paranoia due to "When the LogManager class is loaded into memory the
    // default initialization procedure is initiated" in LogManager Doxygen.
    if (LogManager::getLoggerRepository()->isConfigured()) {
        throw std::logic_error("log4cxx heisen-configuration detected");
    }

    // Separate logic for rank zero and higher ranks.  Broadcasting of
    // buflen serves to enforce rank zero logic running first.  Only rank
    // zero hits the filesystem searching for configuration files.
    if (worldrank == 0) {

        std::string configpath;
        Pool pool;
        File file;

        // Lookup logic for configpath based on log4cxx source:
        //   1) Check environment variable LOG4CXX_CONFIGURATION
        //   2) Check environment variable log4j.configuration
        //   3) Check existence of log4cxx.properties
        //   4) Check existence of log4j.properties
        // See log4cxx/trunk/src/main/cpp/defaultconfigurator.cpp.  Notice that
        // XML-based files are skipped as log4cxx seems to have no in-memory
        // DOMConfigurator::doConfigure(...) functionality.
        if (getenv(log4cxx_config_envvar)) {
            configpath = getenv(log4cxx_config_envvar);
        }
        if (configpath.empty() && getenv(log4j_config_envvar)) {
            configpath = getenv(log4j_config_envvar);
        }
        if (configpath.empty()) {
            if (file.setPath("log4cxx.properties").exists(pool)) {
                configpath = "log4cxx.properties";
            } else if (file.setPath("log4j.properties").exists(pool)) {
                configpath = "log4j.properties";
            }
        }

        // Overwrite default configuration using file contents, if possible.
        std::vector<unsigned char> buf;
        buf.assign(default_log4cxx_config,
                   default_log4cxx_config + sizeof(default_log4cxx_config));
        if (!configpath.empty() && file.setPath(configpath).exists(pool)) {
            // Hideous ifdef to workaround log4cxx's inconsistent char usage.
            buf.resize(file.length(pool));
            FileInputStream fis(file);
#if CHAR_MIN == 0  /* "char" denotes "unsigned char" */
            ByteBuffer bb(reinterpret_cast<char *>(&buf.front()), buf.size());
            fis.read(bb);
#else              /* "char" denotes "signed char" */
            std::vector<char> tmp(buf.size());
            ByteBuffer bb(&tmp.front(), tmp.size());
            fis.read(bb);
            // Not the least bit ASCII 128+ safe but throws if issue arises.
            std::transform(tmp.begin(), tmp.end(), buf.begin(),
                           boost::numeric::converter<unsigned char,char>());
#endif
        }

        // Configure logging subsystem using the buffer contents
        InputStreamPtr bais(new ByteArrayInputStream(buf));
        Properties props;
        PropertyConfigurator pc;
        try {
            props.load(bais);
            pc.doConfigure(props, LogManager::getLoggerRepository());
        } catch (Exception &e) {
            BasicConfigurator::configure(); // Fallback
            buf.push_back('\0');
            WARN("Logging configuration from " << configpath << "invalid: "
                 << e.what() << "\n" <<  &buf.front());
            buf.pop_back();
        }
        initialize_loggers_within_comm_world(worldrank);

        // Broadcast configuration buffer length and contents to other ranks
        int buflen = buf.size();
        SUZERAIN_MPICHKR(MPI_Bcast(&buflen, 1, MPI_INT, 0, MPI_COMM_WORLD));
        SUZERAIN_MPICHKR(MPI_Bcast(&buf.front(), buflen,
                                   MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD));

    } else {

        // Higher ranks receive buffer length, buffer, and initialize logging
        int buflen;
        SUZERAIN_MPICHKR(MPI_Bcast(&buflen, 1, MPI_INT, 0, MPI_COMM_WORLD));
        std::vector<unsigned char> buf(buflen);
        SUZERAIN_MPICHKR(MPI_Bcast(&buf.front(), buflen,
                                   MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD));
        InputStreamPtr bais(new ByteArrayInputStream(buf));
        Properties props;
        PropertyConfigurator pc;
        try {
            props.load(bais);
            pc.doConfigure(props, LogManager::getLoggerRepository());
        } catch (Exception & /* ignored */) {
            BasicConfigurator::configure(); // Fallback
        }
        initialize_loggers_within_comm_world(worldrank);

    }

    // Assert that all ranks indeed have initialized logging
    assert(LogManager::getLoggerRepository()->isConfigured());
}

} // end namespace logger
