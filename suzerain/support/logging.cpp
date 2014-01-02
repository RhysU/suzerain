//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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
 * @copydoc logging.hpp
 */

// HAVE_CONFIG_H, suzerain/config.h, and suzerain/common.hpp not present
// as this code should be useful in a broader, non-Suzerain context.

#include "logging.hpp"

#include <cassert>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// GCC 4.6.3 dislikes some Log4cxx constructs.  Mute -fpermissive to workaround.
// Warning-only idea taken from http://stackoverflow.com/questions/10932479
// The right thing to do is to fix log4cxx, but it seems to be a failed state.
//
// If relevant, push "GCC diagnostic warning -fpermissive"
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-fpermissive"
#endif

#include <log4cxx/basicconfigurator.h>
#include <log4cxx/file.h>
#include <log4cxx/helpers/bytearrayinputstream.h>
#include <log4cxx/helpers/bytebuffer.h>
#include <log4cxx/helpers/exception.h>
#include <log4cxx/helpers/fileinputstream.h>
#include <log4cxx/helpers/loglog.h>
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

// Matching pop for "GCC diagnostic warning -fpermissive"
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402
#pragma GCC diagnostic pop
#endif

using namespace log4cxx;
using namespace log4cxx::helpers;

/**
 * The Subversive Appender-Sifting HierarchyEventListener (SubversiveASHEL
 * (TM)) "blocks" adding Appenders outside a \c target portion of the
 * Logger Hierarchy named at construction time.  See the addAppenderEvent()
 * implementation for the necessary evil.
 */
class SubversiveASHEL : public virtual spi::HierarchyEventListener,
                        public virtual ObjectImpl
{
public:
    DECLARE_LOG4CXX_OBJECT(SubversiveASHEL)
    BEGIN_LOG4CXX_CAST_MAP()
            LOG4CXX_CAST_ENTRY(SubversiveASHEL)
            LOG4CXX_CAST_ENTRY(HierarchyEventListener)
    END_LOG4CXX_CAST_MAP()

    virtual void removeAppenderEvent(
        const LoggerPtr&, const AppenderPtr&) { /* NOP */ }

    virtual void addAppenderEvent(
        const LoggerPtr &logger, const AppenderPtr &appender);

    virtual void setTarget(const LogString& target);

protected:

    LogString target;
};

IMPLEMENT_LOG4CXX_OBJECT(SubversiveASHEL)

void SubversiveASHEL::setTarget(const LogString& target)
{
    if (!target.length())
        throw IllegalArgumentException("target.length() == 0");

    this->target = target;
}

void SubversiveASHEL::addAppenderEvent(const LoggerPtr &logger,
                                       const AppenderPtr &appender)
{
    if (!this->target.length())
        throw RuntimeException("SubversiveASHEL::setTarget() not called");

    // Do nothing if logger is at or below this->target.
    LoggerPtr ancestor = logger;
    do {
        if (0 == this->target.compare(ancestor->getName())) return;
    } while (ancestor = ancestor->getParent());

    // Suppress adding of appenders outside of this->target's subloggers
    const_cast<LoggerPtr&>(logger)->removeAppender(appender);

    // Add appenders intended for RootLogger to this->target
    if (0 == logger->getParent()) {
        Logger::getLogger(target)->addAppender(appender);
    }
}

/** http://old.nabble.com/Static-destruction-fiasco--td31026705.html */
static struct Log4cxxWorkaroundsType {
    Log4cxxWorkaroundsType() {
        APRInitializer::initialize();
    }
} workarounds;

namespace suzerain {

namespace support {

namespace logging {

logger_type rankzero;

logger_type allranks;

interceptor::interceptor(::std::ostream& stream)
    : s(stream), b(s.rdbuf(this))
{
}

interceptor::~interceptor()
{
    s.rdbuf(b);
}

::std::streamsize
interceptor::xsputn(const char *m, ::std::streamsize c)
{
    allranks->forcedLog(::log4cxx::Level::getAll(), ::std::string(m, c));
    return c;
}

void initialize(MPI_Comm, const char * const default_conf)
{
    static const char log4j_config_envvar[]   = "log4j.configuration";
    static const char log4cxx_config_envvar[] = "LOG4CXX_CONFIGURATION";

    // Programmatically enforce execution before log4cxx auto-configuration.
    if (LogManager::getLoggerRepository()->isConfigured()) {
        throw std::logic_error(
            "Apache log4cxx subsystem initialized before logger::initialize()."
            "\nThis is a CODING ERROR and will drastically hurt scalability.");
    }

    // Ensure MPI has been initialized
    {
        int flag;
        MPI_Initialized(&flag);
        assert(flag);
        if (!flag) throw std::logic_error("MPI_Init has not been invoked");
    }

    // Obtain our rank within MPI_COMM_WORLD
    int worldrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);

    // Build a rank-specific name for the "allranks" LoggerPtr instance
    LogString worldrankname;
    {
        std::ostringstream oss;
        oss << 'r' << worldrank;
        worldrankname = oss.str();
    }

    // For ranks > 0 create a HierarchyEventListener managed by the repository
    if (worldrank > 0) {
        SubversiveASHEL *p = new SubversiveASHEL();
        p->setTarget(worldrankname);
        LogManager::getLoggerRepository()->addHierarchyEventListener(p);
    }

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
        // DOMConfigurator::doConfigure(...) functionality and hitting disk
        // in any fashion slows down our startup for large rank counts.
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
        if (default_conf) {
            buf.assign(default_conf, default_conf + strlen(default_conf));
        }
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
            // Not the least bit ASCII 128+ safe but throws if issue arises
            // (Boost-less std::transform + boost::numeric::converter).
            std::vector<char>::const_iterator first = tmp.begin();
            std::vector<unsigned char>::iterator result = buf.begin();
            while (first != tmp.end()) {
                if (*first < 0) throw std::runtime_error(
                    "Logging configuration file must be ASCII 128");
                *result++ = static_cast<unsigned char>(*first++);
            }
#endif
        }

        // Broadcast configuration buffer length
        int buflen = buf.size();
        MPI_Bcast(&buflen, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (buflen) {

            // Configure logging subsystem using the nontrivial buffer contents
            InputStreamPtr bais(new ByteArrayInputStream(buf));
            Properties props;
            PropertyConfigurator pc;
            try {
                props.load(bais);
                pc.doConfigure(props, LogManager::getLoggerRepository());
            } catch (Exception &e) {
                BasicConfigurator::configure(); // Fallback
                buf.push_back('\0');
                LOG4CXX_WARN(Logger::getRootLogger(),
                    "Logging configuration from "
                    << (configpath.size() ? configpath : "default_conf")
                    << "invalid: " << e.what() << "\n" <<  &buf.front());
                buf.pop_back();
            }

            // Broadcast configuration buffer now that rank zero has finished
            MPI_Bcast(&buf.front(), buflen, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        } else {
            // No buffer contents so use a vanilla log4cxx configuration
            BasicConfigurator::configure();
        }

    } else {

        // Higher ranks receive buffer length
        int buflen;
        MPI_Bcast(&buflen, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (buflen) {

            // Configure logging subsystem using incoming buffer contents
            std::vector<unsigned char> buf(buflen);
            MPI_Bcast(&buf.front(), buflen, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
            InputStreamPtr bais(new ByteArrayInputStream(buf));
            Properties props;
            PropertyConfigurator pc;
            try {
                props.load(bais);
                pc.doConfigure(props, LogManager::getLoggerRepository());
            } catch (Exception & /* ignored */) {
                BasicConfigurator::configure(); // Fallback
            }

        } else {
            // No buffer contents so use a vanilla log4cxx configuration
            BasicConfigurator::configure();
        }

    }

    // Assert that all ranks indeed have initialized logging
    assert(LogManager::getLoggerRepository()->isConfigured());

    // "rankzero" is a disguise for the log4cxx::RootLogger where
    // SubversiveASHEL::addAppenderEvent() ensures nothing is emitted
    // on non-root ranks while still leaving isLevelEnabledFor() intact.
    rankzero = Logger::getRootLogger();

    // "allranks" emits events on all ranks by having SubversiveASHEL co-opt
    // appenders from RootLogger on non-root ranks.
    allranks = Logger::getLogger(worldrankname);

    // On non-zero ranks the RootLogger may have no appenders (which is fine).
    // However, it causes a worrisome one-time message like the following
    //     log4cxx: No appender could be found for logger (root).
    //     log4cxx: Please initialize the log4cxx system properly.
    // to appear once from each rank.  Temporarily disable log4cxx's LogLog,
    // emit a message into the ether, and turn LogLog back on again.
    if (worldrank > 0 && rankzero->getAllAppenders().empty()) {
        LogLog::setQuietMode(true);
        rankzero->fatal("Sorry, baby, but I had to crash that Honda.");
        LogLog::setQuietMode(false);
    }
}

} // end namespace logger

} // end namespace support

} // end namespace suzerain
