//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_SUPPORT_LOGGING_HPP
#define SUZERAIN_SUPPORT_LOGGING_HPP

/** @file
 * MPI-aware logging built atop log4cxx
 */

#include <streambuf>
#include <ostream>
#include <string>
#include <mpi.h>
#include <log4cxx/logger.h>

/** @name Logging macros emitting messages only from rank zero.
 *  @{
 */

/** \def TRACE0
 * Log a message at trace level on rank zero.
 * With one argument,  log the provided message using name "root".
 * With two arguments, the arguments are the name and message.
 */

/** \def DEBUG0
 * Log a message at debug level on rank zero.
 * With one argument,  log the provided message using name "root".
 * With two arguments, the arguments are the name and message.
 */

/** \def INFO0
 * Log a message at info level on rank zero.
 * With one argument,  log the provided message using name "root".
 * With two arguments, the arguments are the name and message.
 */

/** \def WARN0
 * Log a message at warn level on rank zero.
 * With one argument,  log the provided message using name "root".
 * With two arguments, the arguments are the name and message.
 */

/** \def ERROR0
 * Log a message at error level on rank zero.
 * With one argument,  log the provided message using name "root".
 * With two arguments, the arguments are the name and message.
 */

/** \def FATAL0
 * Log a message at fatal level on rank zero.
 * With one argument,  log the provided message using name "root".
 * With two arguments, the arguments are the name and message.
 */

/* @} */


/** @name Logging macros emitting messages from all ranks.
 *  @{
 */

/** \def TRACE
 * Log a message at trace level from all ranks with a rank-specific name.
 * With one argument,  log the provided message.
 * With two arguments, the arguments are the name suffix and message.
 */

/** \def DEBUG
 * Log a message at debug level from all ranks with a rank-specific name.
 * With one argument,  log the provided message.
 * With two arguments, the arguments are the name suffix and message.
 */

/** \def INFO
 * Log a message at info level from all ranks with a rank-specific name.
 * With one argument,  log the provided message.
 * With two arguments, the arguments are the name suffix and message.
 */

/** \def WARN
 * Log a message at warn level from all ranks with a rank-specific name.
 * With one argument,  log the provided message.
 * With two arguments, the arguments are the name suffix and message.
 */

/** \def ERROR
 * Log a message at error level from all ranks with a rank-specific name.
 * With one argument,  log the provided message.
 * With two arguments, the arguments are the name suffix and message.
 */

/** \def FATAL
 * Log a message at fatal level from all ranks with a rank-specific name.
 * With one argument,  log the provided message.
 * With two arguments, the arguments are the name suffix and message.
 */

/* @} */


/** @name Macros checking if particular log levels are enabled on rank zero.
 *  @{
 */

/** \def TRACE0_ENABLED
 * Is logging enabled at trace level on rank zero?
 * With no arguments, check enabled status using name "root".
 * With one argument, check enabled status using provided name.
 */

/** \def DEBUG0_ENABLED
 * Is logging enabled at debug level on rank zero?
 * With no arguments, check enabled status using name "root".
 * With one argument, check enabled status using provided name.
 */

/** \def INFO0_ENABLED
 * Is logging enabled at info level on rank zero?
 * With no arguments, check enabled status using name "root".
 * With one argument, check enabled status using provided name.
 */

/** \def WARN0_ENABLED
 * Is logging enabled at warn level on rank zero?
 * With no arguments, check enabled status using name "root".
 * With one argument, check enabled status using provided name.
 */

/** \def ERROR0_ENABLED
 * Is logging enabled at error level on rank zero?
 * With no arguments, check enabled status using name "root".
 * With one argument, check enabled status using provided name.
 */

/** \def FATAL0_ENABLED
 * Is logging enabled at fatal level on rank zero?
 * With no arguments, check enabled status using name "root".
 * With one argument, check enabled status using provided name.
 */

/* @} */


/** @name Macros checking if particular log levels are enabled on all ranks.
 *  @{
 */

/** \def TRACE_ENABLED
 * Is logging enabled at trace level on all ranks?
 * With no arguments, check enabled status using rank-specific name.
 * With one argument, check enabled status using provided name suffix.
 */

/** \def DEBUG_ENABLED
 * Is logging enabled at debug level on all ranks?
 * With no arguments, check enabled status using rank-specific name.
 * With one argument, check enabled status using provided name suffix.
 */

/** \def INFO_ENABLED
 * Is logging enabled at info level on all ranks?
 * With no arguments, check enabled status using rank-specific name.
 * With one argument, check enabled status using provided name suffix.
 */

/** \def WARN_ENABLED
 * Is logging enabled at warn level on all ranks?
 * With no arguments, check enabled status using rank-specific name.
 * With one argument, check enabled status using provided name suffix.
 */

/** \def ERROR_ENABLED
 * Is logging enabled at error level on all ranks?
 * With no arguments, check enabled status using rank-specific name.
 * With one argument, check enabled status using provided name suffix.
 */

/** \def FATAL_ENABLED
 * Is logging enabled at fatal level on all ranks?
 * With no arguments, check enabled status using rank-specific name.
 * With one argument, check enabled status using provided name suffix.
 */

/* @} */


/** @name Macros enabling logging at particular log levels on rank zero.
 *  @{
 */

/** \def TRACE0_ENABLE
 * Enable logging at trace level on rank zero.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def DEBUG0_ENABLE
 * Enable logging at debug level on rank zero.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def INFO0_ENABLE
 * Enable logging at info level on rank zero.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def WARN0_ENABLE
 * Enable logging at warn level on rank zero.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def ERROR0_ENABLE
 * Enable logging at error level on rank zero.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def FATAL0_ENABLE
 * Enable logging at fatal level on rank zero.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/* @} */


/** @name Macros enabling logging at particular log levels on all ranks.
 *  @{
 */

/** \def TRACE_ENABLE
 * Enable logging at trace level on all ranks.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def DEBUG_ENABLE
 * Enable logging at trace level on all ranks.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def INFO_ENABLE
 * Enable logging at trace level on all ranks.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def WARN_ENABLE
 * Enable logging at trace level on all ranks.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def ERROR_ENABLE
 * Enable logging at trace level on all ranks.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/** \def FATAL_ENABLE
 * Enable logging at trace level on all ranks.
 * With no arguments, enable for name "root".
 * With one argument, enable for provided name.
 */

/* @} */

namespace suzerain {

namespace support {

/**
 * Namespace to house MPI-aware logging abstraction built atop log4cxx.
 */
namespace logging {

/**
 * Initialize the logging infrastructure.
 * Must be called after \c MPI_Init and before any logging statements.
 *
 * @param comm Must be MPI_COMM_WORLD.  Present as a placeholder
 *             to clearly communicate initialization order requirements
 *             relative to \c MPI_Init.
 * @param default_conf (Optional) Default configuration to use when neither
 *             \c log4j.configuration nor \c LOG4CXX_CONFIGURATION
 *             is found in the environment.  If \c NULL or zero-length,
 *             log4cxx::BasicConfigurator is used.
 */
void initialize(MPI_Comm comm, const char *default_conf = NULL);

/** Smart pointer to some particular logger instance */
typedef ::log4cxx::LoggerPtr logger_type;

/**
 * Logger that emits messages only on rank zero.  Use the TRACE0, DEBUG0,
 * INFO0, WARN0, ERROR0, and FATAL0 macros to employ this logger.
 */
extern logger_type rankzero;

/**
 * Logger that emits messages on all ranks with a rank-dependent name.  Use the
 * TRACE, DEBUG, INFO, WARN, ERROR, and FATAL macros to employ this logger.
 */
extern logger_type allranks;

/**
 * RAII class that redirects all output from a stream through #allranks.
 * Implementation adapted from http://stackoverflow.com/questions/533038.
 * Implementation inlined for brevity, not because virtual inlining works.
 */
class interceptor : public ::std::streambuf {

public:
    interceptor(::std::ostream& stream);
    ~interceptor();

protected:
    virtual ::std::streamsize xsputn(const char *m, ::std::streamsize c);

private:
    ::std::ostream   &s;
    ::std::streambuf *b;

};

/** Retrieve a particular logger by name */
template <typename S>
inline logger_type get_logger(const S& name)
{
    return ::log4cxx::Logger::getLogger(name);
}

/**
 * Retrieve a particular logger by name with
 * possibility of returning default initialized logger_types.
 */
inline logger_type get_logger(const char * const name)
{
    return name ? ::log4cxx::Logger::getLogger(name) : logger_type();
}

/** Idempotent overload of get_logger() to aid macro definitions */
inline logger_type& get_logger(logger_type& logger)
{
    return logger;
}

/** Retrieve a descendent of a particular logger by name */
template <typename S, typename T>
inline logger_type get_logger(const S& name, const T& child)
{
    std::string s(name);
    s.append(1, '.');
    s.append(child);
    return ::log4cxx::Logger::getLogger(s);
}

/** Retrieve a descendent of a particular logger by name */
template <typename T>
inline logger_type get_logger(logger_type& logger, const T& child)
{
    std::string s;
    logger->getName(s);
    s.append(1, '.');
    s.append(child);
    return ::log4cxx::Logger::getLogger(s);
}

/** Idempotent overload of get_logger() to aid macro definitions */
inline logger_type& get_logger(logger_type&, logger_type& child)
{
    return child;
}

} // end namespace logging

} // end namespace support

} // end namespace suzerain

// Implementation details supporting variadic macro dispatching
#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Adopted from comp.std.c posting by Laurent Deniau on 16 Jan 2006
#define LOGGING_PP_NARG(...)  LOGGING_PP_NARG_(__VA_ARGS__,LOGGING_PP_RSEQ_N())
#define LOGGING_PP_NARG_(...) LOGGING_PP_ARG_N(__VA_ARGS__)
#define LOGGING_PP_ARG_N(_1, _2, _3, N, ...) N
#define LOGGING_PP_RSEQ_N() 3, 2, 1, 0

// Allows invoking macros like prefix_ARITY# for some number #
#define LOGGING_CONCAT(a,b,c) a##b##c
#define LOGGING_DISPATCH(prefix,count,...) \
        LOGGING_CONCAT(prefix,_ARITY,count)(__VA_ARGS__)

// Dispatch to macro level0_ARITY# for # in {1, 2}
#define TRACE0(...) \
        LOGGING_DISPATCH(TRACE0,LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define DEBUG0(...) \
        LOGGING_DISPATCH(DEBUG0,LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define INFO0(...) \
        LOGGING_DISPATCH(INFO0, LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define WARN0(...) \
        LOGGING_DISPATCH(WARN0, LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define ERROR0(...) \
        LOGGING_DISPATCH(ERROR0,LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define FATAL0(...) \
        LOGGING_DISPATCH(FATAL0,LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)

// Dispatch to macro level_ARITY# for # in {1, 2}
#define TRACE(...) \
        LOGGING_DISPATCH(TRACE,LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define DEBUG(...) \
        LOGGING_DISPATCH(DEBUG,LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define INFO(...) \
        LOGGING_DISPATCH(INFO, LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define WARN(...) \
        LOGGING_DISPATCH(WARN, LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define ERROR(...) \
        LOGGING_DISPATCH(ERROR,LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)
#define FATAL(...) \
        LOGGING_DISPATCH(FATAL,LOGGING_PP_NARG(__VA_ARGS__),__VA_ARGS__)

// Dispatch to macro level_ENABLED_ARITY# for # in {1,2}
// Auxiliary dummy argument added since arity zero breaks LOGGING_PP_NARG
#define TRACE0_ENABLED(...) \
        LOGGING_DISPATCH(TRACE0_ENABLED,LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define DEBUG0_ENABLED(...) \
        LOGGING_DISPATCH(DEBUG0_ENABLED,LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define INFO0_ENABLED(...) \
        LOGGING_DISPATCH(INFO0_ENABLED, LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define WARN0_ENABLED(...) \
        LOGGING_DISPATCH(WARN0_ENABLED, LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define ERROR0_ENABLED(...) \
        LOGGING_DISPATCH(ERROR0_ENABLED,LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define FATAL0_ENABLED(...) \
        LOGGING_DISPATCH(FATAL0_ENABLED,LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)

// Dispatch to macro level_ENABLED_ARITY# for # in {1,2}
// Auxiliary dummy argument added since arity zero breaks LOGGING_PP_NARG
#define TRACE_ENABLED(...) \
        LOGGING_DISPATCH(TRACE_ENABLED,LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define DEBUG_ENABLED(...) \
        LOGGING_DISPATCH(DEBUG_ENABLED,LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define INFO_ENABLED(...) \
        LOGGING_DISPATCH(INFO_ENABLED, LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define WARN_ENABLED(...) \
        LOGGING_DISPATCH(WARN_ENABLED, LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define ERROR_ENABLED(...) \
        LOGGING_DISPATCH(ERROR_ENABLED,LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define FATAL_ENABLED(...) \
        LOGGING_DISPATCH(FATAL_ENABLED,LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)

// Dispatch to macro level_ENABLE_ARITY# for # in {1,2}
// Auxiliary dummy argument added since arity zero breaks LOGGING_PP_NARG
#define TRACE0_ENABLE(...) \
        LOGGING_DISPATCH(TRACE0_ENABLE,LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define DEBUG0_ENABLE(...) \
        LOGGING_DISPATCH(DEBUG0_ENABLE,LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define INFO0_ENABLE(...) \
        LOGGING_DISPATCH(INFO0_ENABLE, LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define WARN0_ENABLE(...) \
        LOGGING_DISPATCH(WARN0_ENABLE, LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define ERROR0_ENABLE(...) \
        LOGGING_DISPATCH(ERROR0_ENABLE,LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)
#define FATAL0_ENABLE(...) \
        LOGGING_DISPATCH(FATAL0_ENABLE,LOGGING_PP_NARG("dummy",##__VA_ARGS__),__VA_ARGS__)

// Dispatch to macro level_ENABLE_ARITY# for # in {1,2}
// Auxiliary dummy argument added since arity zero breaks LOGGING_PP_NARG
#define TRACE_ENABLE(...) \
        LOGGING_DISPATCH(TRACE_ENABLE,LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define DEBUG_ENABLE(...) \
        LOGGING_DISPATCH(DEBUG_ENABLE,LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define INFO_ENABLE(...) \
        LOGGING_DISPATCH(INFO_ENABLE, LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define WARN_ENABLE(...) \
        LOGGING_DISPATCH(WARN_ENABLE, LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define ERROR_ENABLE(...) \
        LOGGING_DISPATCH(ERROR_ENABLE,LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)
#define FATAL_ENABLE(...) \
        LOGGING_DISPATCH(FATAL_ENABLE,LOGGING_PP_NARG("dummy",##__VA_ARGS__),##__VA_ARGS__)

// Logging macros emitting messages as "root" only from rank zero.
#define TRACE0_ARITY1(m) do{LOG4CXX_TRACE(::suzerain::support::logging::rankzero,m)}while(0)
#define DEBUG0_ARITY1(m) do{LOG4CXX_DEBUG(::suzerain::support::logging::rankzero,m)}while(0)
#define INFO0_ARITY1(m)  do{LOG4CXX_INFO( ::suzerain::support::logging::rankzero,m)}while(0)
#define WARN0_ARITY1(m)  do{LOG4CXX_WARN( ::suzerain::support::logging::rankzero,m)}while(0)
#define ERROR0_ARITY1(m) do{LOG4CXX_ERROR(::suzerain::support::logging::rankzero,m)}while(0)
#define FATAL0_ARITY1(m) do{LOG4CXX_FATAL(::suzerain::support::logging::rankzero,m)}while(0)

// Logging macros emitting messages from named loggers on rank zero.
#define TRACE0_ARITY2(l,m) do{LOG4CXX_TRACE(::suzerain::support::logging::get_logger(l),m)}while(0)
#define DEBUG0_ARITY2(l,m) do{LOG4CXX_DEBUG(::suzerain::support::logging::get_logger(l),m)}while(0)
#define INFO0_ARITY2(l,m)  do{LOG4CXX_INFO( ::suzerain::support::logging::get_logger(l),m)}while(0)
#define WARN0_ARITY2(l,m)  do{LOG4CXX_WARN( ::suzerain::support::logging::get_logger(l),m)}while(0)
#define ERROR0_ARITY2(l,m) do{LOG4CXX_ERROR(::suzerain::support::logging::get_logger(l),m)}while(0)
#define FATAL0_ARITY2(l,m) do{LOG4CXX_FATAL(::suzerain::support::logging::get_logger(l),m)}while(0)

// Logging macros to check if particular logging levels are enabled.
// Arity is off-by-one to handle limitation of LOGGING_PP_NARG
#define TRACE0_ENABLED_ARITY1() (::suzerain::support::logging::rankzero->isTraceEnabled())
#define DEBUG0_ENABLED_ARITY1() (::suzerain::support::logging::rankzero->isDebugEnabled())
#define INFO0_ENABLED_ARITY1()  (::suzerain::support::logging::rankzero->isInfoEnabled())
#define WARN0_ENABLED_ARITY1()  (::suzerain::support::logging::rankzero->isWarnEnabled())
#define ERROR0_ENABLED_ARITY1() (::suzerain::support::logging::rankzero->isErrorEnabled())
#define FATAL0_ENABLED_ARITY1() (::suzerain::support::logging::rankzero->isFatalEnabled())

// Logging macros to check if named logging levels are enabled.
// Arity is off-by-one to handle limitation of LOGGING_PP_NARG
#define TRACE0_ENABLED_ARITY2(l) (::suzerain::support::logging::get_logger(l)->isTraceEnabled())
#define DEBUG0_ENABLED_ARITY2(l) (::suzerain::support::logging::get_logger(l)->isDebugEnabled())
#define INFO0_ENABLED_ARITY2(l)  (::suzerain::support::logging::get_logger(l)->isInfoEnabled())
#define WARN0_ENABLED_ARITY2(l)  (::suzerain::support::logging::get_logger(l)->isWarnEnabled())
#define ERROR0_ENABLED_ARITY2(l) (::suzerain::support::logging::get_logger(l)->isErrorEnabled())
#define FATAL0_ENABLED_ARITY2(l) (::suzerain::support::logging::get_logger(l)->isFatalEnabled())

// Logging macros to enable logging at and above at particular log levels.
// Arity is off-by-one to handle limitation of LOGGING_PP_NARG
#define TRACE0_ENABLE_ARITY1() do{::suzerain::support::logging::rankzero->setLevel(::log4cxx::Level::getTrace());}while(0)
#define DEBUG0_ENABLE_ARITY1() do{::suzerain::support::logging::rankzero->setLevel(::log4cxx::Level::getDebug());}while(0)
#define INFO0_ENABLE_ARITY1()  do{::suzerain::support::logging::rankzero->setLevel(::log4cxx::Level::getInfo());}while(0)
#define WARN0_ENABLE_ARITY1()  do{::suzerain::support::logging::rankzero->setLevel(::log4cxx::Level::getWarn());}while(0)
#define ERROR0_ENABLE_ARITY1() do{::suzerain::support::logging::rankzero->setLevel(::log4cxx::Level::getError());}while(0)
#define FATAL0_ENABLE_ARITY1() do{::suzerain::support::logging::rankzero->setLevel(::log4cxx::Level::getFatal());}while(0)

// Logging macros to enable logging at and above particular log levels.
// Arity is off-by-one to handle limitation of LOGGING_PP_NARG
#define TRACE0_ENABLE_ARITY2(l) do{::suzerain::support::logging::get_logger(l)->setLevel(::log4cxx::Level::getTrace());}while(0)
#define DEBUG0_ENABLE_ARITY2(l) do{::suzerain::support::logging::get_logger(l)->setLevel(::log4cxx::Level::getDebug());}while(0)
#define INFO0_ENABLE_ARITY2(l)  do{::suzerain::support::logging::get_logger(l)->setLevel(::log4cxx::Level::getInfo());}while(0)
#define WARN0_ENABLE_ARITY2(l)  do{::suzerain::support::logging::get_logger(l)->setLevel(::log4cxx::Level::getWarn());}while(0)
#define ERROR0_ENABLE_ARITY2(l) do{::suzerain::support::logging::get_logger(l)->setLevel(::log4cxx::Level::getError());}while(0)
#define FATAL0_ENABLE_ARITY2(l) do{::suzerain::support::logging::get_logger(l)->setLevel(::log4cxx::Level::getFatal());}while(0)

// Logging macros to enable logging at and above at particular log levels.
// Arity is off-by-one to handle limitation of LOGGING_PP_NARG
#define TRACE_ENABLE_ARITY1() do{::suzerain::support::logging::allranks->setLevel(::log4cxx::Level::getTrace());}while(0)
#define DEBUG_ENABLE_ARITY1() do{::suzerain::support::logging::allranks->setLevel(::log4cxx::Level::getDebug());}while(0)
#define INFO_ENABLE_ARITY1()  do{::suzerain::support::logging::allranks->setLevel(::log4cxx::Level::getInfo());}while(0)
#define WARN_ENABLE_ARITY1()  do{::suzerain::support::logging::allranks->setLevel(::log4cxx::Level::getWarn());}while(0)
#define ERROR_ENABLE_ARITY1() do{::suzerain::support::logging::allranks->setLevel(::log4cxx::Level::getError());}while(0)
#define FATAL_ENABLE_ARITY1() do{::suzerain::support::logging::allranks->setLevel(::log4cxx::Level::getFatal());}while(0)

// Logging macros to enable logging at and above particular log levels.
// Arity is off-by-one to handle limitation of LOGGING_PP_NARG
#define TRACE_ENABLE_ARITY2(l) do{::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->setLevel(::log4cxx::Level::getTrace());}while(0)
#define DEBUG_ENABLE_ARITY2(l) do{::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->setLevel(::log4cxx::Level::getDebug());}while(0)
#define INFO_ENABLE_ARITY2(l)  do{::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->setLevel(::log4cxx::Level::getInfo());}while(0)
#define WARN_ENABLE_ARITY2(l)  do{::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->setLevel(::log4cxx::Level::getWarn());}while(0)
#define ERROR_ENABLE_ARITY2(l) do{::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->setLevel(::log4cxx::Level::getError());}while(0)
#define FATAL_ENABLE_ARITY2(l) do{::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->setLevel(::log4cxx::Level::getFatal());}while(0)

// Logging macros emitting messages from all ranks with rank-specific names.
#define TRACE_ARITY1(m) do{LOG4CXX_TRACE(::suzerain::support::logging::allranks,m)}while(0)
#define DEBUG_ARITY1(m) do{LOG4CXX_DEBUG(::suzerain::support::logging::allranks,m)}while(0)
#define INFO_ARITY1(m)  do{LOG4CXX_INFO( ::suzerain::support::logging::allranks,m)}while(0)
#define WARN_ARITY1(m)  do{LOG4CXX_WARN( ::suzerain::support::logging::allranks,m)}while(0)
#define ERROR_ARITY1(m) do{LOG4CXX_ERROR(::suzerain::support::logging::allranks,m)}while(0)
#define FATAL_ARITY1(m) do{LOG4CXX_FATAL(::suzerain::support::logging::allranks,m)}while(0)

// Logging macros emitting messages from named loggers on all ranks
#define TRACE_ARITY2(l,m) do{LOG4CXX_TRACE(::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l),m)}while(0)
#define DEBUG_ARITY2(l,m) do{LOG4CXX_DEBUG(::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l),m)}while(0)
#define INFO_ARITY2(l,m)  do{LOG4CXX_INFO( ::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l),m)}while(0)
#define WARN_ARITY2(l,m)  do{LOG4CXX_WARN( ::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l),m)}while(0)
#define ERROR_ARITY2(l,m) do{LOG4CXX_ERROR(::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l),m)}while(0)
#define FATAL_ARITY2(l,m) do{LOG4CXX_FATAL(::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l),m)}while(0)

// Logging macros to check if particular logging levels are enabled.
// Arity is off-by-one to handle limitation of LOGGING_PP_NARG
#define TRACE_ENABLED_ARITY1() (::suzerain::support::logging::allranks->isTraceEnabled())
#define DEBUG_ENABLED_ARITY1() (::suzerain::support::logging::allranks->isDebugEnabled())
#define INFO_ENABLED_ARITY1()  (::suzerain::support::logging::allranks->isInfoEnabled())
#define WARN_ENABLED_ARITY1()  (::suzerain::support::logging::allranks->isWarnEnabled())
#define ERROR_ENABLED_ARITY1() (::suzerain::support::logging::allranks->isErrorEnabled())
#define FATAL_ENABLED_ARITY1() (::suzerain::support::logging::allranks->isFatalEnabled())

// Logging macros to check if named logging levels are enabled.
// Arity is off-by-one to handle limitation of LOGGING_PP_NARG
#define TRACE_ENABLED_ARITY2(l) (::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->isTraceEnabled())
#define DEBUG_ENABLED_ARITY2(l) (::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->isDebugEnabled())
#define INFO_ENABLED_ARITY2(l)  (::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->isInfoEnabled())
#define WARN_ENABLED_ARITY2(l)  (::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->isWarnEnabled())
#define ERROR_ENABLED_ARITY2(l) (::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->isErrorEnabled())
#define FATAL_ENABLED_ARITY2(l) (::suzerain::support::logging::get_logger(::suzerain::support::logging::allranks,l)->isFatalEnabled())

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // SUZERAIN_SUPPORT_LOGGING_HPP
