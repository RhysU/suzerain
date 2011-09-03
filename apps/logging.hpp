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
 * logging.hpp: MPI-aware logging tools built atop log4cxx
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <mpi.h>
#include <log4cxx/logger.h>

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
 */
void initialize(MPI_Comm comm);

/** Smart pointer to some particular logger instance */
typedef ::log4cxx::LoggerPtr logger_type;

/** Retrieve a particular logger by name */
inline logger_type get_logger(const std::string& name)
{
    return ::log4cxx::Logger::getLogger(name);
}

/** Retrieve a particular logger by name */
inline logger_type get_logger(const char *name)
{
    return ::log4cxx::Logger::getLogger(name);
}

/** Idempotent overload of get_logger() to aid macro definitions */
inline logger_type& get_logger(logger_type& logger) { return logger; }

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

} // end namespace logging


/** @name Logging macros emitting messages as "root" only from rank zero.
 *  @{
 */

/** Log message \c m at trace level on rank zero. */
#define TRACE0(m) do{LOG4CXX_TRACE(::logging::rankzero,m)}while(0)
/** Log message \c m at debug level on rank zero. */
#define DEBUG0(m) do{LOG4CXX_DEBUG(::logging::rankzero,m)}while(0)
/** Log message \c m at info level on rank zero. */
#define INFO0(m)  do{LOG4CXX_INFO( ::logging::rankzero,m)}while(0)
/** Log message \c m at warn level on rank zero. */
#define WARN0(m)  do{LOG4CXX_WARN( ::logging::rankzero,m)}while(0)
/** Log message \c m at error level on rank zero. */
#define ERROR0(m) do{LOG4CXX_ERROR(::logging::rankzero,m)}while(0)
/** Log message \c m at fatal level on rank zero. */
#define FATAL0(m) do{LOG4CXX_FATAL(::logging::rankzero,m)}while(0)

/* @} */


/** @name Logging macros emitting messages from all ranks.
 *  @{
 */

/** Log message \c m at trace level from all ranks. */
#define TRACE(m)  do{LOG4CXX_TRACE(::logging::allranks,m)}while(0)
/** Log message \c m at debug level from all ranks. */
#define DEBUG(m)  do{LOG4CXX_DEBUG(::logging::allranks,m)}while(0)
/** Log message \c m at info level from all ranks. */
#define INFO(m)   do{LOG4CXX_INFO( ::logging::allranks,m)}while(0)
/** Log message \c m at warn level from all ranks. */
#define WARN(m)   do{LOG4CXX_WARN( ::logging::allranks,m)}while(0)
/** Log message \c m at error level from all ranks. */
#define ERROR(m)  do{LOG4CXX_ERROR(::logging::allranks,m)}while(0)
/** Log message \c m at fatal level from all ranks. */
#define FATAL(m)  do{LOG4CXX_FATAL(::logging::allranks,m)}while(0)

/* @} */


/** @name Macros to check if particular logging levels are enabled.
 *  @{
 */

/** Return true if trace level logging is enabled. */
#define TRACE_ENABLED  (::logging::rankzero->isTraceEnabled())
/** Return true if debug level logging is enabled. */
#define DEBUG_ENABLED  (::logging::rankzero->isDebugEnabled())
/** Return true if info level logging is enabled. */
#define INFO_ENABLED   (::logging::rankzero->isInfoEnabled())
/** Return true if warn level logging is enabled. */
#define WARN_ENABLED   (::logging::rankzero->isWarnEnabled())
/** Return true if error level logging is enabled. */
#define ERROR_ENABLED  (::logging::rankzero->isErrorEnabled())
/** Return true if fatal level logging is enabled. */
#define FATAL_ENABLED  (::logging::rankzero->isFatalEnabled())

/* @} */

/** @name Logging macros emitting messages from named loggers on rank zero.
 *  @{
 */

/** Log message \c m at trace level on rank zero using logger \c l. */
#define LTRACE(l,m) do{LOG4CXX_TRACE(::logging::get_logger(l),m)}while(0)
/** Log message \c m at debug level on rank zero using logger \c l. */
#define LDEBUG(l,m) do{LOG4CXX_DEBUG(::logging::get_logger(l),m)}while(0)
/** Log message \c m at info level on rank zero using logger \c l. */
#define LINFO(l,m)  do{LOG4CXX_INFO( ::logging::get_logger(l),m)}while(0)
/** Log message \c m at warn level on rank zero using logger \c l. */
#define LWARN(l,m)  do{LOG4CXX_WARN( ::logging::get_logger(l),m)}while(0)
/** Log message \c m at error level on rank zero using logger \c l. */
#define LERROR(l,m) do{LOG4CXX_ERROR(::logging::get_logger(l),m)}while(0)
/** Log message \c m at fatal level on rank zero using logger \c l. */
#define LFATAL(l,m) do{LOG4CXX_FATAL(::logging::get_logger(l),m)}while(0)

/* @} */


/** @name Macros to check if particular named logging levels are enabled.
 *  @{
 */

/** Return true if trace level logging is enabled for logger \c l. */
#define LTRACE_ENABLED(l) (::logging::get_logger(l)->isTraceEnabled())
/** Return true if debug level logging is enabled for logger \c l. */
#define LDEBUG_ENABLED(l) (::logging::get_logger(l)->isDebugEnabled())
/** Return true if info level logging is enabled for logger \c l. */
#define LINFO_ENABLED(l)  (::logging::get_logger(l)->isInfoEnabled())
/** Return true if warn level logging is enabled for logger \c l. */
#define LWARN_ENABLED(l)  (::logging::get_logger(l)->isWarnEnabled())
/** Return true if error level logging is enabled for logger \c l. */
#define LERROR_ENABLED(l) (::logging::get_logger(l)->isErrorEnabled())
/** Return true if fatal level logging is enabled for logger \c l. */
#define LFATAL_ENABLED(l) (::logging::get_logger(l)->isFatalEnabled())

/* @} */

#endif // LOGGING_HPP
