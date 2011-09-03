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
#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <mpi.h>
#include <log4cxx/logger.h>

/**
 * Namespace to house MPI-aware logging abstraction built atop log4cxx.
 */
namespace logging {

/** Smart pointer to a particular logger instance */
typedef ::log4cxx::LoggerPtr logger_type;

extern logger_type rankzero;

/** Logger enabled on all ranks */
extern logger_type all;

/**
 * Initialize the logging infrastructure.
 * Must be called after \c MPI_Init and before any logging statements.
 *
 * @param comm Must be MPI_COMM_WORLD.  Present as a placeholder
 *             to communicate initialization order requirements
 *             relative to \c MPI_Init.
 */
void initialize(MPI_Comm comm);

} // end namespace logging

// Logging macros that log from all ranks using "r12345"
#define TRACE(m)       do{LOG4CXX_TRACE(::logging::all,m)}while(0)
#define DEBUG(m)       do{LOG4CXX_DEBUG(::logging::all,m)}while(0)
#define INFO(m)        do{LOG4CXX_INFO( ::logging::all,m)}while(0)
#define WARN(m)        do{LOG4CXX_WARN( ::logging::all,m)}while(0)
#define ERROR(m)       do{LOG4CXX_ERROR(::logging::all,m)}while(0)
#define FATAL(m)       do{LOG4CXX_FATAL(::logging::all,m)}while(0)

// Macros useful for checking if particular levels are enabled
#define INFO_ENABLED   (::logging::rankzero->isInfoEnabled())
#define TRACE_ENABLED  (::logging::rankzero->isTraceEnabled())
#define DEBUG_ENABLED  (::logging::rankzero->isDebugEnabled())

// Logging macros that log only from MPI rank 0 using "root"
#define TRACE0(m)     do{LOG4CXX_TRACE(::logging::rankzero,m)}while(0)
#define DEBUG0(m)     do{LOG4CXX_DEBUG(::logging::rankzero,m)}while(0)
#define INFO0(m)      do{LOG4CXX_INFO( ::logging::rankzero,m)}while(0)
#define WARN0(m)      do{LOG4CXX_WARN( ::logging::rankzero,m)}while(0)
#define ERROR0(m)     do{LOG4CXX_ERROR(::logging::rankzero,m)}while(0)
#define FATAL0(m)     do{LOG4CXX_FATAL(::logging::rankzero,m)}while(0)

// Logging macros taking one-off logger names for infrequent use
// Provided so that configuration can pull off particular outputs, e.g. L2
#define TRACEDUB(n,m) do{LOG4CXX_TRACE(::log4cxx::Logger::getLogger(n),m)}while(0)
#define DEBUGDUB(n,m) do{LOG4CXX_DEBUG(::log4cxx::Logger::getLogger(n),m)}while(0)
#define INFODUB(n,m)  do{LOG4CXX_INFO( ::log4cxx::Logger::getLogger(n),m)}while(0)
#define WARNDUB(n,m)  do{LOG4CXX_WARN( ::log4cxx::Logger::getLogger(n),m)}while(0)
#define ERRORDUB(n,m) do{LOG4CXX_ERROR(::log4cxx::Logger::getLogger(n),m)}while(0)
#define FATALDUB(n,m) do{LOG4CXX_FATAL(::log4cxx::Logger::getLogger(n),m)}while(0)

// Macros for checking if particular levels are enabled for named loggers
#define INFODUB_ENABLED(n)  (::log4cxx::Logger::getLogger(n)->isInfoEnabled())
#define TRACEDUB_ENABLED(n) (::log4cxx::Logger::getLogger(n)->isTraceEnabled())
#define DEBUGDUB_ENABLED(n) (::log4cxx::Logger::getLogger(n)->isDebugEnabled())


#endif // LOGGER_HPP
