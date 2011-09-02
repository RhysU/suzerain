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
 * logger.hpp: logging tools build atop log4cxx
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <log4cxx/logger.h>

/**
 * Namespace to house logging abstraction details.
 */
namespace logger {

/** Logger enabled only on rank zero */
extern ::log4cxx::LoggerPtr rankzero;

/** Logger enabled on all ranks */
extern ::log4cxx::LoggerPtr all;

/**
 * Initialize the logging infrastructure.
 * Must be called after \c MPI_Init and before any logging statements.
 *
 * @param comm Must be MPI_COMM_WORLD.  Present as a placeholder
 *             to communicate initialization order requirements
 *             relative to \c MPI_Init.
 */
void initialize(MPI_Comm comm);

} // end namespace logger

// Logging macros that log only from MPI rank 0
#define TRACE0(msg)     do{LOG4CXX_TRACE(::logger::rankzero,msg)}while(0)
#define DEBUG0(msg)     do{LOG4CXX_DEBUG(::logger::rankzero,msg)}while(0)
#define INFO0(msg)      do{LOG4CXX_INFO( ::logger::rankzero,msg)}while(0)
#define WARN0(msg)      do{LOG4CXX_WARN( ::logger::rankzero,msg)}while(0)
#define ERROR0(msg)     do{LOG4CXX_ERROR(::logger::rankzero,msg)}while(0)
#define FATAL0(msg)     do{LOG4CXX_FATAL(::logger::rankzero,msg)}while(0)

// Logging macros that log from all ranks
#define TRACE(msg)     do{LOG4CXX_TRACE(::logger::all,msg)}while(0)
#define DEBUG(msg)     do{LOG4CXX_DEBUG(::logger::all,msg)}while(0)
#define INFO(msg)      do{LOG4CXX_INFO( ::logger::all,msg)}while(0)
#define WARN(msg)      do{LOG4CXX_WARN( ::logger::all,msg)}while(0)
#define ERROR(msg)     do{LOG4CXX_ERROR(::logger::all,msg)}while(0)
#define FATAL(msg)     do{LOG4CXX_FATAL(::logger::all,msg)}while(0)

// Logging macros taking one-off logger names for infrequent use
// Provided so that configuration can pull off particular outputs, e.g. L2
#define TRACEDUB(nick,msg)     do{LOG4CXX_TRACE(::log4cxx::Logger::getLogger(nick),msg)}while(0)
#define DEBUGDUB(nick,msg)     do{LOG4CXX_DEBUG(::log4cxx::Logger::getLogger(nick),msg)}while(0)
#define INFODUB(nick,msg)      do{LOG4CXX_INFO( ::log4cxx::Logger::getLogger(nick),msg)}while(0)
#define WARNDUB(nick,msg)      do{LOG4CXX_WARN( ::log4cxx::Logger::getLogger(nick),msg)}while(0)
#define ERRORDUB(nick,msg)     do{LOG4CXX_ERROR(::log4cxx::Logger::getLogger(nick),msg)}while(0)
#define FATALDUB(nick,msg)     do{LOG4CXX_FATAL(::log4cxx::Logger::getLogger(nick),msg)}while(0)

// Be very cautious about relying on "_ENABLED" macros across multiple ranks
// Correct design is to add a rank-specific Filter possibly to AppenderSkeleton
#define INFO0_ENABLED          (::logger::rankzero->isInfoEnabled())
#define TRACE0_ENABLED         (::logger::rankzero->isTraceEnabled())
#define DEBUG0_ENABLED         (::logger::rankzero->isDebugEnabled())
#define INFO_ENABLED           (::logger::all->isInfoEnabled())
#define TRACE_ENABLED          (::logger::all->isTraceEnabled())
#define DEBUG_ENABLED          (::logger::all->isDebugEnabled())
#define INFODUB_ENABLED(nick)  (::log4cxx::Logger::getLogger(nick)->isInfoEnabled())
#define TRACEDUB_ENABLED(nick) (::log4cxx::Logger::getLogger(nick)->isTraceEnabled())
#define DEBUGDUB_ENABLED(nick) (::log4cxx::Logger::getLogger(nick)->isDebugEnabled())

#endif // LOGGER_HPP
