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

namespace detail {

// Global logging instance
extern ::log4cxx::LoggerPtr loggerptr;

// World rank determined within initialize method
extern int worldrank;

} // end namespace detail

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

// Logging macros to hide log4cxx specifics
#define TRACE(expr)     LOG4CXX_TRACE(::logger::detail::loggerptr,expr)
#define DEBUG(expr)     LOG4CXX_DEBUG(::logger::detail::loggerptr,expr)
#define INFO(expr)      LOG4CXX_INFO( ::logger::detail::loggerptr,expr)
#define WARN(expr)      LOG4CXX_WARN( ::logger::detail::loggerptr,expr)
#define ERROR(expr)     LOG4CXX_ERROR(::logger::detail::loggerptr,expr)
#define FATAL(expr)     LOG4CXX_FATAL(::logger::detail::loggerptr,expr)
#define INFO_ENABLED    (::logger::detail::loggerptr->isInfoEnabled())
#define TRACE_ENABLED   (::logger::detail::loggerptr->isTraceEnabled())
#define DEBUG_ENABLED   (::logger::detail::loggerptr->isDebugEnabled())

// Logging macros to log only from MPI rank 0
#define TRACE0(expr) \
    do{ if (!::logger::detail::worldrank) TRACE(expr); }while(0)
#define DEBUG0(expr) \
    do{ if (!::logger::detail::worldrank) DEBUG(expr); }while(0)
#define INFO0(expr)  \
    do{ if (!::logger::detail::worldrank) INFO( expr); }while(0)
#define WARN0(expr)  \
    do{ if (!::logger::detail::worldrank) WARN( expr); }while(0)
#define ERROR0(expr) \
    do{ if (!::logger::detail::worldrank) ERROR(expr); }while(0)
#define FATAL0(expr) \
    do{ if (!::logger::detail::worldrank) FATAL(expr); }while(0)

#endif // LOGGER_HPP
