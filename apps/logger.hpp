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

// Global logging instance usable by anyone
extern log4cxx::LoggerPtr logger;

// Logging macros to hide log4cxx specifics
#define TRACE(expr)     LOG4CXX_TRACE(logger,expr)
#define DEBUG(expr)     LOG4CXX_DEBUG(logger,expr)
#define INFO(expr)      LOG4CXX_INFO( logger,expr)
#define WARN(expr)      LOG4CXX_WARN( logger,expr)
#define ERROR(expr)     LOG4CXX_ERROR(logger,expr)
#define FATAL(expr)     LOG4CXX_FATAL(logger,expr)
#define INFO_ENABLED    (logger->isInfoEnabled())
#define TRACE_ENABLED   (logger->isTraceEnabled())
#define DEBUG_ENABLED   (logger->isDebugEnabled())

// Logging macros to log only from MPI rank 0
// Not intended for use within inner loops
#define TRACE0(expr) \
    do{ if (!::suzerain::mpi::comm_rank(MPI_COMM_WORLD)) TRACE(expr); }while(0)
#define DEBUG0(expr) \
    do{ if (!::suzerain::mpi::comm_rank(MPI_COMM_WORLD)) DEBUG(expr); }while(0)
#define INFO0(expr)  \
    do{ if (!::suzerain::mpi::comm_rank(MPI_COMM_WORLD)) INFO( expr); }while(0)
#define WARN0(expr)  \
    do{ if (!::suzerain::mpi::comm_rank(MPI_COMM_WORLD)) WARN( expr); }while(0)
#define ERROR0(expr) \
    do{ if (!::suzerain::mpi::comm_rank(MPI_COMM_WORLD)) ERROR(expr); }while(0)
#define FATAL0(expr) \
    do{ if (!::suzerain::mpi::comm_rank(MPI_COMM_WORLD)) FATAL(expr); }while(0)

void name_logger_within_comm_world();

#endif // LOGGER_HPP
