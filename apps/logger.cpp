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
 * logger.cpp: logging tools build atop log4cxx
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/mpi.hpp>
#include "logger.hpp"

// Will likely be replaced with a more descriptive instance in main()
log4cxx::LoggerPtr logger = log4cxx::Logger::getRootLogger();

// Failure to destroy logger prior to running static instance destructors
// causes segfaults.  Concoct an atexit callback specifically to destroy
// anything pointed to by logger prior to static instance destructors. 
static struct DestructLoggerRegistration {
    DestructLoggerRegistration() { atexit(&destruct_logger); }
    static void destruct_logger() { logger = 0; }
} destructLoggerRegistration;

void name_logger_within_comm_world()
{
    logger = log4cxx::Logger::getLogger(
            suzerain::mpi::comm_rank_identifier(MPI_COMM_WORLD));

    // Log only warnings and above from ranks 1 and higher when not debugging
    if (    suzerain::mpi::comm_rank(MPI_COMM_WORLD) > 0
         && !logger->isDebugEnabled()) {
        logger->setLevel(log4cxx::Level::getWarn());
    }
}
