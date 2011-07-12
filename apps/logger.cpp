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
// Beware http://old.nabble.com/Static-destruction-fiasco--td31026705.html
log4cxx::LoggerPtr logger(log4cxx::Logger::getRootLogger());

void name_logger_within_comm_world()
{
    using log4cxx::LevelPtr;
    using log4cxx::Level;
    using log4cxx::Logger;
    using suzerain::mpi::comm_rank;
    using suzerain::mpi::comm_rank_identifier;

    logger = Logger::getLogger(comm_rank_identifier(MPI_COMM_WORLD));

    // If debugging is disabled, ensure level is at least INFO on ranks 1+
    LevelPtr level = logger->getLevel();
    if (comm_rank(MPI_COMM_WORLD) > 0 && !logger->isDebugEnabled()) {
        if (!level || !level->isGreaterOrEqual(Level::getInfo())) {
            level = Level::getInfo();
        }
    }
    logger->setLevel(level);
}
