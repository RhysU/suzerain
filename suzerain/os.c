/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * os.c: POSIX-based operating system utilities
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <suzerain/config.h>
#include <suzerain/common.h>
#include <suzerain/error.h>
#pragma hdrstop
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int suzerain_fpipe(int *rfd, int rflags, FILE **w, int wflags)
{
    int pipefd[2]; /* read/write are pipefd[0/1] respectively */
    int flags;

    /* Reset errno */
    errno = 0;

    /* Create the pipe */
    if (pipe(pipefd) == -1) {
        SUZERAIN_ERROR_VAL(strerror(errno), SUZERAIN_FAILURE, errno);
    }

    /* Set the flags on the read side */
    flags = fcntl(pipefd[0], F_GETFL, 0);
    if (flags == -1) {
        SUZERAIN_ERROR_VAL(strerror(errno), SUZERAIN_FAILURE, errno);
    }
    if (fcntl(pipefd[0], F_SETFL, flags | rflags) == -1) {
        SUZERAIN_ERROR_VAL(strerror(errno), SUZERAIN_FAILURE, errno);
    }

    /* Set the flags on the write side */
    flags = fcntl(pipefd[1], F_GETFL, 1);
    if (flags == -1) {
        SUZERAIN_ERROR_VAL(strerror(errno), SUZERAIN_FAILURE, errno);
    }
    if (fcntl(pipefd[1], F_SETFL, flags | wflags) == -1) {
        SUZERAIN_ERROR_VAL(strerror(errno), SUZERAIN_FAILURE, errno);
    }

    /* Get the file descriptor for the read side */
    *rfd = pipefd[0];

    /* Get the file descriptor for the write side */
    *w = fdopen(pipefd[1], "a");
    if (*w == NULL) {
        SUZERAIN_ERROR_VAL(strerror(errno), SUZERAIN_FAILURE, errno);
    }

    return SUZERAIN_SUCCESS;
}
