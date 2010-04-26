/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * os.h: POSIX-based operating system utilities
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_OS_H__
#define __SUZERAIN_OS_H__

#include <suzerain/common.h>
#include <fcntl.h>

/** @file
 * Provides utilities atop POSIX methods.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Create a POSIX pipe, set the provided flags on the read and write
 * portions, and create a <tt>FILE*</tt> handle to the write end
 * suitable for use with \c fprintf and friends.
 *
 * @param[out] rfd File descriptor attached to the pipe's read side.
 * @param[in] rflags Flags to set for the new pipe's read side.
 * @param[out] w Handle return for the new pipe's write side.
 * @param[in] wflags Flags to set for the new pipe's write side.
 *
 * @see The manual pages for <tt>fcntl(2)</tt>, <tt>pipe(2)</tt>,
 *      and <tt>pipe(7)</tt> for more information on creating pipes
 *      and the effect of their flags.
 * @see The manual pages for <tt>close(2)</tt> and <tt>fclose(3)</tt>
 *      for information on how to close \c rfd and <tt>*write</tt> after
 *      you have finished using them.
 *
 * @return Zero on successful completion.  On non-zero returns \c
 *         SUZERAIN_FAILURE and \c errno is set appropriately.
 */
int suzerain_fpipe(int *rfd, int rflags, FILE **w, int wflags);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_OS_H__ */
