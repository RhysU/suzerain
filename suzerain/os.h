/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 The PECOS Development Team
 *
 * This file is part of Suzerain.
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * os.h: POSIX-based operating system utilities
 * $Id$
 */

#ifndef __SUZERAIN_OS_H
#define __SUZERAIN_OS_H

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

/**
 * Given a signal name, e.g. "SIGTERM" or "TERM", return the signal's number
 * defined in \c signal.h.
 *
 * @param name Name of the signal to look up.
 *
 * @return The signal number on success or ::SUZERAIN_FAILURE on error.
 */
int suzerain_signal_number(const char * name);

/**
 * Given a signal number return a human-readable name following \c signal.h.
 * For example, 1 will return "SIGHUP" on many systems.
 *
 * @param signum Number of the signal to look up.
 *
 * @return The signal name on success or NULL on failure.
 */
const char * suzerain_signal_name(int signum);

/**
 * Obtain a "good" temporary directory name which is likely to be on
 * non-networked disk.  Uses the first of the following environment variables
 * found which, after leading spaces are removed, has non-zero length:
 * \li \c TMPDIR
 * \li \c TMP
 * \li \c TEMPDIR
 * \li \c TEMP
 * If none of these possibilities give a non-empty result, the method returns
 * the string "/tmp".
 *
 * @return The name of a temporary directory.
 */
const char * suzerain_temporary_directory();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_OS_H */
