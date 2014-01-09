/***************************************************************************
// Adapted from the GNU Scientific Library
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$
***************************************************************************/

#ifndef UNDERLING_ERROR_H
#define UNDERLING_ERROR_H

#include <stdio.h>
#include <underling/visibility.h>

/** @file
 * Provides standardized error numbers and error handling routines.  Error
 * reporting follows the design and conventions used in the <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL) <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Error-Handling.html">
 * error handling routines</a>.  Much of underling's error handling is a direct
 * copy of GSL's API and source code.  Notable exceptions are the MPI error
 * handling macros which are an improved copy of ideas found in <a
 * href="http://www.mcs.anl.gov/petsc/">PETSc</a>.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Standardized error status codes used throughout underling.
 * Where possible these codes are numerically equivalent to
 * <a href="http://www.gnu.org/software/gsl/manual/html_node/Error-Codes.html">
 * GSL's error codes</a>.
 *
 * Note that \ref UNDERLING_SUCCESS is zero to allow code like
 * <code>if (!status) { some_error_handling() }</code>.
 */
enum underling_status {
    UNDERLING_SUCCESS  =  0, /**< Success */
    UNDERLING_EFAULT   =  3, /**< Invalid pointer */
    UNDERLING_EINVAL   =  4, /**< Invalid argument supplied by user */
    UNDERLING_EFAILED  =  5, /**< Generic failure */
    UNDERLING_ESANITY  =  7, /**< Sanity check failed - shouldn't happen */
    UNDERLING_ENOMEM   =  8  /**< Memory allocation failed */
};

/**
 * Calls the error handler last set using underling_set_error_handler
 * when invoked.  This is the entry point to the error handling system.
 *
 * The default behavior is to log the error to the stream specified using
 * underling_set_stream. The functions underling_set_stream,
 * underling_set_stream_handler, and underling_set_error_handler can be used to
 * modify this behavior.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param underling_errno Error code to report.  Should be one of
 *      underling_status if at all possible.
 *
 * @see Most clients should not use this function directly; instead use one of
 *      the convenience macros: UNDERLING_ERROR, UNDERLING_ERROR_VAL,
 *      UNDERLING_ERROR_VOID, UNDERLING_ERROR_NULL
 */
void
underling_error(const char * reason,
                const char * file,
                int line,
                int underling_errno) UNDERLING_API;

/**
 * Print an error message to the current error stream.
 * If a underling_stream_handler_t has been specified, it is used.
 * If a stream has been set using underling_set_stream, it is used.
 * Lastly, the routine prints the error message to standard error.
 *
 * @param label Label used to identify the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param reason Reason for the error.
 */
void
underling_stream_printf(const char *label,
                        const char *file,
                        int line,
                        const char *reason) UNDERLING_API;

/**
 * Look up a human-readable error message for the given error status.
 *
 * @param underling_errno Error code to look up.
 *
 * @return A message suitable for use in logging or error messages.
 */
const char *
underling_strerror(const int underling_errno) UNDERLING_API;

/**
 * Defines the function prototype necessary for an error handler.
 * Error handlers should be reentrant safe if possible.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param underling_errno Error code to report.
 *
 * @see underling_set_error_handler
 */
typedef void underling_error_handler_t(const char * reason,
                                       const char * file,
                                       int line,
                                       int underling_errno);

/**
 * Defines the function prototype necessary for a stream handler.
 * Stream handlers should be reentrant safe if possible.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param underling_errno Error code to report.
 *
 * @see underling_set_stream_handler
 */
typedef void underling_stream_handler_t(const char * label,
                                        const char * file,
                                        int line,
                                        const char * reason);

/**
 * Sets the current error handler for the process.
 * Invoked by underling_error when an error occurs.
 *
 * @param new_handler New error handler to use.
 *
 * @return the previous error handler in use.
 */
underling_error_handler_t *
underling_set_error_handler(underling_error_handler_t * new_handler)
    UNDERLING_API;

/**
 * An error handler implementation that disables all error reporting.
 * Primarily intended for use in test environments.
 *
 * @return the previous error handler in use.
 */
underling_error_handler_t *
underling_set_error_handler_off(void) UNDERLING_API;

/**
 * Sets the current stream handler for the process.
 * Used by the default error handling behavior, and possibly by
 * other custom error handling routines.
 *
 * @param new_handler New stream handler to use.
 *
 * @return the previous stream handler in use.
 */
underling_stream_handler_t *
underling_set_stream_handler(underling_stream_handler_t * new_handler)
    UNDERLING_API;

/**
 * Set the default stream for error message display.  Default
 * behavior is to use stderr.
 *
 * @param new_stream New stream to use.
 *
 * @return the previous stream in use.
 */
FILE *
underling_set_stream(FILE * new_stream) UNDERLING_API;

/**
 * Invokes underling_error and returns the value \c underling_errno.
 * Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param underling_errno Error status to report and returned from the current
 *      function.
 */
#define UNDERLING_ERROR(reason, underling_errno) \
       do { \
       underling_error (reason, __FILE__, __LINE__, underling_errno) ; \
       return underling_errno ; \
       } while (0)

/**
 * Invokes underling_error using \c underling_errno and returns the value \c
 * value.  Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param underling_errno Error status to report.
 * @param value Value to return from the current function.
 */
#define UNDERLING_ERROR_VAL(reason, underling_errno, value) \
       do { \
       underling_error (reason, __FILE__, __LINE__, underling_errno) ; \
       return value ; \
       } while (0)

/**
 * Invokes underling_error using \c underling_errno and returns from the
 * current function.  Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param underling_errno Error status to report.
 */
#define UNDERLING_ERROR_VOID(reason, underling_errno) \
       do { \
       underling_error (reason, __FILE__, __LINE__, underling_errno) ; \
       return ; \
       } while (0)

/**
 * Invokes underling_error using \c underling_errno and returns NULL
 * from the current function.  Useful for out-of-memory conditions.
 * Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param underling_errno Error status to report.
 */
#define UNDERLING_ERROR_NULL(reason, underling_errno) \
        UNDERLING_ERROR_VAL(reason, underling_errno, 0)

/**
 * Invokes underling_error using \c underling_errno but \em does \em not return
 * from the current function.  Automatically provides file and line
 * information.
 *
 * @param reason Message to report.
 * @param underling_errno Error status to report.
 */
#define UNDERLING_ERROR_REPORT(reason, underling_errno) \
       do { \
       underling_error (reason, __FILE__, __LINE__, underling_errno) ; \
       } while (0)

/** \cond INTERNAL */
/* Internal helper macro for implementing UNDERLING_MPICHKx macros */
#define UNDERLING_MPICHKx_TEMPLATE(underling_error_macro,stmt) \
    do { \
        const int _chk_stat = (stmt); \
        if (_chk_stat != MPI_SUCCESS) { \
            char _chk_reason[255]; \
            char _chk_mpistring[MPI_MAX_ERROR_STRING]; \
            int _chk_len; \
            const int _chk_string_stat \
                = MPI_Error_string(_chk_stat,_chk_mpistring,&_chk_len); \
            snprintf(_chk_reason, sizeof(_chk_reason)/sizeof(_chk_reason[0]), \
                    "Encountered MPI error code %d: %s", _chk_stat, \
                    (_chk_string_stat == MPI_SUCCESS) \
                    ? _chk_mpistring : "UNKNOWN"); \
            underling_error_macro(_chk_reason, UNDERLING_EFAILED); \
        } \
    } while(0)
/** \endcond */

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * UNDERLING_ERROR.  Any relevant message is looked up using \c
 * MPI_Error_string and reported.  \c UNDERLING_EFAILED is the return value
 * provided to \c UNDERLING_ERROR.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>underling/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRQ for the original inspiration for this macro.
 */
#define UNDERLING_MPICHKQ(stmt) \
    UNDERLING_MPICHKx_TEMPLATE(UNDERLING_ERROR,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * UNDERLING_ERROR_NULL.  Any relevant message is looked up using \c
 * MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>underling/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRQ for the original inspiration for this macro.
 */
#define UNDERLING_MPICHKN(stmt) \
    UNDERLING_MPICHKx_TEMPLATE(UNDERLING_ERROR_NULL,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * UNDERLING_ERROR_VOID.  Any relevant message is looked up using \c
 * MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>underling/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRV for the original inspiration for this macro.
 */
#define UNDERLING_MPICHKV(stmt) \
    UNDERLING_MPICHKx_TEMPLATE(UNDERLING_ERROR_VOID,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * UNDERLING_ERROR_REPORT. The current function \em continues \em executing.
 * Any relevant message is looked up using \c MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>underling/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRV for the original inspiration for this macro.
 */
#define UNDERLING_MPICHKR(stmt) \
    UNDERLING_MPICHKx_TEMPLATE(UNDERLING_ERROR_REPORT,stmt)

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* UNDERLING_ERROR_H__ */
