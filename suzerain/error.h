/*--------------------------------------------------------------------------
 *
 * Functionality adopted from the GNU Scientific Library.
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * Copyright (C) 2010, 2011 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
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
 */

#ifndef SUZERAIN_ERROR_H
#define SUZERAIN_ERROR_H

/** @file
 * Provides standardized error numbers and error handling routines.  Error
 * reporting follows the design and conventions used in the <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL) <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Error-Handling.html">
 * error handling routines</a>.  Much of suzerain's error code is a direct copy
 * of GSL's API and source code.  Notable exceptions are the MPI error handling
 * macros which are an improved copy of ideas found in <a
 * href="http://www.mcs.anl.gov/petsc/">PETSc</a>.
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Standardized error status codes used throughout Suzerain.
 * Where possible these codes are numerically equivalent to
 * <a href="http://www.gnu.org/software/gsl/manual/html_node/Error-Codes.html">
 * GSL's error codes</a>.
 */
enum suzerain_error_status {
    SUZERAIN_FAILURE  = -1, /**< Failure */
    SUZERAIN_CONTINUE = -2, /**< Iteration has not converged */
    SUZERAIN_SUCCESS  =  0, /**< Success */
    SUZERAIN_EDOM     =  1, /**< Input domain error */
    SUZERAIN_ERANGE   =  2, /**< Output range error */
    SUZERAIN_EFAULT   =  3, /**< Invalid pointer */
    SUZERAIN_EINVAL   =  4, /**< Invalid argument supplied by user */
    SUZERAIN_EFAILED  =  5, /**< Generic failure */
    SUZERAIN_ESANITY  =  7, /**< Sanity check failed - shouldn't happen */
    SUZERAIN_ENOMEM   =  8, /**< Memory allocation failed */
    SUZERAIN_EBADFUNC =  9, /**< Problem with user-supplied function */
    SUZERAIN_EMAXITER = 11, /**< Exceeded max number of iterations */
    SUZERAIN_EZERODIV = 12, /**< Tried to divide by zero */
    SUZERAIN_EROUND   = 18, /**< Failed because of roundoff error */
    SUZERAIN_EBADLEN  = 19, /**< Matrix or vector lengths are not conformant */
    SUZERAIN_ESING    = 21, /**< Apparent singularity detected */
    SUZERAIN_EDIVERGE = 22, /**< Integral or series is divergent */
    SUZERAIN_EUNIMPL  = 24  /**< Requested feature not (yet) implemented */
};

/**
 * Calls the error handler last set using suzerain_set_error_handler
 * when invoked.  This is the entry point to the error handling system.
 *
 * The default behavior is to log the error to the stream specified using
 * suzerain_set_stream. The functions suzerain_set_stream,
 * suzerain_set_stream_handler, and suzerain_set_error_handler can be used to
 * modify this behavior.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param suzerain_errno Error code to report.  Should be one of
 *      suzerain_error_status if at all possible.
 *
 * @see Most clients should not use this function directly; instead use one of
 *      the convenience macros: SUZERAIN_ERROR, SUZERAIN_ERROR_VAL,
 *      SUZERAIN_ERROR_VOID, SUZERAIN_ERROR_NULL
 */
void suzerain_error(const char * reason,
                    const char * file,
                    int line,
                    int suzerain_errno);

/**
 * Print an error message to the current error stream.
 * If a suzerain_stream_handler_t has been specified, it is used.
 * If a stream has been set using suzerain_set_stream, it is used.
 * Lastly, the routine prints the error message to standard error.
 *
 * @param label Label used to identify the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param reason Reason for the error.
 */
void suzerain_stream_printf(const char *label,
                            const char *file,
                            int line,
                            const char *reason);

/**
 * Look up a human-readable error message for the given error status.
 *
 * @param suzerain_errno Error code to look up.
 *
 * @return A message suitable for use in logging or error messages.
 */
const char * suzerain_strerror(const int suzerain_errno);

/**
 * Defines the function prototype necessary for an error handler.
 * Error handlers should be reentrant safe if possible.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param suzerain_errno Error code to report.
 *
 * @see suzerain_set_error_handler
 */
typedef void suzerain_error_handler_t(const char * reason,
                                      const char * file,
                                      int line,
                                      int suzerain_errno);

/**
 * Defines the function prototype necessary for a stream handler.
 * Stream handlers should be reentrant safe if possible.
 *
 * @param reason Reason for the error.
 * @param file File in which the error was reported.
 * @param line Line at which the error was reported.
 * @param suzerain_errno Error code to report.
 *
 * @see suzerain_set_stream_handler
 */
typedef void suzerain_stream_handler_t(const char * label,
                                       const char * file,
                                       int line,
                                       const char * reason);

/**
 * Sets the current error handler for the process.
 * Invoked by suzerain_error when an error occurs.
 *
 * @param new_handler New error handler to use.
 *
 * @return the previous error handler in use.
 */
suzerain_error_handler_t *
suzerain_set_error_handler(suzerain_error_handler_t * new_handler);

/**
 * An error handler implementation that disables all error reporting.
 * Primarily intended for use in test environments.
 *
 * @return the previous error handler in use.
 */
suzerain_error_handler_t *
suzerain_set_error_handler_off(void);

/**
 * Sets the current stream handler for the process.
 * Used by the default error handling behavior, and possibly by
 * other custom error handling routines.
 *
 * @param new_handler New stream handler to use.
 *
 * @return the previous stream handler in use.
 */
suzerain_stream_handler_t *
suzerain_set_stream_handler(suzerain_stream_handler_t * new_handler);

/**
 * Set the default stream for error message display.  Default
 * behavior is to use stderr.
 *
 * @param new_stream New stream to use.
 *
 * @return the previous stream in use.
 */
FILE * suzerain_set_stream(FILE * new_stream);

/**
 * Invokes suzerain_error and returns the value \c suzerain_errno.
 * Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param suzerain_errno Error status to report and returned from the current
 *      function.
 */
#define SUZERAIN_ERROR(reason, suzerain_errno) \
       do { \
       suzerain_error (reason, __FILE__, __LINE__, suzerain_errno) ; \
       return suzerain_errno ; \
       } while (0)

/**
 * Invokes suzerain_error using \c suzerain_errno and returns the value \c
 * value.  Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param suzerain_errno Error status to report.
 * @param value Value to return from the current function.
 */
#define SUZERAIN_ERROR_VAL(reason, suzerain_errno, value) \
       do { \
       suzerain_error (reason, __FILE__, __LINE__, suzerain_errno) ; \
       return value ; \
       } while (0)

/**
 * Invokes suzerain_error using \c suzerain_errno and returns from the current
 * function.  Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param suzerain_errno Error status to report.
 */
#define SUZERAIN_ERROR_VOID(reason, suzerain_errno) \
       do { \
       suzerain_error (reason, __FILE__, __LINE__, suzerain_errno) ; \
       return ; \
       } while (0)

/**
 * Invokes suzerain_error using \c suzerain_errno and returns NULL
 * from the current function.  Useful for out-of-memory conditions.
 * Automatically provides file and line information.
 *
 * @param reason Message to report.
 * @param suzerain_errno Error status to report.
 */
#define SUZERAIN_ERROR_NULL(reason, suzerain_errno) \
        SUZERAIN_ERROR_VAL(reason, suzerain_errno, 0)

/**
 * Invokes suzerain_error using \c suzerain_errno but \em does \em not return
 * from the current function.  Automatically provides file and line
 * information.
 *
 * @param reason Message to report.
 * @param suzerain_errno Error status to report.
 */
#define SUZERAIN_ERROR_REPORT(reason, suzerain_errno) \
       do { \
       suzerain_error (reason, __FILE__, __LINE__, suzerain_errno) ; \
       } while (0)

/** Shorthand for using \c SUZERAIN_ERROR_VAL for unimplemented logic. */
#define SUZERAIN_ERROR_VAL_UNIMPLEMENTED(value) \
        SUZERAIN_ERROR_VAL("Unimplemented logic!", SUZERAIN_EUNIMPL, value)

/** Shorthand for using \c SUZERAIN_ERROR_VOID for unimplemented logic. */
#define SUZERAIN_ERROR_VOID_UNIMPLEMENTED() \
        SUZERAIN_ERROR_VOID("Unimplemented logic!", SUZERAIN_EUNIMPL)

/** Shorthand for using \c SUZERAIN_ERROR_NULL for unimplemented logic. */
#define SUZERAIN_ERROR_NULL_UNIMPLEMENTED() \
        SUZERAIN_ERROR_NULL("Unimplemented logic!", SUZERAIN_EUNIMPL)

/** Shorthand for using \c SUZERAIN_ERROR_REPORT for unimplemented logic. */
#define SUZERAIN_ERROR_REPORT_UNIMPLEMENTED() \
        SUZERAIN_ERROR_REPORT("Unimplemented logic!", SUZERAIN_EUNIMPL)

#ifndef SUZERAIN_PARSED_BY_DOXYGEN
/* Internal helper macro for implementing SUZERAIN_MPICHKx macros */
#define SUZERAIN_MPICHKx_TEMPLATE(suzerain_error_macro,stmt) \
    do { \
        const int _chk_stat = (stmt); \
        if (_chk_stat != MPI_SUCCESS) { \
            char _chk_reason[384]; \
            char _chk_mpistring[MPI_MAX_ERROR_STRING]; \
            int _chk_len; \
            const int _chk_string_stat \
                = MPI_Error_string(_chk_stat,_chk_mpistring,&_chk_len); \
            snprintf(_chk_reason, sizeof(_chk_reason), \
                    "Encountered MPI error code %d: %s", _chk_stat, \
                    (_chk_string_stat == MPI_SUCCESS) \
                    ? _chk_mpistring : "UNKNOWN"); \
            suzerain_error_macro(_chk_reason, SUZERAIN_EFAILED); \
        } \
    } while(0)
#endif /* SUZERAIN_PARSED_BY_DOXYGEN */

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * SUZERAIN_ERROR.  Any relevant message is looked up using \c MPI_Error_string
 * and reported.  \c SUZERAIN_EFAILED is the return value provided to \c
 * SUZERAIN_ERROR.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>suzerain/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRQ for the original inspiration for this macro.
 */
#define SUZERAIN_MPICHKQ(stmt) \
    SUZERAIN_MPICHKx_TEMPLATE(SUZERAIN_ERROR,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * SUZERAIN_ERROR_NULL.  Any relevant message is looked up using \c
 * MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>suzerain/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRQ for the original inspiration for this macro.
 */
#define SUZERAIN_MPICHKN(stmt) \
    SUZERAIN_MPICHKx_TEMPLATE(SUZERAIN_ERROR_NULL,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * SUZERAIN_ERROR_VOID.  Any relevant message is looked up using \c
 * MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>suzerain/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRV for the original inspiration for this macro.
 */
#define SUZERAIN_MPICHKV(stmt) \
    SUZERAIN_MPICHKx_TEMPLATE(SUZERAIN_ERROR_VOID,stmt)

/**
 * Executes \c stmt once handling any resulting MPI error per \c
 * SUZERAIN_ERROR_REPORT. The current function \em continues \em executing. Any
 * relevant message is looked up using \c MPI_Error_string and reported.
 *
 * @param stmt Statement, presumably an MPI call, to be executed.
 * @note <tt>suzerain/mpi.h</tt> must be <tt>include</tt>d for the macro
 *       expansion to compile correctly.
 * @warning This macro expands to a not insignificant amount of code.
 *          It should not be used in performance critical regions.
 * @see PETSc's CHKERRV for the original inspiration for this macro.
 */
#define SUZERAIN_MPICHKR(stmt) \
    SUZERAIN_MPICHKx_TEMPLATE(SUZERAIN_ERROR_REPORT,stmt)

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_ERROR_H */
