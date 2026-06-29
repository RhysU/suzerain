/*--------------------------------------------------------------------------
 *
 * Code copyright and licensing details follow on routine-specific basis.
 *
 *--------------------------------------------------------------------------
 */

/** @file
 * @copydoc pre_gsl.h
 */

#include <suzerain/pre_gsl.h>

#include <gsl/gsl_ieee_utils.h>
#include <stdio.h>
#include <unistd.h>

#include <suzerain/common.h>

/*
 * License identical to gsl_ieee_env_setup() from GSL.
 * Copyright (c) 2012-2014 Rhys Ulerich
 */
void
mpi_gsl_ieee_env_setup(const int rank)
{
    // Used for write(2)-based error messages
    static const char newline[] = "\n";

    // Flush stdout, stderr
    if (fflush(stdout))
        perror("mpi_gsl_ieee_env_setup fflush(stdout) before redirect");
    if (fflush(stderr))
        perror("mpi_gsl_ieee_env_setup fflush(stderr) before redirect");

    // Save stdout, stderr so we may restore them later
    int stdout_copy, stderr_copy;
    if ((stdout_copy = dup(fileno(stdout))) < 0)
        perror("mpi_gsl_ieee_env_setup error duplicating stdout");
    if ((stderr_copy = dup(fileno(stderr))) < 0)
        perror("mpi_gsl_ieee_env_setup error duplicating stderr");

    // On non-root processes redirect stdout, stderr to /dev/null
    if (rank) {
        if (!freopen("/dev/null", "a", stdout))
            perror("mpi_gsl_ieee_env_setup redirecting stdout");
        if (!freopen("/dev/null", "a", stderr))
            perror("mpi_gsl_ieee_env_setup redirecting stderr");
    }

    // Invoke gsl_ieee_env_setup on all ranks.
    gsl_ieee_env_setup();

    // Flush stdout, stderr again.
    // Error messages sent to stderr_copy on a best-effort basis
    if (fflush(stdout)) {
        static const char pre[] = "mpi_gsl_ieee_env_setup fflush(stdout) after redirect: ";
        const char *msg = strerror(errno);
        fsync(stderr_copy);
        write(stderr_copy, pre,     sizeof(pre));
        write(stderr_copy, msg,     strlen(msg));
        write(stderr_copy, newline, strlen(newline));
        fsync(stderr_copy);
    }
    if (fflush(stderr)) {
        static const char pre[] = "mpi_gsl_ieee_env_setup fflush(stderr) after redirect: ";
        const char *msg = strerror(errno);
        fsync(stderr_copy);
        write(stderr_copy, pre,     sizeof(pre));
        write(stderr_copy, msg,     strlen(msg));
        write(stderr_copy, newline, strlen(newline));
        fsync(stderr_copy);
    }

    // Restore stdout, stderr
    // Error messages sent to stderr_copy on a best-effort basis
    if (dup2(stdout_copy, fileno(stdout)) < 0) {
        static const char pre[] = "mpi_gsl_ieee_env_setup reopening stdout: ";
        const char *msg = strerror(errno);
        fsync(stderr_copy);
        write(stderr_copy, pre,     sizeof(pre));
        write(stderr_copy, msg,     strlen(msg));
        write(stderr_copy, newline, strlen(newline));
        fsync(stderr_copy);
    }
    if (dup2(stderr_copy, fileno(stderr)) < 0) {
        static const char pre[] = "mpi_gsl_ieee_env_setup reopening stderr: ";
        const char *msg = strerror(errno);
        fsync(stderr_copy);
        write(stderr_copy, pre,     sizeof(pre));
        write(stderr_copy, msg,     strlen(msg));
        write(stderr_copy, newline, strlen(newline));
        fsync(stderr_copy);
    }

    // Close saved versions of stdout, stderr
    if (close(stdout_copy))
        perror("mpi_gsl_ieee_env_setup closing stdout_copy");
    if (close(stderr_copy))
        perror("mpi_gsl_ieee_env_setup closing stderr_copy");

    // Clear any errors that may have occurred on stdout, stderr
    clearerr(stdout);
    clearerr(stderr);
}
