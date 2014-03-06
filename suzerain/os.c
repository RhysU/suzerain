/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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

/** @file
 * @copydoc os.h
 */

#include <suzerain/os.h>

#include <suzerain/common.h>
#include <suzerain/countof.h>
#include <suzerain/error.h>

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

int suzerain_signal_number(const char * name)
{
    // Avoid issues with NULLs and zero-length values
    if (!name || !*name) return SUZERAIN_FAILURE;

    // Skip past any leading "SIG" prefix
    if (strncasecmp(name, "SIG", 3) == 0) name += 3;

    // Use string length to reduce the number of comparisons necessary
    switch (strlen(name)) {
    case 2:
#ifdef SIGIO
        if (strcasecmp(name, "IO") == 0) return SIGIO;
#endif
        return SUZERAIN_FAILURE;

    case 3:
#ifdef SIGBUS
        if (strcasecmp(name, "BUS") == 0) return SIGBUS;
#endif
#ifdef SIGCLD
        if (strcasecmp(name, "CLD") == 0) return SIGCLD;
#endif
#ifdef SIGFPE
        if (strcasecmp(name, "FPE") == 0) return SIGFPE;
#endif
#ifdef SIGHUP
        if (strcasecmp(name, "HUP") == 0) return SIGHUP;
#endif
#ifdef SIGILL
        if (strcasecmp(name, "ILL") == 0) return SIGILL;
#endif
#ifdef SIGINT
        if (strcasecmp(name, "INT") == 0) return SIGINT;
#endif
#ifdef SIGIOT
        if (strcasecmp(name, "IOT") == 0) return SIGIOT;
#endif
#ifdef SIGPWR
        if (strcasecmp(name, "PWR") == 0) return SIGPWR;
#endif
#ifdef SIGSYS
        if (strcasecmp(name, "SYS") == 0) return SIGSYS;
#endif
#ifdef SIGURG
        if (strcasecmp(name, "URG") == 0) return SIGURG;
#endif
        return SUZERAIN_FAILURE;

    case 4:
#ifdef SIGABRT
        if (strcasecmp(name, "ABRT") == 0) return SIGABRT;
#endif
#ifdef SIGALRM
        if (strcasecmp(name, "ALRM") == 0) return SIGALRM;
#endif
#ifdef SIGCHLD
        if (strcasecmp(name, "CHLD") == 0) return SIGCHLD;
#endif
#ifdef SIGCONT
        if (strcasecmp(name, "CONT") == 0) return SIGCONT;
#endif
#ifdef SIGKILL
        if (strcasecmp(name, "KILL") == 0) return SIGKILL;
#endif
#ifdef SIGPIPE
        if (strcasecmp(name, "PIPE") == 0) return SIGPIPE;
#endif
#ifdef SIGPOLL
        if (strcasecmp(name, "POLL") == 0) return SIGPOLL;
#endif
#ifdef SIGPROF
        if (strcasecmp(name, "PROF") == 0) return SIGPROF;
#endif
#ifdef SIGQUIT
        if (strcasecmp(name, "QUIT") == 0) return SIGQUIT;
#endif
#ifdef SIGSEGV
        if (strcasecmp(name, "SEGV") == 0) return SIGSEGV;
#endif
#ifdef SIGSTOP
        if (strcasecmp(name, "STOP") == 0) return SIGSTOP;
#endif
#ifdef SIGTERM
        if (strcasecmp(name, "TERM") == 0) return SIGTERM;
#endif
#ifdef SIGTRAP
        if (strcasecmp(name, "TRAP") == 0) return SIGTRAP;
#endif
#ifdef SIGTSTP
        if (strcasecmp(name, "TSTP") == 0) return SIGTSTP;
#endif
#ifdef SIGTTIN
        if (strcasecmp(name, "TTIN") == 0) return SIGTTIN;
#endif
#ifdef SIGTTOU
        if (strcasecmp(name, "TTOU") == 0) return SIGTTOU;
#endif
#ifdef SIGUSR1
        if (strcasecmp(name, "USR1") == 0) return SIGUSR1;
#endif
#ifdef SIGUSR2
        if (strcasecmp(name, "USR2") == 0) return SIGUSR2;
#endif
#ifdef SIGXCPU
        if (strcasecmp(name, "XCPU") == 0) return SIGXCPU;
#endif
#ifdef SIGXFSZ
        if (strcasecmp(name, "XFSZ") == 0) return SIGXFSZ;
#endif
        return SUZERAIN_FAILURE;

    case 5:
#ifdef SIGWINCH
        if (strcasecmp(name, "WINCH") == 0) return SIGWINCH;
#endif
        return SUZERAIN_FAILURE;

    case 6:
#ifdef SIGSTKFLT
        if (strcasecmp(name, "STKFLT") == 0) return SIGSTKFLT;
#endif
#ifdef SIGVTALRM
        if (strcasecmp(name, "VTALRM") == 0) return SIGVTALRM;
#endif
        return SUZERAIN_FAILURE;

    default:
        return SUZERAIN_FAILURE;
    }
}

const char * suzerain_signal_name(int signum)
{
    switch (signum) {
#ifdef SIGHUP
    case SIGHUP: return "SIGHUP";
#endif
#ifdef SIGINT
    case SIGINT: return "SIGINT";
#endif
#ifdef SIGQUIT
    case SIGQUIT: return "SIGQUIT";
#endif
#ifdef SIGILL
    case SIGILL: return "SIGILL";
#endif
#ifdef SIGTRAP
    case SIGTRAP: return "SIGTRAP";
#endif
#ifdef SIGABRT
    case SIGABRT: return "SIGABRT";
#endif
#if defined SIGIOT && defined SIGABRT && SIGIOT != SIGABRT
    case SIGIOT:
        return "SIGIOT";
#endif
#ifdef SIGBUS
    case SIGBUS: return "SIGBUS";
#endif
#ifdef SIGFPE
    case SIGFPE: return "SIGFPE";
#endif
#ifdef SIGKILL
    case SIGKILL: return "SIGKILL";
#endif
#ifdef SIGUSR1
    case SIGUSR1: return "SIGUSR1";
#endif
#ifdef SIGSEGV
    case SIGSEGV: return "SIGSEGV";
#endif
#ifdef SIGUSR2
    case SIGUSR2: return "SIGUSR2";
#endif
#ifdef SIGPIPE
    case SIGPIPE: return "SIGPIPE";
#endif
#ifdef SIGALRM
    case SIGALRM: return "SIGALRM";
#endif
#ifdef SIGTERM
    case SIGTERM: return "SIGTERM";
#endif
#ifdef SIGSTKFLT
    case SIGSTKFLT: return "SIGSTKFLT";
#endif
#if defined SIGCLD && defined SIGCHLD && SIGCLD != SIGCHLD
    case SIGCLD: return "SIGCLD";
#endif
#ifdef SIGCHLD
    case SIGCHLD: return "SIGCHLD";
#endif
#ifdef SIGCONT
    case SIGCONT: return "SIGCONT";
#endif
#ifdef SIGSTOP
    case SIGSTOP: return "SIGSTOP";
#endif
#ifdef SIGTSTP
    case SIGTSTP: return "SIGTSTP";
#endif
#ifdef SIGTTIN
    case SIGTTIN: return "SIGTTIN";
#endif
#ifdef SIGTTOU
    case SIGTTOU: return "SIGTTOU";
#endif
#ifdef SIGURG
    case SIGURG: return "SIGURG";
#endif
#ifdef SIGXCPU
    case SIGXCPU: return "SIGXCPU";
#endif
#ifdef SIGXFSZ
    case SIGXFSZ: return "SIGXFSZ";
#endif
#ifdef SIGVTALRM
    case SIGVTALRM: return "SIGVTALRM";
#endif
#ifdef SIGPROF
    case SIGPROF: return "SIGPROF";
#endif
#ifdef SIGWINCH
    case SIGWINCH: return "SIGWINCH";
#endif
#if defined SIGPOLL && defined SIGIO && SIGPOLL != SIGIO
    case SIGPOLL: return "SIGPOLL";
#endif
#ifdef SIGIO
    case SIGIO: return "SIGIO";
#endif
#ifdef SIGPWR
    case SIGPWR: return "SIGPWR";
#endif
#ifdef SIGSYS
    case SIGSYS: return "SIGSYS";
#endif
    default: return NULL;
    }
}

const char * suzerain_temporary_directory() {

    const char *var[] = { "TMPDIR", "TMP", "TEMPDIR", "TEMP" };

    const char * s = NULL;
    size_t i;
    for (i = 0; i < SUZERAIN_COUNTOF(var); ++i) {
        s = getenv(var[i]);               // Retrieve from environment
        if (s) {
            while (isspace(s[0])) ++s;    // Skip any leading whitespace
            if (!s[0]) s = NULL;          // Skip zero-length answers
        }
        if (s) break;                     // Stop on something good
    }
    if (!s) s = "/tmp";                   // Default if nothing good

    return s;
}
