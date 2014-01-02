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

#include <suzerain/os.h>

#include <suzerain/common.h>

static
int
test_suzerain_fpipe( void )
{
    int rfd;
    FILE *w;
    int error;

    error = suzerain_fpipe(&rfd, O_NONBLOCK, &w, O_NONBLOCK);
    if (error) {
        fprintf(stderr,
                "suzerain_fpipe returned %d: %s\n", error, strerror(errno));
        return EXIT_FAILURE;
    }

    const char msg[] = "Mr. Watson -- come here -- I want to see you.\n";
    const int  len   = strlen(msg);

    /* Call fprintf on the write side of the buffer */
    error = fprintf(w, msg);
    if (error < 0) {
        fprintf(stderr,
                "fprintf reported error %d writing on suzerain_fpipe\n",
                error);
        return EXIT_FAILURE;
    }
    error = fflush(w);
    if (error) {
        fprintf(stderr,
                "fflush reported error %d writing on suzerain_fpipe: %s\n",
                error, strerror(errno));
        return EXIT_FAILURE;
    }
    error = fclose(w);
    if (error) {
        fprintf(stderr,
                "fclose reported %d closing suzerain_fpipe write side: %s\n",
                error, strerror(errno));
        return EXIT_FAILURE;
    }

    /* Check that the message came across on the read side */
    char buffer[len];
    int pos = 0;
    while (pos < len) {
        const ssize_t amount = read(rfd, buffer, len - pos);
        if (amount <= 0) {
            break;
        }
        pos += amount;
    }
    if ((pos != len) || (0 != strncmp(msg, buffer, len))) {
        fprintf(stderr, "Received incorrect data from suzerain_fpipe\n");
        fprintf(stderr, "Expected: (%d)[[[%s]]]\n", len, msg);
        buffer[len-1] = 0;
        fprintf(stderr, "Received: (%d)[[[%s]]]\n", pos, buffer);
        return EXIT_FAILURE;
    }
    error = close(rfd);
    if (error) {
        fprintf(stderr,
                "fclose reported %d closing suzerain_fpipe write side: %s",
                error, strerror(errno));
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

static
int
test_suzerain_signal_number( void )
{
    int success = 1;

    // The ISO C standard only requires the signal names SIGABRT, SIGFPE,
    // SIGILL, SIGINT, SIGSEGV, and SIGTERM to be defined.

    success |= (suzerain_signal_number("SIGABRT") == SIGABRT);
    success |= (suzerain_signal_number("sigabrt") == SIGABRT);
    success |= (suzerain_signal_number("ABRT") == SIGABRT);
    success |= (suzerain_signal_number("Abrt") == SIGABRT);
    success |= (suzerain_signal_number("abrt") == SIGABRT);

    success |= (suzerain_signal_number("SIGFPE") == SIGFPE);
    success |= (suzerain_signal_number("sigfpe") == SIGFPE);
    success |= (suzerain_signal_number("FPE") == SIGFPE);
    success |= (suzerain_signal_number("Fpe") == SIGFPE);
    success |= (suzerain_signal_number("fpe") == SIGFPE);

    success |= (suzerain_signal_number("SIGILL") == SIGILL);
    success |= (suzerain_signal_number("sigill") == SIGILL);
    success |= (suzerain_signal_number("ILL") == SIGILL);
    success |= (suzerain_signal_number("Ill") == SIGILL);
    success |= (suzerain_signal_number("ill") == SIGILL);

    success |= (suzerain_signal_number("SIGINT") == SIGINT);
    success |= (suzerain_signal_number("sigINT") == SIGINT);
    success |= (suzerain_signal_number("INT") == SIGINT);
    success |= (suzerain_signal_number("Int") == SIGINT);
    success |= (suzerain_signal_number("INT") == SIGINT);

    success |= (suzerain_signal_number("SIGSEGV") == SIGSEGV);
    success |= (suzerain_signal_number("sigsegv") == SIGSEGV);
    success |= (suzerain_signal_number("SEGV") == SIGSEGV);
    success |= (suzerain_signal_number("Segv") == SIGSEGV);
    success |= (suzerain_signal_number("segv") == SIGSEGV);

    success |= (suzerain_signal_number("SIGTERM") == SIGTERM);
    success |= (suzerain_signal_number("sigterm") == SIGTERM);
    success |= (suzerain_signal_number("TERM") == SIGTERM);
    success |= (suzerain_signal_number("Term") == SIGTERM);
    success |= (suzerain_signal_number("term") == SIGTERM);

    return success == 1 ? EXIT_SUCCESS : EXIT_FAILURE;
}

static
int
test_suzerain_signal_name( void )
{
    int success = 1;

    // The ISO C standard only requires the signal names SIGABRT, SIGFPE,
    // SIGILL, SIGINT, SIGSEGV, and SIGTERM to be defined.

    success |= !strcmp(suzerain_signal_name(SIGABRT), "SIGABRT");
    success |= !strcmp(suzerain_signal_name(SIGFPE),  "SIGFPE");
    success |= !strcmp(suzerain_signal_name(SIGILL),  "SIGILL");
    success |= !strcmp(suzerain_signal_name(SIGINT),  "SIGINT");
    success |= !strcmp(suzerain_signal_name(SIGSEGV), "SIGSEGV");
    success |= !strcmp(suzerain_signal_name(SIGTERM), "SIGTERM");

    return success == 1 ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char *argv[])
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);
    int status = 0;

    status |= test_suzerain_fpipe();
    status |= test_suzerain_signal_number();
    status |= test_suzerain_signal_name();

    if (status) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
