#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <suzerain/os.h>

int
test_suzerain_fpipe()
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

int main(int argc, char *argv[])
{
    int status = 0;

    status |= test_suzerain_fpipe();

    if (status) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
