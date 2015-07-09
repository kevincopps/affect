/* #include <stdio.h> */
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_LEN 256

int main( int argc, char *argv[] ) {
    char buffer[MAX_LEN];
    int err_pipe[2];
    int saved_stderr;

    /* printf("length of buffer string %lu\n", strlen(buffer)); */

    saved_stderr = dup(STDERR_FILENO); /* save stderr for display later */

    if (pipe(err_pipe) != 0)
        exit(1);

    dup2(err_pipe[1], STDERR_FILENO);   /* redirect stderr to the pipe */
    close(err_pipe[1]);

    fprintf(stderr, "writing to pipe and saving for later\n");
    fflush(stderr);
    read(err_pipe[0], buffer, MAX_LEN); /* read from pipe into buffer */

    dup2(saved_stderr, STDERR_FILENO);  /* reconnect stdout for testing */

    fprintf(stderr, "writing to regular stderr\n");

    /* printf("length of buffer string %lu\n", strlen(buffer)); */

    printf("%s\n", buffer);

    return 0;
}