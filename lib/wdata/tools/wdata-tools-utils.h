
#define cppmallocl(pointer, size, type)                                     \
    if ((pointer = (type *)malloc((size) * sizeof(type))) == NULL)          \
    {                                                                       \
        fprintf(stderr, "error: cannot malloc()! Exiting!\n");              \
        fprintf(stderr, "error: file=`%s`, line=%d\n", __FILE__, __LINE__); \
        exit(1);                                                          \
    }

#define file_operationl(cmd)                                             \
    {                                                                    \
        ierr = cmd;                                                      \
        if (ierr)                                                        \
        {                                                                \
            fprintf(stderr, "FILE ERROR:: cannot execute: %s\n", #cmd);  \
            fprintf(stderr, "file=`%s`, line=%d\n", __FILE__, __LINE__); \
            fprintf(stderr, "Error=%d\nExiting!\n", ierr);               \
            exit(1);                                       \
        }                                                                \
    }

void copy_txt(const char *src, const char *trg) {
    FILE *file1, *file2;
    char buffer[1024]; // Buffer to hold lines of text
    printf("# COPYING: %s <-- %s\n", trg, src);

    // Open src in append mode
    file1 = fopen(src, "r");
    if (file1 == NULL) {
        printf("# COULD NOT FILE %s FOR COPY. SKIPPIG!\n", src);
        return;
    }

    // Open trg in read mode
    file2 = fopen(trg, "w");
    if (file2 == NULL) {
        printf("# COULD NOT WRITE FILE %s. SKIPPING!\n", trg);
        fclose(file1); // Ensure file1 is closed before exiting
        return;
    }

    // Read content from filename2 and append it to filename1
    while (fgets(buffer, sizeof(buffer), file1) != NULL) {
        fputs(buffer, file2);
    }

    // Close both files
    fclose(file1);
    fclose(file2);

    return;
}
