#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <errno.h>
#include "io.h"


extern char *program_invocation_name;
extern char *program_invocation_short_name;

void Print_(const PrintType kPrintType, const char *format, ...){
    va_list args;
    va_start (args, format);

    if(kPrintType == MSG){
        fprintf(stdout, "  ");
        vfprintf(stdout, format, args);
    } else if(kPrintType == ERR){
        fprintf(stdout, "\033[31m  ");
        vfprintf(stdout, format, args);
        fprintf(stdout, "\033[0m");
    } else {
        fprintf(stdout, "\033[33m  ");
        vfprintf(stdout, format, args);
        fprintf(stdout, "\033[0m");
    }
    fflush(stdout);
    va_end(args);
}

void Debug_(const char *kFilename, const char *kFunction, const int kLine, const char *format, ...){
    va_list args;
    va_start (args, format);

    time_t timer;
    char timestamp[25];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);
    strftime(timestamp, 25, "%b %d %H:%M:%S", tm_info);

    Print_(DBG,"%s:%s:%s:%d: ", timestamp, kFilename, kFunction, kLine);
    Print_(DBG, format, args);

    va_end(args);
}

void Log_(const char *kFilename, const char *kFunction, const int kLine, const PrintType kPrintType, const char *format, ...){
    va_list args;
    va_start (args, format);

    char filename[256];
    if(kPrintType == MSG)
        sprintf(filename, "%s%s", program_invocation_name, ".log");
    else if(kPrintType == ERR)
        sprintf(filename, "%s%s", program_invocation_name, ".err");
    else
        sprintf(filename, "%s%s", program_invocation_name, ".dbg");

    FILE *file = fopen(filename, "a");
    if(file){
        time_t timer;
        char timestamp[25];
        struct tm* tm_info;

        time(&timer);
        tm_info = localtime(&timer);
        strftime(timestamp, 25, "%b %d %H:%M:%S", tm_info);

        fprintf(file,"%s:", timestamp);
        if(kPrintType != MSG)
            fprintf(file, "%s:%s:%d: ", kFilename, kFunction, kLine);

        vfprintf(file, format, args);

        fclose(file);
    }
    #ifdef DEBUG
    else fprintf(stderr, "Could not open log file\n");
    #endif
    va_end(args);
}


