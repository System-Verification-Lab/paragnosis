#ifndef BNMC_INCLUDE_IO_H_
#define BNMC_INCLUDE_IO_H_

enum PrintType { MSG = 1, ERR = 2, DBG = 3};


void Log_(const char*, const char*, const int, const PrintType, const char*, ...);
void Debug_(const char*, const char*, const int, const char*, ...);
void Print_(const PrintType, const char*, ...);

#define Debug(...) Debug_(__FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)
#define Log(...)   Log_(__FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)
#define Print(...) Print_( __VA_ARGS__)

#endif

