#ifndef BNC_DEFINES_H
#define BNC_DEFINES_H

//#define MEMORY_SAFE
#define SIMPLE_ALLOCATION
#define VERBOSE
//#define INTERMEDIATE_DOT
//#define WITHOUT_CACHING
//#define TIMED

#if __GNUC__
    #if __x86_64__ || __ppc64__
        #define ENV64BIT
    #else
        #define ENV32BIT
    #endif
#else
    #error "Must define either ENV32BIT or ENV64BIT"
#endif

#endif
