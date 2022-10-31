#ifndef CUSTOM_STRING_EXCEPTION_H
#define CUSTOM_STRING_EXCEPTION_H

#include "stringf.h"

class string_exception : public std::string, public std::exception {
    public:
        template< typename... argv >
        string_exception( const char* format, argv... args ) {
            std::string::operator=(stringf(format, args...));
        };
        const char* what() const throw() { return c_str(); } ;
};

#define create_exception(T) class T : public string_exception { public: using string_exception::string_exception; };
#define create_derived_exception(T,U) class T : public U { public: using U::U; };

#endif

