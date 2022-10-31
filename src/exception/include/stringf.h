#ifndef STRINGF_H
#define STRINGF_H

#include <string>

template< typename... argv >
std::string stringf( const char* format, argv... args ) {
    const size_t kSize = std::snprintf( NULL, 0, format, args... );

    std::string output;
    output.resize(kSize);
    std::snprintf( &(output[0]), kSize+1, format, args... );
    return std::move(output);
}

#endif

