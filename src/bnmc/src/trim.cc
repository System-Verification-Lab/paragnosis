#include "trim.h"
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

// trim from start (in place)
inline void ltrim(std::string &s){
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
inline void rtrim(std::string &s){
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
void Trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

char* Trim(char *str){
    int last = 0;
    char *s = str;
    while(s[0] != '\0' && (s[0] == ' ' || s[0] == '\t' || s[0] == '\n')) s++;
    while(s[last] != '\0') last++;
    while(last-1 >= 0 && (s[last-1] == ' ' || s[last-1] == '\t' || s[last-1] == '\n')){
        s[last-1] = '\0';
        last--;
    }

    return s;
}

// void Cli::Trim(std::string &str){
//    // trim white  spaces
//    size_t endpos = str.find_last_not_of(" \t");
//    if( std::string::npos != endpos )
//        str = str.substr( 0, endpos+1 );
//
//    size_t startpos = str.find_first_not_of(" \t");
//    if( std::string::npos != startpos )
//        str = str.substr( startpos );
//}


