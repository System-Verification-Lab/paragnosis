#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "stringsep.h"

using namespace std;

string_sep_t string_sep(string str, char delimiter) {
    string_sep_t internal;
    stringstream ss(str);
    string tok;

    while(getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}


