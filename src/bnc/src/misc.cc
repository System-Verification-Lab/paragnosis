#include "misc.h"
#include <bn-to-cnf/bayesnet.h>
#include <iostream>
#include <cstdio>
#include <memory>

std::string exec(const char* cmd) {
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
}

const char* get_bayesnet_filename(bayesnet *bn){
    return bn->get_filename();
}

const char* get_bayesnet_filename(cnf *f){
   return get_bayesnet_filename(f->get_bayesnet());
}

