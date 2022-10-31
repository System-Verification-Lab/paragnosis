#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include "ace.h"
#include "options.h"
#include <stdlib.h>
#include <exception/stringf.h>
#include <regex>
#include "exceptions.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

std::string abs_path(std::string &rel_path){
    char *real_path = realpath(rel_path.c_str(), NULL);
    std::string ret_real_path = real_path;
    free(real_path);
    return ret_real_path;
};

std::string exec(const char* cmd) {
    char buffer[128];
    std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try {
        while (!feof(pipe)) {
            if (fgets(buffer, 128, pipe) != NULL)
                result += buffer;
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return result;
}

#ifdef ACE

using namespace bnmc;

probability_t ace_posterior(const Evidence &kEvidence,Timer &t){
    std::string query = kEvidence.GetQueryString();
    std::string bn = manager.files.get_filename(file_t::BN);
    std::string hugin = abs_path(bn);
    std::string cmd = stringf("%s --network %s --query '%s'",TOSTRING(ACE),hugin.c_str(),query.c_str());

    std::string output = exec(cmd.c_str());

    std::regex prob_rgx("Prob : ([0-9.e\-]+)");
    std::regex time_rgx("Time : ([0-9.]+)");
    std::smatch match;

    try {
        probability_t p = -1;
        if (std::regex_search(output, match, prob_rgx)){
            p = std::stod(match.str(1));
            if(p == -1) {
                printf("OUTPUT:\n");
                printf("%s\n",output.c_str());
                printf("===========================\n");
            }
            if (std::regex_search(output, match, time_rgx)){
                double time = std::stod(match.str(1));
                t.Add<Timer::Microseconds>(time*Timer::Microseconds::period::den);
            } else {
                printf("OUTPUT:\n");
                printf("%s\n",output.c_str());
                printf("===========================\n");
                throw ModelCounterException("Could not get Ace timing");
            }
        } else {
            printf("OUTPUT:\n");
            printf("%s\n",output.c_str());
            printf("===========================\n");
            throw ModelCounterException("Could not get Ace probability");
        }
        return p;
    } catch(...) {
        // add average time on failed queries
        t.Add<Timer::Microseconds>(t.GetTotal<Timer::Microseconds>()/bnmc::manager.total_query_count);
        throw;
    }
}

#endif
