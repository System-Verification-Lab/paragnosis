#ifndef BNMC_INCLUDE_CLI_H_
#define BNMC_INCLUDE_CLI_H_

#include <string>
#include <vector>
#include <unordered_map>

namespace bnmc {

class Cli;

//using Command = int (*)(void*);
typedef int (Cli::*CommandPtr)(void*);
typedef int (Command)(void*);

class Cli {
    public:
        Cli();
        int RunCommand(std::string&);
        void AddCommand(std::string, CommandPtr);
        void Greet();
        void Prompt();
        Command Quit;

    protected:
        static char** TabCompletion(const char*, int ,int);
        static char* TabCompletionGenerator(const char* text, int state);

        std::string arguments_;
        static std::vector< std::string > command_names_;
        std::unordered_map< std::string, CommandPtr> commands_;
};

}

#endif
