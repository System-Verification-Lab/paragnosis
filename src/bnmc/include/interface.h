#ifndef BNMC_INCLUDE_INTERFACE_H_
#define BNMC_INCLUDE_INTERFACE_H_

#include <unordered_map>
#include <string>
#include <bnc/files.h>
#include "modelcounter.h"
#include "types.h"
#include "cli.h"

namespace bnmc {

class Interface : public Cli {
    public:
        Interface();
        void Init();
        void Prompt();
        enum class ExhaustiveType { COMPARE, VERIFY, ENUMERATE };
    private:
        void AddCommands();
        void Read(file_t);

        void Exhaustive(const ExhaustiveType, const int, const int, const int);
        int Query(const ModelType);
        int Posteriors(const ModelType);

        int Evid();

        probability_t Query(const ModelType, Evidence&, Timer&);

        // cli commands
        Command Help;
        Command Assignments;
        Command Query;
        Command Logfile;
        Command Evid;
        Command Posteriors;
        Command Load;
        Command Compare;

        #ifdef ACE
        Command InitAce;
        #endif

        Command Verify;

        ModelCounter<ModelType::WPBDD> wpbdd_;
        ModelCounter<ModelType::PWPBDD> pwpbdd_;
        ModelCounter<ModelType::MULTIGRAPH> multigraph_;
        ModelCounter<ModelType::TDMULTIGRAPH> tdmultigraph_;
};

}

#endif

